// SPDX-License-Identifier: Apache-2.0
//! Elliptic-curve arithmetic on Montgomery curves over GF(p²).
//!
//! Ported from `src/ec/ref/{include/ec.h, lvlx/*.c}` in the C reference.
//! x-only Kummer-line arithmetic, Jacobian-point arithmetic, 2ⁿ-isogeny
//! evaluation, and deterministic 2ᶠ-torsion basis generation.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use crate::gf::*;
use crate::mp::{self, Digit, LOG2RADIX, RADIX};
use crate::precomp::{NWORDS_ORDER, BITS};

pub mod jac;
pub mod xisog;
pub mod isog_chains;
pub mod basis;
pub mod biextension;

pub use jac::*;
pub use xisog::*;
pub use isog_chains::*;
pub use basis::*;

// ===========================================================================
// Types (ec.h)
// ===========================================================================

/// Projective point on the Kummer line E/±1 in Montgomery (X:Z) coordinates.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcPoint {
    pub x: Fp2,
    pub z: Fp2,
}

/// Projective point in Jacobian (X:Y:Z) coordinates representing (X/Z², Y/Z³).
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct JacPoint {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
}

/// Addition components: (P+Q) = (u-v : w), (P-Q) = (u+v : w).
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct AddComponents {
    pub u: Fp2,
    pub v: Fp2,
    pub w: Fp2,
}

/// A torsion-subgroup basis (P, Q, P-Q) in x-only form.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcBasis {
    pub p: EcPoint,
    pub q: EcPoint,
    pub pmq: EcPoint,
}

/// A Montgomery curve y² = x³ + (A/C)x² + x in projective form (A : C),
/// optionally caching the normalized A24 = ((A+2C)/4C : 1).
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcCurve {
    pub a: Fp2,
    pub c: Fp2,
    pub a24: EcPoint,
    pub is_a24_computed_and_normalized: bool,
}

/// A 2ⁿ-isogeny: domain curve, kernel generator, and walk length.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcIsogEven {
    pub curve: EcCurve,
    pub kernel: EcPoint,
    pub length: u32,
}

/// Isomorphism (X:Z) ↦ (Nx·X + Nz·Z : D·Z) between Montgomery curves.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcIsom {
    pub nx: Fp2,
    pub nz: Fp2,
    pub d: Fp2,
}

/// KPS state for a degree-2 isogeny.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcKps2 {
    pub k: EcPoint,
}

/// KPS state for a degree-4 isogeny.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct EcKps4 {
    pub k: [EcPoint; 3],
}

// ===========================================================================
// Inline helpers from ec.h
// ===========================================================================

#[inline]
pub fn copy_point(p: &mut EcPoint, q: &EcPoint) {
    *p = *q;
}

#[inline]
pub fn copy_basis(b1: &mut EcBasis, b0: &EcBasis) {
    *b1 = *b0;
}

#[inline]
pub fn copy_curve(e1: &mut EcCurve, e2: &EcCurve) {
    *e1 = *e2;
}

/// Compute A24 = (A+2C : 4C), or copy the cached normalized value.
#[inline]
pub fn ac_to_a24(a24: &mut EcPoint, e: &EcCurve) {
    if e.is_a24_computed_and_normalized {
        copy_point(a24, &e.a24);
        return;
    }
    fp2_add(&mut a24.z, &e.c, &e.c);
    let z = a24.z;
    fp2_add(&mut a24.x, &e.a, &z);
    fp2_add(&mut a24.z, &z, &z);
}

/// Recover (A : C) = (2·(A24.x·2 - A24.z) : A24.z) from A24 = (A+2C : 4C).
#[inline]
pub fn a24_to_ac(e: &mut EcCurve, a24: &EcPoint) {
    let mut t = Fp2::default();
    fp2_add(&mut t, &a24.x, &a24.x);
    fp2_sub(&mut e.a, &t, &a24.z);
    let a = e.a;
    fp2_add(&mut e.a, &a, &a);
    fp2_copy(&mut e.c, &a24.z);
}

// ===========================================================================
// ec.c port
// ===========================================================================

/// Initialize a point as the identity (1 : 0).
pub fn ec_point_init(p: &mut EcPoint) {
    fp2_set_one(&mut p.x);
    fp2_set_zero(&mut p.z);
}

/// Initialize a curve as E₀ with (A : C) = (0 : 1).
pub fn ec_curve_init(e: &mut EcCurve) {
    fp2_set_zero(&mut e.a);
    fp2_set_one(&mut e.c);
    ec_point_init(&mut e.a24);
    e.is_a24_computed_and_normalized = false;
}

/// Constant-time select: Q ← P1 if option==0, else P2 (option must be 0 or all-ones).
#[inline]
pub fn select_point(q: &mut EcPoint, p1: &EcPoint, p2: &EcPoint, option: Digit) {
    let ctl = option as u32;
    fp2_select(&mut q.x, &p1.x, &p2.x, ctl);
    fp2_select(&mut q.z, &p1.z, &p2.z, ctl);
}

/// Constant-time conditional swap: if option is all-ones swap P↔Q.
#[inline]
pub fn cswap_points(p: &mut EcPoint, q: &mut EcPoint, option: Digit) {
    let ctl = option as u32;
    fp2_cswap(&mut p.x, &mut q.x, ctl);
    fp2_cswap(&mut p.z, &mut q.z, ctl);
}

/// Normalize (X : Z) → (X/Z : 1) in place.
pub fn ec_normalize_point(p: &mut EcPoint) {
    fp2_inv(&mut p.z);
    let (x, z) = (p.x, p.z);
    fp2_mul(&mut p.x, &x, &z);
    fp2_set_one(&mut p.z);
}

/// Normalize (A : C) → (A/C : 1) in place.
pub fn ec_normalize_curve(e: &mut EcCurve) {
    fp2_inv(&mut e.c);
    let (a, c) = (e.a, e.c);
    fp2_mul(&mut e.a, &a, &c);
    fp2_set_one(&mut e.c);
}

/// Compute and cache the normalized A24 = ((A+2C)/4C : 1).
pub fn ec_curve_normalize_a24(e: &mut EcCurve) {
    if !e.is_a24_computed_and_normalized {
        let mut a24 = EcPoint::default();
        ac_to_a24(&mut a24, e);
        ec_normalize_point(&mut a24);
        e.a24 = a24;
        e.is_a24_computed_and_normalized = true;
    }
    debug_assert!(fp2_is_one(&e.a24.z) != 0);
}

/// Normalize both (A : C) and A24 in place.
pub fn ec_normalize_curve_and_a24(e: &mut EcCurve) {
    if fp2_is_one(&e.c) == 0 {
        ec_normalize_curve(e);
    }
    if !e.is_a24_computed_and_normalized {
        let a = e.a;
        fp2_add_one(&mut e.a24.x, &a);
        let t = e.a24.x;
        fp2_add_one(&mut e.a24.x, &t);
        fp_copy(&mut e.a24.x.im, &e.a.im);
        let t = e.a24.x;
        fp2_half(&mut e.a24.x, &t);
        let t = e.a24.x;
        fp2_half(&mut e.a24.x, &t);
        fp2_set_one(&mut e.a24.z);
        e.is_a24_computed_and_normalized = true;
    }
}

#[inline]
pub fn ec_is_zero(p: &EcPoint) -> u32 {
    fp2_is_zero(&p.z)
}

#[inline]
pub fn ec_has_zero_coordinate(p: &EcPoint) -> u32 {
    fp2_is_zero(&p.x) | fp2_is_zero(&p.z)
}

/// Projective equality on (X:Z): P = Q iff both zero, or neither zero and Px·Qz = Qx·Pz.
pub fn ec_is_equal(p: &EcPoint, q: &EcPoint) -> u32 {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let l_zero = ec_is_zero(p);
    let r_zero = ec_is_zero(q);
    fp2_mul(&mut t0, &p.x, &q.z);
    fp2_mul(&mut t1, &p.z, &q.x);
    let lr_equal = fp2_is_equal(&t0, &t1);
    // Faithful port of C's `(~l_zero & ~r_zero * lr_equal)`: `*` binds tighter than `&`.
    (l_zero & r_zero) | (!l_zero & (!r_zero).wrapping_mul(lr_equal))
}

/// Nonzero iff P is exact 2-torsion (nonzero and 2P = 0).
pub fn ec_is_two_torsion(p: &EcPoint, e: &EcCurve) -> u32 {
    if ec_is_zero(p) != 0 {
        return 0;
    }
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    let s = t0;
    fp2_sqr(&mut t0, &s);
    fp2_sub(&mut t1, &p.x, &p.z);
    let s = t1;
    fp2_sqr(&mut t1, &s);
    fp2_sub(&mut t2, &t0, &t1);
    let s = t1;
    fp2_add(&mut t1, &t0, &s);
    let s = t2;
    fp2_mul(&mut t2, &s, &e.a);
    let s = t1;
    fp2_mul(&mut t1, &s, &e.c);
    let s = t1;
    fp2_add(&mut t1, &s, &s);
    fp2_add(&mut t0, &t1, &t2);
    let x_is_zero = fp2_is_zero(&p.x);
    let tmp_is_zero = fp2_is_zero(&t0);
    x_is_zero | tmp_is_zero
}

/// Nonzero iff P is exact 4-torsion.
pub fn ec_is_four_torsion(p: &EcPoint, e: &EcCurve) -> u32 {
    let mut test = EcPoint::default();
    xdbl_a24(&mut test, p, &e.a24, e.is_a24_computed_and_normalized);
    ec_is_two_torsion(&test, e)
}

/// Nonzero iff (P, Q) form a full 4-basis.
pub fn ec_is_basis_four_torsion(b: &EcBasis, e: &EcCurve) -> u32 {
    let mut p2 = EcPoint::default();
    let mut q2 = EcPoint::default();
    xdbl_a24(&mut p2, &b.p, &e.a24, e.is_a24_computed_and_normalized);
    xdbl_a24(&mut q2, &b.q, &e.a24, e.is_a24_computed_and_normalized);
    ec_is_two_torsion(&p2, e) & ec_is_two_torsion(&q2, e) & !ec_is_equal(&p2, &q2)
}

/// 1 if A is a valid Montgomery coefficient (A ≠ ±2), else 0.
pub fn ec_curve_verify_a(a: &Fp2) -> i32 {
    let mut t = Fp2::default();
    fp2_set_one(&mut t);
    let one = t.re;
    fp_add(&mut t.re, &one, &one);
    if fp2_is_equal(a, &t) != 0 {
        return 0;
    }
    let two = t.re;
    fp_neg(&mut t.re, &two);
    if fp2_is_equal(a, &t) != 0 {
        return 0;
    }
    1
}

/// Initialize a curve from coefficient A; returns 1 if valid, else 0.
pub fn ec_curve_init_from_a(e: &mut EcCurve, a: &Fp2) -> i32 {
    ec_curve_init(e);
    fp2_copy(&mut e.a, a);
    ec_curve_verify_a(a)
}

/// j-invariant of the Montgomery curve.
pub fn ec_j_inv(j_inv: &mut Fp2, curve: &EcCurve) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    fp2_sqr(&mut t1, &curve.c);
    fp2_sqr(j_inv, &curve.a);
    fp2_add(&mut t0, &t1, &t1);
    let j = *j_inv;
    let s = t0;
    fp2_sub(&mut t0, &j, &s);
    let s = t0;
    fp2_sub(&mut t0, &s, &t1);
    fp2_sub(j_inv, &t0, &t1);
    let s = t1;
    fp2_sqr(&mut t1, &s);
    let j = *j_inv;
    fp2_mul(j_inv, &j, &t1);
    let s = t0;
    fp2_add(&mut t0, &s, &s);
    let s = t0;
    fp2_add(&mut t0, &s, &s);
    fp2_sqr(&mut t1, &t0);
    let s = t0;
    fp2_mul(&mut t0, &s, &t1);
    let s = t0;
    fp2_add(&mut t0, &s, &s);
    let s = t0;
    fp2_add(&mut t0, &s, &s);
    fp2_inv(j_inv);
    let j = *j_inv;
    fp2_mul(j_inv, &t0, &j);
}

/// x-only doubling on E₀ (A=0, C=1).
pub fn xdbl_e0(q: &mut EcPoint, p: &EcPoint) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    let s = t0;
    fp2_sqr(&mut t0, &s);
    fp2_sub(&mut t1, &p.x, &p.z);
    let s = t1;
    fp2_sqr(&mut t1, &s);
    fp2_sub(&mut t2, &t0, &t1);
    let s = t1;
    fp2_add(&mut t1, &s, &s);
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_add(&mut q.z, &t1, &t2);
    let s = q.z;
    fp2_mul(&mut q.z, &s, &t2);
}

/// x-only doubling, computing A24 from (A:C) on the fly.
pub fn xdbl(q: &mut EcPoint, p: &EcPoint, ac: &EcPoint) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    let s = t0;
    fp2_sqr(&mut t0, &s);
    fp2_sub(&mut t1, &p.x, &p.z);
    let s = t1;
    fp2_sqr(&mut t1, &s);
    fp2_sub(&mut t2, &t0, &t1);
    fp2_add(&mut t3, &ac.z, &ac.z);
    let s = t1;
    fp2_mul(&mut t1, &s, &t3);
    let s = t1;
    fp2_add(&mut t1, &s, &s);
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_add(&mut t0, &t3, &ac.x);
    let s = t0;
    fp2_mul(&mut t0, &s, &t2);
    let s = t0;
    fp2_add(&mut t0, &s, &t1);
    fp2_mul(&mut q.z, &t0, &t2);
}

/// x-only doubling using precomputed A24 = (A+2C : 4C) (or normalized).
pub fn xdbl_a24(q: &mut EcPoint, p: &EcPoint, a24: &EcPoint, a24_normalized: bool) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    let s = t0;
    fp2_sqr(&mut t0, &s);
    fp2_sub(&mut t1, &p.x, &p.z);
    let s = t1;
    fp2_sqr(&mut t1, &s);
    fp2_sub(&mut t2, &t0, &t1);
    if !a24_normalized {
        let s = t1;
        fp2_mul(&mut t1, &s, &a24.z);
    }
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_mul(&mut t0, &t2, &a24.x);
    let s = t0;
    fp2_add(&mut t0, &s, &t1);
    fp2_mul(&mut q.z, &t0, &t2);
}

/// Differential addition: R = P + Q given PQ = P - Q.
pub fn xadd(r: &mut EcPoint, p: &EcPoint, q: &EcPoint, pq: &EcPoint) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut t2, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    let s = t0;
    fp2_mul(&mut t0, &s, &t3);
    let s = t1;
    fp2_mul(&mut t1, &s, &t2);
    fp2_add(&mut t2, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    let s = t2;
    fp2_sqr(&mut t2, &s);
    let s = t3;
    fp2_sqr(&mut t3, &s);
    let s = t2;
    fp2_mul(&mut t2, &pq.z, &s);
    fp2_mul(&mut r.z, &pq.x, &t3);
    fp2_copy(&mut r.x, &t2);
}

/// Simultaneous doubling and differential addition: R = 2P, S = P+Q.
pub fn xdbladd(
    r: &mut EcPoint,
    s: &mut EcPoint,
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    a24: &EcPoint,
    a24_normalized: bool,
) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_sqr(&mut r.x, &t0);
    fp2_sub(&mut t2, &q.x, &q.z);
    fp2_add(&mut s.x, &q.x, &q.z);
    let st0 = t0;
    fp2_mul(&mut t0, &st0, &t2);
    fp2_sqr(&mut r.z, &t1);
    let st1 = t1;
    fp2_mul(&mut t1, &st1, &s.x);
    let rx = r.x;
    fp2_sub(&mut t2, &rx, &r.z);
    if !a24_normalized {
        let rz = r.z;
        fp2_mul(&mut r.z, &rz, &a24.z);
    }
    let (rx, rz) = (r.x, r.z);
    fp2_mul(&mut r.x, &rx, &rz);
    fp2_mul(&mut s.x, &a24.x, &t2);
    fp2_sub(&mut s.z, &t0, &t1);
    let (rz, sx) = (r.z, s.x);
    fp2_add(&mut r.z, &rz, &sx);
    fp2_add(&mut s.x, &t0, &t1);
    let rz = r.z;
    fp2_mul(&mut r.z, &rz, &t2);
    let sz = s.z;
    fp2_sqr(&mut s.z, &sz);
    let sx = s.x;
    fp2_sqr(&mut s.x, &sx);
    let sz = s.z;
    fp2_mul(&mut s.z, &sz, &pq.x);
    let sx = s.x;
    fp2_mul(&mut s.x, &sx, &pq.z);
}

/// Montgomery ladder: Q = k·P over `kbits` bits.
pub fn xmul(q: &mut EcPoint, p: &EcPoint, k: &[Digit], kbits: i32, curve: &EcCurve) {
    let mut a24 = EcPoint::default();
    if !curve.is_a24_computed_and_normalized {
        fp2_add(&mut a24.x, &curve.c, &curve.c);
        let s = a24.x;
        fp2_add(&mut a24.z, &s, &s);
        let s = a24.x;
        fp2_add(&mut a24.x, &s, &curve.a);
    } else {
        fp2_copy(&mut a24.x, &curve.a24.x);
        fp2_copy(&mut a24.z, &curve.a24.z);
        debug_assert!(fp2_is_one(&a24.z) != 0);
    }

    let mut r0 = EcPoint::default();
    let mut r1 = EcPoint::default();
    ec_point_init(&mut r0);
    fp2_copy(&mut r1.x, &p.x);
    fp2_copy(&mut r1.z, &p.z);

    let mut prevbit: u32 = 0;
    for i in (0..kbits).rev() {
        let bit = ((k[(i as u32 >> LOG2RADIX) as usize] >> (i as u32 & (RADIX - 1))) & 1) as u32;
        let swap = bit ^ prevbit;
        prevbit = bit;
        let mask = 0u64.wrapping_sub(swap as u64);
        cswap_points(&mut r0, &mut r1, mask);
        let (in0, in1) = (r0, r1);
        xdbladd(&mut r0, &mut r1, &in0, &in1, p, &a24, true);
    }
    let mask = 0u64.wrapping_sub(prevbit as u64);
    cswap_points(&mut r0, &mut r1, mask);

    fp2_copy(&mut q.x, &r0.x);
    fp2_copy(&mut q.z, &r0.z);
}

/// Biscalar Montgomery ladder: S = k·P + l·Q. Returns 1 on success, 0 on
/// invalid inputs (zero coordinates, etc).
pub fn xdblmul(
    s_out: &mut EcPoint,
    p: &EcPoint,
    k: &[Digit],
    q: &EcPoint,
    l: &[Digit],
    pq: &EcPoint,
    kbits: i32,
    curve: &EcCurve,
) -> i32 {
    if (ec_has_zero_coordinate(p) | ec_has_zero_coordinate(q) | ec_has_zero_coordinate(pq)) != 0 {
        return 0;
    }

    let bitk0 = k[0] & 1;
    let bitl0 = l[0] & 1;
    let maskk = 0u64.wrapping_sub(bitk0);
    let maskl = 0u64.wrapping_sub(bitl0);
    let mut sigma = [bitk0 ^ 1, bitl0 ^ 1];
    let evens = sigma[0] + sigma[1];
    let mevens = 0u64.wrapping_sub(evens & 1);

    sigma[0] &= mevens;
    sigma[1] = (sigma[1] & mevens) | (1 & !mevens);

    let mut one = [0 as Digit; NWORDS_ORDER];
    one[0] = 1;
    let mut k_t = [0 as Digit; NWORDS_ORDER];
    let mut l_t = [0 as Digit; NWORDS_ORDER];
    mp::mp_sub(&mut k_t, k, &one, NWORDS_ORDER);
    mp::mp_sub(&mut l_t, l, &one, NWORDS_ORDER);
    let kt = k_t;
    mp::select_ct(&mut k_t, &kt, k, maskk, NWORDS_ORDER);
    let lt = l_t;
    mp::select_ct(&mut l_t, &lt, l, maskl, NWORDS_ORDER);

    let mut r = [0 as Digit; 2 * BITS];
    let mut pre_sigma: Digit = 0;
    for i in 0..kbits as usize {
        let maskk = 0u64.wrapping_sub(sigma[0] ^ pre_sigma);
        mp::swap_ct(&mut k_t, &mut l_t, maskk, NWORDS_ORDER);

        let (bs1_ip1, bs2_ip1) = if i == kbits as usize - 1 {
            (0, 0)
        } else {
            (
                mp::mp_shiftr(&mut k_t, 1, NWORDS_ORDER),
                mp::mp_shiftr(&mut l_t, 1, NWORDS_ORDER),
            )
        };
        let bs1_i = k_t[0] & 1;
        let bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        pre_sigma = sigma[0];
        let maskk = 0u64.wrapping_sub(r[2 * i + 1]);
        let mut temp = [0 as Digit; 1];
        let s0 = [sigma[0]];
        let s1 = [sigma[1]];
        mp::select_ct(&mut temp, &s0, &s1, maskk, 1);
        let mut s1_out = [0 as Digit; 1];
        mp::select_ct(&mut s1_out, &s1, &s0, maskk, 1);
        sigma[1] = s1_out[0];
        sigma[0] = temp[0];
    }

    let mut rr = [EcPoint::default(); 3];
    ec_point_init(&mut rr[0]);
    let maskk = 0u64.wrapping_sub(sigma[0]);
    select_point(&mut rr[1], p, q, maskk);
    select_point(&mut rr[2], q, p, maskk);

    let mut diff1a = EcPoint::default();
    let mut diff1b = EcPoint::default();
    let mut diff2a = EcPoint::default();
    let mut diff2b = EcPoint::default();
    fp2_copy(&mut diff1a.x, &rr[1].x);
    fp2_copy(&mut diff1a.z, &rr[1].z);
    fp2_copy(&mut diff1b.x, &rr[2].x);
    fp2_copy(&mut diff1b.z, &rr[2].z);

    let (r1, r2) = (rr[1], rr[2]);
    xadd(&mut rr[2], &r1, &r2, pq);
    if ec_has_zero_coordinate(&rr[2]) != 0 {
        return 0;
    }
    fp2_copy(&mut diff2a.x, &rr[2].x);
    fp2_copy(&mut diff2a.z, &rr[2].z);
    fp2_copy(&mut diff2b.x, &pq.x);
    fp2_copy(&mut diff2b.z, &pq.z);

    let a_is_zero = fp2_is_zero(&curve.a) != 0;

    let mut t = [EcPoint::default(); 3];
    for i in (0..kbits as usize).rev() {
        let h = r[2 * i] + r[2 * i + 1];
        let maskk = 0u64.wrapping_sub(h & 1);
        select_point(&mut t[0], &rr[0], &rr[1], maskk);
        let maskk = 0u64.wrapping_sub(h >> 1);
        let t0 = t[0];
        select_point(&mut t[0], &t0, &rr[2], maskk);
        let t0 = t[0];
        if a_is_zero {
            xdbl_e0(&mut t[0], &t0);
        } else {
            debug_assert!(fp2_is_one(&curve.a24.z) != 0);
            xdbl_a24(&mut t[0], &t0, &curve.a24, true);
        }

        let maskk = 0u64.wrapping_sub(r[2 * i + 1]);
        select_point(&mut t[1], &rr[0], &rr[1], maskk);
        select_point(&mut t[2], &rr[1], &rr[2], maskk);

        cswap_points(&mut diff1a, &mut diff1b, maskk);
        let (t1, t2) = (t[1], t[2]);
        xadd(&mut t[1], &t1, &t2, &diff1a);
        xadd(&mut t[2], &rr[0], &rr[2], &diff2a);

        let maskk = 0u64.wrapping_sub(h & 1);
        cswap_points(&mut diff2a, &mut diff2b, maskk);

        rr[0] = t[0];
        rr[1] = t[1];
        rr[2] = t[2];
    }

    select_point(s_out, &rr[0], &rr[1], mevens);
    let maskk = 0u64.wrapping_sub(bitk0 & bitl0);
    let s0 = *s_out;
    select_point(s_out, &s0, &rr[2], maskk);
    1
}

/// 3-point ladder: R = P + m·Q. Returns 1 on success.
pub fn ec_ladder3pt(
    r: &mut EcPoint,
    m: &[Digit],
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    e: &EcCurve,
) -> i32 {
    debug_assert!(e.is_a24_computed_and_normalized);
    if fp2_is_one(&e.a24.z) == 0 {
        return 0;
    }
    if ec_has_zero_coordinate(pq) != 0 {
        return 0;
    }

    let mut x0 = *q;
    let mut x1 = *p;
    let mut x2 = *pq;

    for i in 0..NWORDS_ORDER {
        let mut t: Digit = 1;
        for _ in 0..RADIX {
            let mask = ((t & m[i] == 0) as u64).wrapping_neg();
            cswap_points(&mut x1, &mut x2, mask);
            let (in0, in1) = (x0, x1);
            xdbladd(&mut x0, &mut x1, &in0, &in1, &x2, &e.a24, true);
            cswap_points(&mut x1, &mut x2, mask);
            t <<= 1;
        }
    }
    *r = x1;
    1
}

/// Single doubling, choosing the cheapest formula given curve normalization.
pub fn ec_dbl(res: &mut EcPoint, p: &EcPoint, curve: &EcCurve) {
    if curve.is_a24_computed_and_normalized {
        debug_assert!(fp2_is_one(&curve.a24.z) != 0);
        xdbl_a24(res, p, &curve.a24, true);
    } else {
        let ac = EcPoint { x: curve.a, z: curve.c };
        xdbl(res, p, &ac);
    }
}

/// n-fold doubling.
pub fn ec_dbl_iter(res: &mut EcPoint, n: i32, p: &EcPoint, curve: &mut EcCurve) {
    if n == 0 {
        copy_point(res, p);
        return;
    }
    if n > 50 {
        ec_curve_normalize_a24(curve);
    }
    if curve.is_a24_computed_and_normalized {
        debug_assert!(fp2_is_one(&curve.a24.z) != 0);
        xdbl_a24(res, p, &curve.a24, true);
        for _ in 0..n - 1 {
            let s = *res;
            xdbl_a24(res, &s, &curve.a24, true);
        }
    } else {
        let ac = EcPoint { x: curve.a, z: curve.c };
        xdbl(res, p, &ac);
        for _ in 0..n - 1 {
            let s = *res;
            xdbl(res, &s, &ac);
        }
    }
}

/// n-fold doubling of an entire basis.
pub fn ec_dbl_iter_basis(res: &mut EcBasis, n: i32, b: &EcBasis, curve: &mut EcCurve) {
    ec_dbl_iter(&mut res.p, n, &b.p, curve);
    ec_dbl_iter(&mut res.q, n, &b.q, curve);
    ec_dbl_iter(&mut res.pmq, n, &b.pmq, curve);
}

/// Scalar multiplication wrapper.
pub fn ec_mul(res: &mut EcPoint, scalar: &[Digit], kbits: i32, p: &EcPoint, curve: &mut EcCurve) {
    if kbits > 50 {
        ec_curve_normalize_a24(curve);
    }
    xmul(res, p, scalar, kbits, curve);
}

/// Biscalar multiplication: res = scalarP·P + scalarQ·Q. Returns 1 on success.
pub fn ec_biscalar_mul(
    res: &mut EcPoint,
    scalar_p: &[Digit],
    scalar_q: &[Digit],
    kbits: i32,
    pq: &EcBasis,
    curve: &EcCurve,
) -> i32 {
    if fp2_is_zero(&pq.pmq.z) != 0 {
        return 0;
    }
    if kbits == 1 {
        if ec_is_two_torsion(&pq.p, curve) == 0
            || ec_is_two_torsion(&pq.q, curve) == 0
            || ec_is_two_torsion(&pq.pmq, curve) == 0
        {
            return 0;
        }
        let bp = scalar_p[0] & 1;
        let bq = scalar_q[0] & 1;
        match (bp, bq) {
            (0, 0) => ec_point_init(res),
            (1, 0) => copy_point(res, &pq.p),
            (0, 1) => copy_point(res, &pq.q),
            (1, 1) => copy_point(res, &pq.pmq),
            _ => unreachable!(),
        }
        return 1;
    }
    let mut e = *curve;
    if fp2_is_zero(&curve.a) == 0 {
        ec_curve_normalize_a24(&mut e);
    }
    xdblmul(res, &pq.p, scalar_p, &pq.q, scalar_q, &pq.pmq, kbits, &e)
}

// ===========================================================================
// Debug / test helpers (from ec.h static inlines)
// ===========================================================================

/// Nonzero iff P has order exactly 2^t.
pub fn test_point_order_twof(p: &EcPoint, e: &EcCurve, t: i32) -> u32 {
    let mut test = *p;
    let mut curve = *e;
    if ec_is_zero(&test) != 0 {
        return 0;
    }
    let pt = test;
    ec_dbl_iter(&mut test, t - 1, &pt, &mut curve);
    if ec_is_zero(&test) != 0 {
        return 0;
    }
    let pt = test;
    ec_dbl(&mut test, &pt, &curve);
    ec_is_zero(&test)
}

/// Nonzero iff all of P, Q, PmQ have order exactly 2^t.
pub fn test_basis_order_twof(b: &EcBasis, e: &EcCurve, t: i32) -> u32 {
    test_point_order_twof(&b.p, e, t)
        & test_point_order_twof(&b.q, e, t)
        & test_point_order_twof(&b.pmq, e, t)
}

/// Nonzero iff Jacobian P has order exactly 2^t.
pub fn test_jac_order_twof(p: &JacPoint, e: &EcCurve, t: i32) -> u32 {
    let mut test = *p;
    if fp2_is_zero(&test.z) != 0 {
        return 0;
    }
    for _ in 0..t - 1 {
        let s = test;
        jac_dbl(&mut test, &s, e);
    }
    if fp2_is_zero(&test.z) != 0 {
        return 0;
    }
    let s = test;
    jac_dbl(&mut test, &s, e);
    fp2_is_zero(&test.z)
}

#[allow(dead_code)]
pub fn ec_point_print(name: &str, p: &EcPoint) {
    if fp2_is_zero(&p.z) != 0 {
        println!("{} = INF", name);
    } else {
        let mut a = p.z;
        fp2_inv(&mut a);
        let s = a;
        fp2_mul(&mut a, &s, &p.x);
        fp2_print(name, &a);
    }
}

#[allow(dead_code)]
pub fn ec_curve_print(name: &str, e: &EcCurve) {
    let mut a = e.c;
    fp2_inv(&mut a);
    let s = a;
    fp2_mul(&mut a, &s, &e.a);
    fp2_print(name, &a);
}

// ===========================================================================
// Tests (port of curve-arith-test.c, basis-gen-test.c)
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precomp::TORSION_EVEN_POWER;

    // Test-only PRNG matching the pattern in gf::tests.
    struct Prng(u64);
    impl Prng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn fp(&mut self) -> Fp {
            let mut buf = [0u8; FP_ENCODED_BYTES];
            for chunk in buf.chunks_mut(8) {
                let x = self.next().to_le_bytes();
                chunk.copy_from_slice(&x[..chunk.len()]);
            }
            let mut a = Fp::default();
            fp_decode_reduce(&mut a, &buf);
            a
        }
        fn fp2(&mut self) -> Fp2 {
            Fp2 { re: self.fp(), im: self.fp() }
        }
    }

    fn projective_is_on_curve(p: &EcPoint, curve: &EcCurve) -> bool {
        let mut t0 = Fp2::default();
        let mut t1 = Fp2::default();
        let mut t2 = Fp2::default();
        fp2_mul(&mut t0, &curve.c, &p.x);
        fp2_mul(&mut t1, &t0, &p.z);
        let s = t1;
        fp2_mul(&mut t1, &s, &curve.a);
        fp2_mul(&mut t2, &curve.c, &p.z);
        let s = t0;
        fp2_sqr(&mut t0, &s);
        let s = t2;
        fp2_sqr(&mut t2, &s);
        let s = t0;
        fp2_add(&mut t0, &s, &t1);
        let s = t0;
        fp2_add(&mut t0, &s, &t2);
        let s = t0;
        fp2_mul(&mut t0, &s, &p.x);
        let s = t0;
        fp2_mul(&mut t0, &s, &p.z);
        fp2_is_square(&t0) != 0 || fp2_is_zero(&t0) != 0
    }

    fn ec_random(prng: &mut Prng, curve: &EcCurve) -> EcPoint {
        let mut p = EcPoint::default();
        fp2_set_one(&mut p.z);
        loop {
            p.x = prng.fp2();
            if projective_is_on_curve(&p, curve) {
                break;
            }
        }
        let z = prng.fp2();
        let px = p.x;
        fp2_mul(&mut p.x, &px, &z);
        p.z = z;
        p
    }

    /// Port of test_extras.c::projective_difference_point.
    fn projective_difference_point(pq: &mut EcPoint, p: &EcPoint, q: &EcPoint, curve: &EcCurve) {
        let mut bxx = Fp2::default();
        let mut bxz = Fp2::default();
        let mut bzz = Fp2::default();
        let mut t0 = Fp2::default();
        let mut t1 = Fp2::default();
        fp2_mul(&mut t0, &p.x, &q.x);
        fp2_mul(&mut t1, &p.z, &q.z);
        fp2_sub(&mut bxx, &t0, &t1);
        let s = bxx;
        fp2_sqr(&mut bxx, &s);
        let s = bxx;
        fp2_mul(&mut bxx, &s, &curve.c);
        fp2_add(&mut bxz, &t0, &t1);
        fp2_mul(&mut t0, &p.x, &q.z);
        fp2_mul(&mut t1, &p.z, &q.x);
        fp2_add(&mut bzz, &t0, &t1);
        let s = bxz;
        fp2_mul(&mut bxz, &s, &bzz);
        fp2_sub(&mut bzz, &t0, &t1);
        let s = bzz;
        fp2_sqr(&mut bzz, &s);
        let s = bzz;
        fp2_mul(&mut bzz, &s, &curve.c);
        let s = bxz;
        fp2_mul(&mut bxz, &s, &curve.c);
        let s = t0;
        fp2_mul(&mut t0, &s, &t1);
        let s = t0;
        fp2_mul(&mut t0, &s, &curve.a);
        let s = t0;
        fp2_add(&mut t0, &s, &s);
        let s = bxz;
        fp2_add(&mut bxz, &s, &t0);

        fp_copy(&mut t0.re, &curve.c.re);
        fp_neg(&mut t0.im, &curve.c.im);
        let s = t0;
        fp2_sqr(&mut t0, &s);
        let s = t0;
        fp2_mul(&mut t0, &s, &curve.c);
        fp_copy(&mut t1.re, &p.z.re);
        fp_neg(&mut t1.im, &p.z.im);
        let s = t1;
        fp2_sqr(&mut t1, &s);
        let s = t0;
        fp2_mul(&mut t0, &s, &t1);
        fp_copy(&mut t1.re, &q.z.re);
        fp_neg(&mut t1.im, &q.z.im);
        let s = t1;
        fp2_sqr(&mut t1, &s);
        let s = t0;
        fp2_mul(&mut t0, &s, &t1);
        let s = bxx;
        fp2_mul(&mut bxx, &s, &t0);
        let s = bxz;
        fp2_mul(&mut bxz, &s, &t0);
        let s = bzz;
        fp2_mul(&mut bzz, &s, &t0);

        fp2_sqr(&mut t0, &bxz);
        fp2_mul(&mut t1, &bxx, &bzz);
        let s = t0;
        fp2_sub(&mut t0, &s, &t1);
        fp2_sqrt(&mut t0);
        fp2_add(&mut pq.x, &bxz, &t0);
        fp2_copy(&mut pq.z, &bzz);
    }

    fn make_e0() -> EcCurve {
        let mut curve = EcCurve::default();
        ec_curve_init(&mut curve);
        ec_curve_normalize_a24(&mut curve);
        curve
    }

    const ITERS: usize = 50;

    #[test]
    fn j_invariant_e0_is_1728() {
        let curve = make_e0();
        let mut j = Fp2::default();
        ec_j_inv(&mut j, &curve);
        let mut expected = Fp2::default();
        fp2_set_small(&mut expected, 1728);
        assert_ne!(fp2_is_equal(&j, &expected), 0, "j(E0) != 1728");
    }

    #[test]
    fn xdbl_xadd() {
        let curve = make_e0();
        let mut prng = Prng(0xC0FFEE);
        for _ in 0..ITERS {
            let mut p = ec_random(&mut prng, &curve);
            let mut q = ec_random(&mut prng, &curve);
            let mut pq = EcPoint::default();
            projective_difference_point(&mut pq, &p, &q, &curve);

            // 2(P+Q) = 2P + 2Q
            let mut r1 = EcPoint::default();
            let mut r2 = EcPoint::default();
            xadd(&mut r1, &p, &q, &pq);
            let s = r1;
            ec_dbl(&mut r1, &s, &curve);
            let (sp, sq, spq) = (p, q, pq);
            ec_dbl(&mut p, &sp, &curve);
            ec_dbl(&mut q, &sq, &curve);
            ec_dbl(&mut pq, &spq, &curve);
            xadd(&mut r2, &p, &q, &pq);
            assert_ne!(ec_is_equal(&r1, &r2), 0, "2(P+Q) != 2P+2Q");

            // (P+Q)+(P-Q) = 2P
            xadd(&mut r1, &p, &q, &pq);
            let sq = q;
            ec_dbl(&mut q, &sq, &curve);
            let sr1 = r1;
            xadd(&mut r1, &sr1, &pq, &q);
            let sp = p;
            ec_dbl(&mut p, &sp, &curve);
            assert_ne!(ec_is_equal(&r1, &p), 0, "(P+Q)+(P-Q) != 2P");
        }
    }

    #[test]
    fn xdbladd_matches_separate() {
        let curve = make_e0();
        let mut a24 = EcPoint::default();
        ac_to_a24(&mut a24, &curve);
        let mut prng = Prng(0xBEEF);
        for _ in 0..ITERS {
            let p = ec_random(&mut prng, &curve);
            let q = ec_random(&mut prng, &curve);
            let mut pq = EcPoint::default();
            projective_difference_point(&mut pq, &p, &q, &curve);

            let mut r1 = EcPoint::default();
            let mut r2 = EcPoint::default();
            xdbladd(&mut r1, &mut r2, &p, &q, &pq, &a24, false);
            let mut sum = EcPoint::default();
            xadd(&mut sum, &p, &q, &pq);
            assert_ne!(ec_is_equal(&r2, &sum), 0, "xDBLADD add wrong");
            let mut dbl = EcPoint::default();
            ec_dbl(&mut dbl, &p, &curve);
            assert_ne!(ec_is_equal(&r1, &dbl), 0, "xDBLADD dbl wrong");
        }
    }

    #[test]
    fn xdbl_variants_agree() {
        let curve = make_e0();
        let mut a24 = EcPoint::default();
        ac_to_a24(&mut a24, &curve);
        let mut a24n = a24;
        ec_normalize_point(&mut a24n);
        let mut prng = Prng(0x1234);
        let z = prng.fp2();
        let mut a24r = EcPoint::default();
        fp2_mul(&mut a24r.x, &a24.x, &z);
        fp2_mul(&mut a24r.z, &a24.z, &z);

        for _ in 0..ITERS {
            let p = ec_random(&mut prng, &curve);
            let ac = EcPoint { x: curve.a, z: curve.c };
            let mut r1 = EcPoint::default();
            let mut r2 = EcPoint::default();
            let mut r3 = EcPoint::default();
            let mut r4 = EcPoint::default();
            xdbl(&mut r1, &p, &ac);
            xdbl_a24(&mut r2, &p, &a24r, false);
            xdbl_a24(&mut r3, &p, &a24n, true);
            xdbl_e0(&mut r4, &p);
            assert_ne!(ec_is_equal(&r1, &r2), 0);
            assert_ne!(ec_is_equal(&r1, &r3), 0);
            assert_ne!(ec_is_equal(&r1, &r4), 0);
        }
    }

    #[test]
    fn zero_identities() {
        let curve = make_e0();
        let mut prng = Prng(0xCAFE);
        let mut zero = EcPoint::default();
        ec_point_init(&mut zero);
        assert_ne!(ec_is_zero(&zero), 0);

        for _ in 0..ITERS {
            let p = ec_random(&mut prng, &curve);
            let mut r = EcPoint::default();
            xadd(&mut r, &zero, &zero, &zero);
            assert_ne!(ec_is_zero(&r), 0, "0+0 != 0");

            ec_dbl(&mut r, &p, &curve);
            let r2 = r;
            xadd(&mut r, &p, &p, &r2);
            assert_ne!(ec_is_zero(&r), 0, "P-P != 0");

            ec_dbl(&mut r, &zero, &curve);
            assert_ne!(ec_is_zero(&r), 0, "2·0 != 0");

            xadd(&mut r, &p, &zero, &p);
            assert_ne!(ec_is_equal(&r, &p), 0, "P+0 != P");
            xadd(&mut r, &zero, &p, &p);
            assert_ne!(ec_is_equal(&r, &p), 0, "0+P != P");

            let mut q = EcPoint::default();
            xdbladd(&mut r, &mut q, &p, &zero, &p, &curve.a24, false);
            assert_ne!(ec_is_equal(&q, &p), 0, "P+0 != P (xdbladd)");
            xdbladd(&mut r, &mut q, &zero, &p, &p, &curve.a24, false);
            assert_ne!(ec_is_equal(&q, &p), 0, "0+P != P (xdbladd)");
            assert_ne!(ec_is_zero(&r), 0, "2·0 != 0 (xdbladd)");
        }
    }

    #[test]
    fn jacobian_laws() {
        let curve = make_e0();
        let mut prng = Prng(0xDEAD);
        let mut jac_zero = JacPoint::default();
        jac_init(&mut jac_zero);

        for _ in 0..ITERS {
            let mut p = ec_random(&mut prng, &curve);
            ec_normalize_point(&mut p);
            let mut q = ec_random(&mut prng, &curve);
            ec_normalize_point(&mut q);

            let mut s = JacPoint::default();
            let mut t = JacPoint::default();
            fp2_copy(&mut s.x, &p.x);
            ec_recover_y(&mut s.y, &s.x, &curve);
            fp2_set_one(&mut s.z);
            fp2_copy(&mut t.x, &q.x);
            ec_recover_y(&mut t.y, &t.x, &curve);
            fp2_set_one(&mut t.z);

            let mut r = JacPoint::default();
            let mut u = JacPoint::default();

            jac_add(&mut r, &jac_zero, &jac_zero, &curve);
            assert_ne!(jac_is_equal(&r, &jac_zero), 0, "0+0 jac");

            jac_dbl(&mut r, &jac_zero, &curve);
            assert_ne!(jac_is_equal(&r, &jac_zero), 0, "2·0 jac");

            jac_neg(&mut r, &s);
            let r2 = r;
            jac_add(&mut r, &s, &r2, &curve);
            assert_ne!(jac_is_equal(&r, &jac_zero), 0, "P-P jac");

            jac_add(&mut r, &s, &jac_zero, &curve);
            assert_ne!(jac_is_equal(&r, &s), 0, "P+0 jac");
            jac_add(&mut r, &jac_zero, &s, &curve);
            assert_ne!(jac_is_equal(&r, &s), 0, "0+P jac");

            jac_dbl(&mut r, &s, &curve);
            jac_add(&mut u, &s, &s, &curve);
            assert_ne!(jac_is_equal(&r, &u), 0, "P+P=2P jac");

            jac_add(&mut r, &t, &s, &curve);
            let st = t;
            jac_add(&mut t, &s, &st, &curve);
            assert_ne!(jac_is_equal(&r, &t), 0, "comm jac");

            // Associativity
            let r2 = r;
            jac_dbl(&mut r, &r2, &curve);
            jac_add(&mut u, &s, &t, &curve);
            let su = u;
            jac_add(&mut u, &su, &r, &curve);
            let r2 = r;
            jac_add(&mut r, &r2, &t, &curve);
            let r2 = r;
            jac_add(&mut r, &r2, &s, &curve);
            assert_ne!(jac_is_equal(&r, &u), 0, "assoc jac");

            // ws roundtrip
            let mut t0 = Fp2::default();
            let mut ao3 = Fp2::default();
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            assert_ne!(jac_is_equal(&s, &r), 0, "ws roundtrip");

            let s2 = s;
            jac_dbl(&mut s, &s2, &curve);
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            assert_ne!(jac_is_equal(&s, &r), 0, "ws roundtrip 2");

            // DBLW
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            let st0 = t0;
            jac_dblw(&mut r, &mut t0, &r2, &st0);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            let s2 = s;
            jac_dbl(&mut s, &s2, &curve);
            assert_ne!(jac_is_equal(&s, &r), 0, "DBLW");
        }
    }

    #[test]
    fn ladder_dbl_matches_iterated() {
        let curve = make_e0();
        let mut prng = Prng(0x1111);
        let mut curve_m = curve;
        for _ in 0..10 {
            let p = ec_random(&mut prng, &curve);
            let four: [Digit; 1] = [4];
            let mut r1 = EcPoint::default();
            ec_mul(&mut r1, &four, 3, &p, &mut curve_m);
            let mut r2 = EcPoint::default();
            ec_dbl(&mut r2, &p, &curve);
            let s = r2;
            ec_dbl(&mut r2, &s, &curve);
            assert_ne!(ec_is_equal(&r1, &r2), 0, "[4]P via ladder != dbl·dbl");
        }
    }

    fn inner_test_generated_basis(basis: &EcBasis, curve: &EcCurve, n: usize) {
        let mut p = basis.p;
        let mut q = basis.q;
        for _ in 0..n - 1 {
            let (sp, sq) = (p, q);
            xdbl_a24(&mut p, &sp, &curve.a24, curve.is_a24_computed_and_normalized);
            xdbl_a24(&mut q, &sq, &curve.a24, curve.is_a24_computed_and_normalized);
        }
        assert_eq!(ec_is_zero(&p), 0, "P not full order");
        assert_eq!(ec_is_zero(&q), 0, "Q not full order");
        assert_eq!(ec_is_equal(&p, &q), 0, "P,Q dependent");
        assert_ne!(fp2_is_zero(&q.x), 0, "Q not above (0,0)");
        let (sp, sq) = (p, q);
        xdbl_a24(&mut p, &sp, &curve.a24, curve.is_a24_computed_and_normalized);
        xdbl_a24(&mut q, &sq, &curve.a24, curve.is_a24_computed_and_normalized);
        assert_ne!(ec_is_zero(&p), 0, "P order != 2^n");
        assert_ne!(ec_is_zero(&q), 0, "Q order != 2^n");
    }

    #[test]
    fn basis_generation_e0() {
        let mut curve = make_e0();
        let mut basis = EcBasis::default();
        ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, TORSION_EVEN_POWER as i32);
        inner_test_generated_basis(&basis, &curve, TORSION_EVEN_POWER);
    }

    #[test]
    fn basis_generation_a6() {
        let mut curve = EcCurve::default();
        ec_curve_init(&mut curve);
        fp2_set_small(&mut curve.a, 6);
        fp2_set_one(&mut curve.c);
        ec_curve_normalize_a24(&mut curve);

        let mut basis = EcBasis::default();
        let hint = ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, TORSION_EVEN_POWER as i32);
        inner_test_generated_basis(&basis, &curve, TORSION_EVEN_POWER);

        let mut basis2 = EcBasis::default();
        let ok = ec_curve_to_basis_2f_from_hint(&mut basis2, &mut curve, TORSION_EVEN_POWER as i32, hint);
        assert_eq!(ok, 1);
        assert_ne!(ec_is_equal(&basis.p, &basis2.p), 0);
        assert_ne!(ec_is_equal(&basis.q, &basis2.q), 0);
        assert_ne!(ec_is_equal(&basis.pmq, &basis2.pmq), 0);

        // partial order
        let mut basis3 = EcBasis::default();
        ec_curve_to_basis_2f_to_hint(&mut basis3, &mut curve, 128);
        inner_test_generated_basis(&basis3, &curve, 128);
    }

    #[test]
    fn lift_basis_roundtrip() {
        let mut curve = EcCurve::default();
        ec_curve_init(&mut curve);
        fp2_set_small(&mut curve.a, 6);
        ec_curve_normalize_a24(&mut curve);

        let mut basis = EcBasis::default();
        ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, TORSION_EVEN_POWER as i32);

        let mut jp = JacPoint::default();
        let mut jq = JacPoint::default();
        let r = lift_basis(&mut jp, &mut jq, &mut basis, &mut curve);
        assert_ne!(r, 0);

        // jac → xz should recover x(P), x(Q)
        let mut xp = EcPoint::default();
        jac_to_xz(&mut xp, &jp);
        assert_ne!(ec_is_equal(&xp, &basis.p), 0);
    }
}
