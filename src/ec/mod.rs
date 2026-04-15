//! Elliptic-curve arithmetic on Montgomery curves over GF(p²).
//!
//! Ported from `src/ec/ref/{include/ec.h, lvlx/*.c}` in the C reference.
//! x-only Kummer-line arithmetic, Jacobian-point arithmetic, 2ⁿ-isogeny
//! evaluation, and deterministic 2ᶠ-torsion basis generation.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use crate::gf::*;
use crate::mp::{self, Digit, LOG2RADIX, RADIX};
use crate::precomp::{BITS, NWORDS_ORDER};

pub mod basis;
pub mod biextension;
pub mod isog_chains;
pub mod jac;
mod tests_mutants;
pub mod xisog;

pub use basis::*;
pub use isog_chains::*;
pub use jac::*;
pub use xisog::*;

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

impl EcPoint {
    /// The Kummer-line identity (1 : 0).
    pub const IDENTITY: Self = Self {
        x: Fp2::ONE,
        z: Fp2::ZERO,
    };
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

impl EcCurve {
    /// The starting curve E₀ : y² = x³ + x, with (A : C) = (0 : 1).
    pub const fn e0() -> Self {
        Self {
            a: Fp2::ZERO,
            c: Fp2::ONE,
            a24: EcPoint::IDENTITY,
            is_a24_computed_and_normalized: false,
        }
    }
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

/// Compute A24 = (A+2C : 4C), or copy the cached normalized value.
#[inline]
pub fn ac_to_a24(a24: &mut EcPoint, e: &EcCurve) {
    if e.is_a24_computed_and_normalized {
        *a24 = e.a24;
        return;
    }
    a24.z = e.c + e.c;
    let z = a24.z;
    a24.x = e.a + z;
    a24.z = z + z;
}

/// Recover (A : C) = (2·(A24.x·2 - A24.z) : A24.z) from A24 = (A+2C : 4C).
#[inline]
pub fn a24_to_ac(e: &mut EcCurve, a24: &EcPoint) {
    e.a = (a24.x.dbl() - a24.z).dbl();
    e.c = a24.z;
}

// ===========================================================================
// ec.c port
// ===========================================================================

/// Constant-time select: Q ← P1 if option==0, else P2 (option must be 0 or all-ones).
#[inline]
pub fn select_point(q: &mut EcPoint, p1: &EcPoint, p2: &EcPoint, option: Digit) {
    let ctl = option as u32;
    q.x = Fp2::select(&p1.x, &p2.x, ctl);
    q.z = Fp2::select(&p1.z, &p2.z, ctl);
}

/// Constant-time conditional swap: if option is all-ones swap P↔Q.
#[inline]
pub fn cswap_points(p: &mut EcPoint, q: &mut EcPoint, option: Digit) {
    let ctl = option as u32;
    Fp2::cswap(&mut p.x, &mut q.x, ctl);
    Fp2::cswap(&mut p.z, &mut q.z, ctl);
}

impl EcPoint {
    /// Normalize (X : Z) → (X/Z : 1) in place.
    pub fn normalize(&mut self) {
        self.x *= self.z.inv();
        self.z = Fp2::ONE;
    }
    #[inline]
    pub fn is_zero_ct(&self) -> u32 {
        self.z.is_zero_ct()
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.z.is_zero()
    }
    #[inline]
    pub fn has_zero_coordinate_ct(&self) -> u32 {
        self.x.is_zero_ct() | self.z.is_zero_ct()
    }
    /// Projective equality: P = Q iff both zero, or neither zero and Px·Qz = Qx·Pz.
    pub fn is_equal_ct(&self, q: &EcPoint) -> u32 {
        let l_zero = self.is_zero_ct();
        let r_zero = q.is_zero_ct();
        let lr_equal = (self.x * q.z).is_equal_ct(&(self.z * q.x));
        // Faithful port of C's `(~l_zero & ~r_zero * lr_equal)`: `*` binds tighter than `&`.
        (l_zero & r_zero) | (!l_zero & (!r_zero).wrapping_mul(lr_equal))
    }
}

impl PartialEq for EcPoint {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal_ct(other) != 0
    }
}
impl Eq for EcPoint {}

impl EcCurve {
    /// Normalize (A : C) → (A/C : 1) in place.
    pub fn normalize(&mut self) {
        self.a *= self.c.inv();
        self.c = Fp2::ONE;
    }
    /// Compute and cache the normalized A24 = ((A+2C)/4C : 1).
    pub fn normalize_a24(&mut self) {
        if !self.is_a24_computed_and_normalized {
            let mut a24 = EcPoint::default();
            ac_to_a24(&mut a24, self);
            a24.normalize();
            self.a24 = a24;
            self.is_a24_computed_and_normalized = true;
        }
        debug_assert!(self.a24.z.is_one());
    }
    /// Normalize both (A : C) and A24 in place.
    pub fn normalize_and_a24(&mut self) {
        if !self.c.is_one() {
            self.normalize();
        }
        if !self.is_a24_computed_and_normalized {
            self.a24.x = Fp2 {
                re: self.a.re + Fp::from_small(2),
                im: self.a.im,
            }
            .half()
            .half();
            self.a24.z = Fp2::ONE;
            self.is_a24_computed_and_normalized = true;
        }
    }
    /// j-invariant.
    pub fn j_inv(&self) -> Fp2 {
        let cc = self.c.square();
        let t0 = self.a.square() - cc.mul_small(3);
        // 256·(A² − 3C²)³ / (C⁴·(A² − 4C²))
        let num = (t0.square() * t0).mul_small(256);
        num * (cc.square() * (t0 - cc)).inv()
    }
    /// True iff A is a valid Montgomery coefficient (A ≠ ±2).
    pub fn verify_a(a: &Fp2) -> bool {
        let two = Fp2::from_small(2);
        *a != two && *a != -two
    }
    /// Construct a curve from coefficient A, or `None` if A is invalid.
    pub fn try_from_a(a: &Fp2) -> Option<Self> {
        Self::verify_a(a).then(|| Self {
            a: *a,
            ..Self::e0()
        })
    }
    /// True iff P has exact 2-torsion (nonzero and 2P = 0).
    pub fn is_two_torsion_ct(&self, p: &EcPoint) -> u32 {
        if p.is_zero() {
            return 0;
        }
        let t0 = (p.x + p.z).square();
        let t1 = (p.x - p.z).square();
        let s = (t0 - t1) * self.a + (t0 + t1) * self.c.dbl();
        p.x.is_zero_ct() | s.is_zero_ct()
    }
    /// True iff P has exact 4-torsion.
    pub fn is_four_torsion_ct(&self, p: &EcPoint) -> u32 {
        let mut test = EcPoint::default();
        xdbl_a24(&mut test, p, &self.a24, self.is_a24_computed_and_normalized);
        self.is_two_torsion_ct(&test)
    }
    /// True iff (P, Q) form a full 4-basis.
    pub fn is_basis_four_torsion(&self, b: &EcBasis) -> bool {
        let mut p2 = EcPoint::default();
        let mut q2 = EcPoint::default();
        xdbl_a24(
            &mut p2,
            &b.p,
            &self.a24,
            self.is_a24_computed_and_normalized,
        );
        xdbl_a24(
            &mut q2,
            &b.q,
            &self.a24,
            self.is_a24_computed_and_normalized,
        );
        (self.is_two_torsion_ct(&p2) & self.is_two_torsion_ct(&q2)) != 0 && p2 != q2
    }
}

/// x-only doubling on E₀ (A=0, C=1).
pub fn xdbl_e0(q: &mut EcPoint, p: &EcPoint) {
    let mut t0 = p.x + p.z;
    t0.square_ip();
    let mut t1 = p.x - p.z;
    t1.square_ip();
    let t2 = t0 - t1;
    t1.dbl_ip();
    q.x = t0 * t1;
    q.z = t1 + t2;
    q.z *= t2;
}

/// x-only doubling, computing A24 from (A:C) on the fly.
pub fn xdbl(q: &mut EcPoint, p: &EcPoint, ac: &EcPoint) {
    let mut t0 = p.x + p.z;
    t0.square_ip();
    let mut t1 = p.x - p.z;
    t1.square_ip();
    let t2 = t0 - t1;
    let t3 = ac.z + ac.z;
    t1 *= t3;
    t1.dbl_ip();
    q.x = t0 * t1;
    t0 = t3 + ac.x;
    t0 *= t2;
    t0 += t1;
    q.z = t0 * t2;
}

/// x-only doubling using precomputed A24 = (A+2C : 4C) (or normalized).
pub fn xdbl_a24(q: &mut EcPoint, p: &EcPoint, a24: &EcPoint, a24_normalized: bool) {
    let mut t0 = p.x + p.z;
    t0.square_ip();
    let mut t1 = p.x - p.z;
    t1.square_ip();
    let t2 = t0 - t1;
    if !a24_normalized {
        t1 *= a24.z;
    }
    q.x = t0 * t1;
    t0 = t2 * a24.x;
    t0 += t1;
    q.z = t0 * t2;
}

/// Differential addition: R = P + Q given PQ = P - Q.
pub fn xadd(r: &mut EcPoint, p: &EcPoint, q: &EcPoint, pq: &EcPoint) {
    let mut t0 = p.x + p.z;
    let mut t1 = p.x - p.z;
    let mut t2 = q.x + q.z;
    let mut t3 = q.x - q.z;
    t0 *= t3;
    t1 *= t2;
    t2 = t0 + t1;
    t3 = t0 - t1;
    t2.square_ip();
    t3.square_ip();
    t2 *= pq.z;
    r.z = pq.x * t3;
    r.x = t2;
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
    let mut t0 = p.x + p.z;
    let mut t1 = p.x - p.z;
    r.x = t0.square();
    let mut t2 = q.x - q.z;
    s.x = q.x + q.z;
    t0 *= t2;
    r.z = t1.square();
    t1 *= s.x;
    let rx = r.x;
    t2 = rx - r.z;
    if !a24_normalized {
        r.z *= a24.z;
    }
    let (rx, rz) = (r.x, r.z);
    r.x = rx * rz;
    s.x = a24.x * t2;
    s.z = t0 - t1;
    let (rz, sx) = (r.z, s.x);
    r.z = rz + sx;
    s.x = t0 + t1;
    r.z *= t2;
    s.z.square_ip();
    s.x.square_ip();
    s.z *= pq.x;
    s.x *= pq.z;
}

/// Montgomery ladder: Q = k·P over `kbits` bits.
pub fn xmul(q: &mut EcPoint, p: &EcPoint, k: &[Digit], kbits: i32, curve: &EcCurve) {
    let mut a24 = EcPoint::default();
    if curve.is_a24_computed_and_normalized {
        a24.x = curve.a24.x;
        a24.z = curve.a24.z;
        debug_assert!(a24.z.is_one());
    } else {
        a24.x = curve.c + curve.c;
        let s = a24.x;
        a24.z = s + s;
        a24.x += curve.a;
    }

    let mut r0 = EcPoint::IDENTITY;
    let mut r1 = *p;

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

    q.x = r0.x;
    q.z = r0.z;
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
) -> Option<EcPoint> {
    if ((p).has_zero_coordinate_ct() | (q).has_zero_coordinate_ct() | (pq).has_zero_coordinate_ct())
        != 0
    {
        return None;
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

    let mut one = [0u64; NWORDS_ORDER];
    one[0] = 1;
    let mut k_t = [0u64; NWORDS_ORDER];
    let mut l_t = [0u64; NWORDS_ORDER];
    mp::mp_sub(&mut k_t, k, &one, NWORDS_ORDER);
    mp::mp_sub(&mut l_t, l, &one, NWORDS_ORDER);
    let kt = k_t;
    mp::select_ct(&mut k_t, &kt, k, maskk, NWORDS_ORDER);
    let lt = l_t;
    mp::select_ct(&mut l_t, &lt, l, maskl, NWORDS_ORDER);

    let mut r = [0u64; 2 * BITS];
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
        let d = (sigma[0] ^ sigma[1]) & maskk;
        sigma[0] ^= d;
        sigma[1] ^= d;
    }

    let mut rr = [EcPoint::default(); 3];
    rr[0] = EcPoint::IDENTITY;
    let maskk = 0u64.wrapping_sub(sigma[0]);
    select_point(&mut rr[1], p, q, maskk);
    select_point(&mut rr[2], q, p, maskk);

    let mut diff1a = rr[1];
    let mut diff1b = rr[2];

    let (r1, r2) = (rr[1], rr[2]);
    xadd(&mut rr[2], &r1, &r2, pq);
    if rr[2].has_zero_coordinate_ct() != 0 {
        return None;
    }
    let mut diff2a = rr[2];
    let mut diff2b = *pq;

    let a_is_zero = curve.a.is_zero();

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
            debug_assert!(curve.a24.z.is_one());
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
    Some(*s_out)
}

/// 3-point ladder: R = P + m·Q.
pub fn ec_ladder3pt(
    m: &[Digit],
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    e: &EcCurve,
) -> Option<EcPoint> {
    debug_assert!(e.is_a24_computed_and_normalized);
    if e.a24.z.is_one_ct() == 0 || (pq).has_zero_coordinate_ct() != 0 {
        return None;
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
    Some(x1)
}

/// Single doubling, choosing the cheapest formula given curve normalization.
pub fn ec_dbl(res: &mut EcPoint, p: &EcPoint, curve: &EcCurve) {
    if curve.is_a24_computed_and_normalized {
        debug_assert!(curve.a24.z.is_one());
        xdbl_a24(res, p, &curve.a24, true);
    } else {
        let ac = EcPoint {
            x: curve.a,
            z: curve.c,
        };
        xdbl(res, p, &ac);
    }
}

/// n-fold doubling.
pub fn ec_dbl_iter(res: &mut EcPoint, n: i32, p: &EcPoint, curve: &mut EcCurve) {
    if n == 0 {
        *res = *p;
        return;
    }
    if n > 50 {
        curve.normalize_a24();
    }
    if curve.is_a24_computed_and_normalized {
        debug_assert!(curve.a24.z.is_one());
        xdbl_a24(res, p, &curve.a24, true);
        for _ in 0..n - 1 {
            let s = *res;
            xdbl_a24(res, &s, &curve.a24, true);
        }
    } else {
        let ac = EcPoint {
            x: curve.a,
            z: curve.c,
        };
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
        curve.normalize_a24();
    }
    xmul(res, p, scalar, kbits, curve);
}

/// Biscalar multiplication: res = scalarP·P + scalarQ·Q.
pub fn ec_biscalar_mul(
    scalar_p: &[Digit],
    scalar_q: &[Digit],
    kbits: i32,
    pq: &EcBasis,
    curve: &EcCurve,
) -> Option<EcPoint> {
    if pq.pmq.z.is_zero() {
        return None;
    }
    if kbits == 1 {
        if curve.is_two_torsion_ct(&pq.p) == 0
            || curve.is_two_torsion_ct(&pq.q) == 0
            || curve.is_two_torsion_ct(&pq.pmq) == 0
        {
            return None;
        }
        return Some(match (scalar_p[0] & 1, scalar_q[0] & 1) {
            (0, 0) => EcPoint::IDENTITY,
            (1, 0) => pq.p,
            (0, 1) => pq.q,
            (1, 1) => pq.pmq,
            _ => unreachable!(),
        });
    }
    let mut e = *curve;
    if curve.a.is_zero_ct() == 0 {
        e.normalize_a24();
    }
    let mut res = EcPoint::default();
    xdblmul(
        &mut res, &pq.p, scalar_p, &pq.q, scalar_q, &pq.pmq, kbits, &e,
    )
}

// ===========================================================================
// Debug / test helpers (from ec.h static inlines)
// ===========================================================================

/// True iff P has order exactly 2^t.
pub fn test_point_order_twof(p: &EcPoint, e: &EcCurve, t: i32) -> bool {
    let mut test = *p;
    let mut curve = *e;
    if test.is_zero_ct() != 0 {
        return false;
    }
    let pt = test;
    ec_dbl_iter(&mut test, t - 1, &pt, &mut curve);
    if test.is_zero_ct() != 0 {
        return false;
    }
    let pt = test;
    ec_dbl(&mut test, &pt, &curve);
    test.is_zero_ct() != 0
}

/// True iff all of P, Q, P−Q have order exactly 2^t.
pub fn test_basis_order_twof(b: &EcBasis, e: &EcCurve, t: i32) -> bool {
    test_point_order_twof(&b.p, e, t)
        && test_point_order_twof(&b.q, e, t)
        && test_point_order_twof(&b.pmq, e, t)
}

/// True iff Jacobian P has order exactly 2^t.
pub fn test_jac_order_twof(p: &JacPoint, e: &EcCurve, t: i32) -> bool {
    let mut test = *p;
    if test.z.is_zero() {
        return false;
    }
    for _ in 0..t - 1 {
        let s = test;
        jac_dbl(&mut test, &s, e);
    }
    if test.z.is_zero() {
        return false;
    }
    let s = test;
    jac_dbl(&mut test, &s, e);
    test.z.is_zero()
}

// ===========================================================================
// Tests (port of curve-arith-test.c, basis-gen-test.c)
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precomp::TORSION_EVEN_POWER;
    use crate::test_util::Prng;

    fn projective_is_on_curve(p: &EcPoint, curve: &EcCurve) -> bool {
        let mut t0 = curve.c * p.x;
        let mut t1 = t0 * p.z;
        t1 *= curve.a;
        let mut t2 = curve.c * p.z;
        t0.square_ip();
        t2.square_ip();
        t0 += t1;
        t0 += t2;
        t0 *= p.x;
        t0 *= p.z;
        t0.is_square() || t0.is_zero()
    }

    fn ec_random(prng: &mut Prng, curve: &EcCurve) -> EcPoint {
        let mut p = EcPoint::default();
        p.z = Fp2::ONE;
        loop {
            p.x = prng.fp2();
            if projective_is_on_curve(&p, curve) {
                break;
            }
        }
        let z = prng.fp2();
        p.x *= z;
        p.z = z;
        p
    }

    /// Port of test_extras.c::projective_difference_point.
    fn projective_difference_point(pq: &mut EcPoint, p: &EcPoint, q: &EcPoint, curve: &EcCurve) {
        let mut t0 = p.x * q.x;
        let mut t1 = p.z * q.z;
        let mut bxx = t0 - t1;
        bxx.square_ip();
        bxx *= curve.c;
        let mut bxz = t0 + t1;
        t0 = p.x * q.z;
        t1 = p.z * q.x;
        let mut bzz = t0 + t1;
        bxz *= bzz;
        bzz = t0 - t1;
        bzz.square_ip();
        bzz *= curve.c;
        bxz *= curve.c;
        t0 *= t1;
        t0 *= curve.a;
        t0.dbl_ip();
        bxz += t0;

        t0.re = curve.c.re;
        t0.im = -curve.c.im;
        t0.square_ip();
        t0 *= curve.c;
        t1.re = p.z.re;
        t1.im = -p.z.im;
        t1.square_ip();
        t0 *= t1;
        t1.re = q.z.re;
        t1.im = -q.z.im;
        t1.square_ip();
        t0 *= t1;
        bxx *= t0;
        bxz *= t0;
        bzz *= t0;

        t0 = bxz.square();
        t1 = bxx * bzz;
        t0 -= t1;
        t0 = t0.sqrt();
        pq.x = bxz + t0;
        pq.z = bzz;
    }

    fn make_e0() -> EcCurve {
        let mut curve = EcCurve::e0();
        curve.normalize_a24();
        curve
    }

    const ITERS: usize = 50;

    #[test]
    fn j_invariant_e0_is_1728() {
        let curve = make_e0();
        let j = curve.j_inv();
        let expected = Fp2::from_small(1728);
        assert_ne!(j.is_equal_ct(&expected), 0, "j(E0) != 1728");
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
            assert_ne!(r1.is_equal_ct(&r2), 0, "2(P+Q) != 2P+2Q");

            // (P+Q)+(P-Q) = 2P
            xadd(&mut r1, &p, &q, &pq);
            let sq = q;
            ec_dbl(&mut q, &sq, &curve);
            let sr1 = r1;
            xadd(&mut r1, &sr1, &pq, &q);
            let sp = p;
            ec_dbl(&mut p, &sp, &curve);
            assert_ne!(r1.is_equal_ct(&p), 0, "(P+Q)+(P-Q) != 2P");
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
            assert_ne!(r2.is_equal_ct(&sum), 0, "xDBLADD add wrong");
            let mut dbl = EcPoint::default();
            ec_dbl(&mut dbl, &p, &curve);
            assert_ne!(r1.is_equal_ct(&dbl), 0, "xDBLADD dbl wrong");
        }
    }

    #[test]
    fn xdbl_variants_agree() {
        let curve = make_e0();
        let mut a24 = EcPoint::default();
        ac_to_a24(&mut a24, &curve);
        let mut a24n = a24;
        a24n.normalize();
        let mut prng = Prng(0x1234);
        let z = prng.fp2();
        let mut a24r = EcPoint::default();
        a24r.x = a24.x * z;
        a24r.z = a24.z * z;

        for _ in 0..ITERS {
            let p = ec_random(&mut prng, &curve);
            let ac = EcPoint {
                x: curve.a,
                z: curve.c,
            };
            let mut r1 = EcPoint::default();
            let mut r2 = EcPoint::default();
            let mut r3 = EcPoint::default();
            let mut r4 = EcPoint::default();
            xdbl(&mut r1, &p, &ac);
            xdbl_a24(&mut r2, &p, &a24r, false);
            xdbl_a24(&mut r3, &p, &a24n, true);
            xdbl_e0(&mut r4, &p);
            assert_ne!(r1.is_equal_ct(&r2), 0);
            assert_ne!(r1.is_equal_ct(&r3), 0);
            assert_ne!(r1.is_equal_ct(&r4), 0);
        }
    }

    #[test]
    fn zero_identities() {
        let curve = make_e0();
        let mut prng = Prng(0xCAFE);
        let zero = EcPoint::IDENTITY;
        assert_ne!(zero.is_zero_ct(), 0);

        for _ in 0..ITERS {
            let p = ec_random(&mut prng, &curve);
            let mut r = EcPoint::default();
            xadd(&mut r, &zero, &zero, &zero);
            assert_ne!(r.is_zero_ct(), 0, "0+0 != 0");

            ec_dbl(&mut r, &p, &curve);
            let r2 = r;
            xadd(&mut r, &p, &p, &r2);
            assert_ne!(r.is_zero_ct(), 0, "P-P != 0");

            ec_dbl(&mut r, &zero, &curve);
            assert_ne!(r.is_zero_ct(), 0, "2·0 != 0");

            xadd(&mut r, &p, &zero, &p);
            assert_ne!(r.is_equal_ct(&p), 0, "P+0 != P");
            xadd(&mut r, &zero, &p, &p);
            assert_ne!(r.is_equal_ct(&p), 0, "0+P != P");

            let mut q = EcPoint::default();
            xdbladd(&mut r, &mut q, &p, &zero, &p, &curve.a24, false);
            assert_ne!(q.is_equal_ct(&p), 0, "P+0 != P (xdbladd)");
            xdbladd(&mut r, &mut q, &zero, &p, &p, &curve.a24, false);
            assert_ne!(q.is_equal_ct(&p), 0, "0+P != P (xdbladd)");
            assert_ne!(r.is_zero_ct(), 0, "2·0 != 0 (xdbladd)");
        }
    }

    #[test]
    fn jacobian_laws() {
        let curve = make_e0();
        let mut prng = Prng(0xDEAD);
        let jac_zero = JacPoint::IDENTITY;

        for _ in 0..ITERS {
            let mut p = ec_random(&mut prng, &curve);
            p.normalize();
            let mut q = ec_random(&mut prng, &curve);
            q.normalize();

            let mut s = JacPoint {
                x: p.x,
                y: ec_recover_y(&p.x, &curve).unwrap(),
                z: Fp2::ONE,
            };
            let t = JacPoint {
                x: q.x,
                y: ec_recover_y(&q.x, &curve).unwrap(),
                z: Fp2::ONE,
            };

            let mut r = JacPoint::default();
            let mut u = JacPoint::default();

            jac_add(&mut r, &jac_zero, &jac_zero, &curve);
            assert!(jac_is_equal(&r, &jac_zero), "0+0 jac");

            jac_dbl(&mut r, &jac_zero, &curve);
            assert!(jac_is_equal(&r, &jac_zero), "2·0 jac");

            r = -s;
            let r2 = r;
            jac_add(&mut r, &s, &r2, &curve);
            assert!(jac_is_equal(&r, &jac_zero), "P-P jac");

            jac_add(&mut r, &s, &jac_zero, &curve);
            assert!(jac_is_equal(&r, &s), "P+0 jac");
            jac_add(&mut r, &jac_zero, &s, &curve);
            assert!(jac_is_equal(&r, &s), "0+P jac");

            jac_dbl(&mut r, &s, &curve);
            jac_add(&mut u, &s, &s, &curve);
            assert!(jac_is_equal(&r, &u), "P+P=2P jac");

            jac_add(&mut r, &t, &s, &curve);
            let mut t2 = JacPoint::default();
            jac_add(&mut t2, &s, &t, &curve);
            assert!(jac_is_equal(&r, &t2), "comm jac");

            // Associativity
            let r2 = r;
            jac_dbl(&mut r, &r2, &curve);
            jac_add(&mut u, &s, &t2, &curve);
            let su = u;
            jac_add(&mut u, &su, &r, &curve);
            let r2 = r;
            jac_add(&mut r, &r2, &t2, &curve);
            let r2 = r;
            jac_add(&mut r, &r2, &s, &curve);
            assert!(jac_is_equal(&r, &u), "assoc jac");

            // ws roundtrip
            let mut t0 = Fp2::default();
            let mut ao3 = Fp2::default();
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            assert!(jac_is_equal(&s, &r), "ws roundtrip");

            let s2 = s;
            jac_dbl(&mut s, &s2, &curve);
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            assert!(jac_is_equal(&s, &r), "ws roundtrip 2");

            // DBLW
            jac_to_ws(&mut r, &mut t0, &mut ao3, &s, &curve);
            let r2 = r;
            let st0 = t0;
            jac_dblw(&mut r, &mut t0, &r2, &st0);
            let r2 = r;
            jac_from_ws(&mut r, &r2, &ao3, &curve);
            let s2 = s;
            jac_dbl(&mut s, &s2, &curve);
            assert!(jac_is_equal(&s, &r), "DBLW");
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
            assert_ne!(r1.is_equal_ct(&r2), 0, "[4]P via ladder != dbl·dbl");
        }
    }

    fn inner_test_generated_basis(basis: &EcBasis, curve: &EcCurve, n: usize) {
        let mut p = basis.p;
        let mut q = basis.q;
        for _ in 0..n - 1 {
            let (sp, sq) = (p, q);
            xdbl_a24(
                &mut p,
                &sp,
                &curve.a24,
                curve.is_a24_computed_and_normalized,
            );
            xdbl_a24(
                &mut q,
                &sq,
                &curve.a24,
                curve.is_a24_computed_and_normalized,
            );
        }
        assert_eq!(p.is_zero_ct(), 0, "P not full order");
        assert_eq!(q.is_zero_ct(), 0, "Q not full order");
        assert_eq!(p.is_equal_ct(&q), 0, "P,Q dependent");
        assert_ne!(q.x.is_zero_ct(), 0, "Q not above (0,0)");
        let (sp, sq) = (p, q);
        xdbl_a24(
            &mut p,
            &sp,
            &curve.a24,
            curve.is_a24_computed_and_normalized,
        );
        xdbl_a24(
            &mut q,
            &sq,
            &curve.a24,
            curve.is_a24_computed_and_normalized,
        );
        assert_ne!(p.is_zero_ct(), 0, "P order != 2^n");
        assert_ne!(q.is_zero_ct(), 0, "Q order != 2^n");
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
        let mut curve = EcCurve::e0();
        curve.a = Fp2::from_small(6);
        curve.c = Fp2::ONE;
        curve.normalize_a24();

        let mut basis = EcBasis::default();
        let hint = ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, TORSION_EVEN_POWER as i32);
        inner_test_generated_basis(&basis, &curve, TORSION_EVEN_POWER);

        let basis2 =
            ec_curve_to_basis_2f_from_hint(&mut curve, TORSION_EVEN_POWER as i32, hint).unwrap();
        assert_ne!(basis.p.is_equal_ct(&basis2.p), 0);
        assert_ne!(basis.q.is_equal_ct(&basis2.q), 0);
        assert_ne!(basis.pmq.is_equal_ct(&basis2.pmq), 0);

        // partial order
        let mut basis3 = EcBasis::default();
        ec_curve_to_basis_2f_to_hint(&mut basis3, &mut curve, 128);
        inner_test_generated_basis(&basis3, &curve, 128);
    }

    #[test]
    fn lift_basis_roundtrip() {
        let mut curve = EcCurve::e0();
        curve.a = Fp2::from_small(6);
        curve.normalize_a24();

        let mut basis = EcBasis::default();
        ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, TORSION_EVEN_POWER as i32);

        let mut jp = JacPoint::default();
        let mut jq = JacPoint::default();
        assert!(lift_basis(&mut jp, &mut jq, &mut basis, &mut curve));

        // jac → xz should recover x(P), x(Q)
        let mut xp = EcPoint::default();
        jac_to_xz(&mut xp, &jp);
        assert_ne!(xp.is_equal_ct(&basis.p), 0);
    }
}
