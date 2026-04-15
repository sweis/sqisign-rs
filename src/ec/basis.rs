// SPDX-License-Identifier: Apache-2.0
//! Deterministic 2ᶠ-torsion basis generation (entangled bases).
//! Port of `ec/ref/lvlx/basis.c`.

use super::*;
use crate::precomp::{
    BASIS_E0_PX, BASIS_E0_QX, P_COFACTOR_FOR_2F, P_COFACTOR_FOR_2F_BITLENGTH, TORSION_EVEN_POWER,
};

/// Recover y from x on y² = x³ + Ax² + x. Returns 0xFFFFFFFF iff x was on-curve.
pub fn ec_recover_y(y: &mut Fp2, px: &Fp2, curve: &EcCurve) -> u32 {
    let mut t0 = Fp2::default();
    fp2_sqr(&mut t0, px);
    fp2_mul(y, &t0, &curve.a);
    fp2_add_ip(y, px);
    fp2_mul_ip(&mut t0, px);
    fp2_add_ip(y, &t0);
    fp2_sqrt_verify(y)
}

/// Deterministic x(P-Q) from x(P), x(Q) per ePrint 2017/518 Prop. 3.
fn difference_point(pq: &mut EcPoint, p: &EcPoint, q: &EcPoint, curve: &EcCurve) {
    let mut bxx = Fp2::default();
    let mut bxz = Fp2::default();
    let mut bzz = Fp2::default();
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();

    fp2_mul(&mut t0, &p.x, &q.x);
    fp2_mul(&mut t1, &p.z, &q.z);
    fp2_sub(&mut bxx, &t0, &t1);
    fp2_sqr_ip(&mut bxx);
    fp2_mul_ip(&mut bxx, &curve.c);
    fp2_add(&mut bxz, &t0, &t1);
    fp2_mul(&mut t0, &p.x, &q.z);
    fp2_mul(&mut t1, &p.z, &q.x);
    fp2_add(&mut bzz, &t0, &t1);
    fp2_mul_ip(&mut bxz, &bzz);
    fp2_sub(&mut bzz, &t0, &t1);
    fp2_sqr_ip(&mut bzz);
    fp2_mul_ip(&mut bzz, &curve.c);
    fp2_mul_ip(&mut bxz, &curve.c);
    fp2_mul_ip(&mut t0, &t1);
    fp2_mul_ip(&mut t0, &curve.a);
    fp2_dbl_ip(&mut t0);
    fp2_add_ip(&mut bxz, &t0);

    fp_copy(&mut t0.re, &curve.c.re);
    fp_neg(&mut t0.im, &curve.c.im);
    fp2_sqr_ip(&mut t0);
    fp2_mul_ip(&mut t0, &curve.c);
    fp_copy(&mut t1.re, &p.z.re);
    fp_neg(&mut t1.im, &p.z.im);
    fp2_sqr_ip(&mut t1);
    fp2_mul_ip(&mut t0, &t1);
    fp_copy(&mut t1.re, &q.z.re);
    fp_neg(&mut t1.im, &q.z.im);
    fp2_sqr_ip(&mut t1);
    fp2_mul_ip(&mut t0, &t1);
    fp2_mul_ip(&mut bxx, &t0);
    fp2_mul_ip(&mut bxz, &t0);
    fp2_mul_ip(&mut bzz, &t0);

    fp2_sqr(&mut t0, &bxz);
    fp2_mul(&mut t1, &bxx, &bzz);
    fp2_sub_ip(&mut t0, &t1);
    fp2_sqrt(&mut t0);
    fp2_add(&mut pq.x, &bxz, &t0);
    fp2_copy(&mut pq.z, &bzz);
}

/// Lift an x-only basis to a Jacobian pair, assuming (A:C)=(A:1) and P=(X:1).
pub fn lift_basis_normalized(
    p: &mut JacPoint,
    q: &mut JacPoint,
    b: &mut EcBasis,
    e: &EcCurve,
) -> u32 {
    debug_assert!(fp2_is_one(&b.p.z) != 0);
    debug_assert!(fp2_is_one(&e.c) != 0);

    fp2_copy(&mut p.x, &b.p.x);
    fp2_copy(&mut q.x, &b.q.x);
    fp2_copy(&mut q.z, &b.q.z);
    fp2_set_one(&mut p.z);
    let px = p.x;
    let ret = ec_recover_y(&mut p.y, &px, e);

    let mut v1 = Fp2::default();
    let mut v2 = Fp2::default();
    let mut v3 = Fp2::default();
    let mut v4 = Fp2::default();
    fp2_mul(&mut v1, &p.x, &q.z);
    let qx = q.x;
    fp2_add(&mut v2, &qx, &v1);
    fp2_sub(&mut v3, &qx, &v1);
    fp2_sqr_ip(&mut v3);
    fp2_mul_ip(&mut v3, &b.pmq.x);
    fp2_add(&mut v1, &e.a, &e.a);
    fp2_mul_ip(&mut v1, &q.z);
    fp2_add_ip(&mut v2, &v1);
    fp2_mul(&mut v4, &p.x, &q.x);
    fp2_add_ip(&mut v4, &q.z);
    fp2_mul_ip(&mut v2, &v4);
    fp2_mul_ip(&mut v1, &q.z);
    fp2_sub_ip(&mut v2, &v1);
    fp2_mul_ip(&mut v2, &b.pmq.z);
    fp2_sub(&mut q.y, &v3, &v2);
    fp2_add(&mut v1, &p.y, &p.y);
    fp2_mul_ip(&mut v1, &q.z);
    fp2_mul_ip(&mut v1, &b.pmq.z);
    let qx = q.x;
    fp2_mul(&mut q.x, &qx, &v1);
    let qz = q.z;
    fp2_mul(&mut q.z, &qz, &v1);

    let qz = q.z;
    fp2_sqr(&mut v1, &qz);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &v1);
    let qx = q.x;
    fp2_mul(&mut q.x, &qx, &qz);
    ret
}

/// Lift an x-only basis to a Jacobian pair, normalizing first.
/// Returns 0 if any input projective coordinate is zero (point at infinity);
/// otherwise the batched inverse would silently zero `e.a`.
pub fn lift_basis(p: &mut JacPoint, q: &mut JacPoint, b: &mut EcBasis, e: &mut EcCurve) -> u32 {
    if (fp2_is_zero(&b.p.z) | fp2_is_zero(&e.c)) != 0 {
        return 0;
    }
    let mut inverses = [b.p.z, e.c];
    fp2_batched_inv(&mut inverses);
    fp2_set_one(&mut b.p.z);
    fp2_set_one(&mut e.c);
    let bx = b.p.x;
    fp2_mul(&mut b.p.x, &bx, &inverses[0]);
    let ea = e.a;
    fp2_mul(&mut e.a, &ea, &inverses[1]);
    lift_basis_normalized(p, q, b, e)
}

/// On-curve test for an affine x (assumes C=1).
fn is_on_curve(x: &Fp2, curve: &EcCurve) -> u32 {
    debug_assert!(fp2_is_one(&curve.c) != 0);
    let mut t0 = Fp2::default();
    fp2_add(&mut t0, x, &curve.a);
    fp2_mul_ip(&mut t0, x);
    let s = t0;
    fp2_add_one(&mut t0, &s);
    fp2_mul_ip(&mut t0, x);
    fp2_is_square(&t0)
}

/// Clear odd cofactor and excess 2-power so that P has order exactly 2ᶠ.
#[inline]
fn clear_cofactor_for_maximal_even_order(p: &mut EcPoint, curve: &mut EcCurve, f: i32) {
    let s = *p;
    ec_mul(
        p,
        &P_COFACTOR_FOR_2F,
        P_COFACTOR_FOR_2F_BITLENGTH as i32,
        &s,
        curve,
    );
    for _ in 0..(TORSION_EVEN_POWER as i32 - f) {
        let s = *p;
        xdbl_a24(p, &s, &curve.a24, curve.is_a24_computed_and_normalized);
    }
}

/// Hard cap on basis-search iterations. Honest curves succeed within a few
/// iterations (failure probability < 2⁻¹²⁸), but an adversarial A on the
/// verify path can make `find_nqr_factor` non-terminating (e.g. Re(A²)=2
/// forces Im(t0)=0 so t0 is always a square in Fp²). The C reference loops
/// unbounded; we cap and return `None` so the caller can reject.
const BASIS_SEARCH_CAP: u16 = 256;

/// Find a NQR factor -1/(1+ib); returns hint b, or `None` if no NQR is found
/// within `BASIS_SEARCH_CAP` candidates (adversarial A).
fn find_nqr_factor(x: &mut Fp2, curve: &EcCurve, start: u8) -> Option<u8> {
    let mut n: u16 = start as u16;
    let mut b = Fp::default();
    let mut tmp = Fp::default();
    let mut z = Fp2::default();
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();

    let cap = n.saturating_add(BASIS_SEARCH_CAP);
    loop {
        let mut qr_b = true;
        while qr_b {
            if n >= cap {
                return None;
            }
            fp_set_small(&mut tmp, (n as u64) * (n as u64) + 1);
            qr_b = fp_is_square(&tmp) != 0;
            n += 1;
        }

        fp_set_small(&mut b, (n - 1) as u64);
        fp2_set_zero(&mut t0);
        fp2_set_one(&mut z);
        fp_copy(&mut z.im, &b);
        fp_copy(&mut t0.im, &b);

        fp2_sqr(&mut t1, &curve.a);
        fp2_mul_ip(&mut t0, &t1);
        fp2_sqr(&mut t1, &z);
        fp2_sub_ip(&mut t0, &t1);
        if fp2_is_square(&t0) == 0 {
            break;
        }
    }

    fp2_copy(x, &z);
    fp2_inv(x);
    fp2_mul_ip(x, &curve.a);
    fp2_neg_ip(x);

    Some(if n <= 128 { (n - 1) as u8 } else { 0 })
}

/// Find smallest n ≥ start with x(P) = n·A on-curve; returns hint n (or 0 on overflow).
///
/// The C reference asserts `!is_square(A)` here. That holds on the signing
/// path (caller picks the branch by `is_square(A)`), but on the verify path
/// `from_hint` reaches this with attacker-chosen A and hint, so we drop the
/// assert; the post-hoc validation in `ec_curve_to_basis_2f_from_hint`
/// rejects the resulting basis if it is malformed.
fn find_na_x_coord(x: &mut Fp2, curve: &EcCurve, start: u8) -> Option<u8> {
    let mut n: u16 = start as u16;
    if n == 1 {
        fp2_copy(x, &curve.a);
    } else {
        fp2_mul_small(x, &curve.a, n as u32);
    }
    let cap = n.saturating_add(BASIS_SEARCH_CAP);
    while is_on_curve(x, curve) == 0 {
        if n >= cap {
            return None;
        }
        fp2_add_ip(x, &curve.a);
        n += 1;
    }
    Some(if n < 128 { n as u8 } else { 0 })
}

/// Precomputed E₀ basis specialised for A=0.
fn ec_basis_e0_2f(pq2: &mut EcBasis, curve: &EcCurve, f: i32) {
    debug_assert!(fp2_is_zero(&curve.a) != 0);
    let mut p = EcPoint::default();
    let mut q = EcPoint::default();
    fp2_copy(&mut p.x, &BASIS_E0_PX);
    fp2_copy(&mut q.x, &BASIS_E0_QX);
    fp2_set_one(&mut p.z);
    fp2_set_one(&mut q.z);

    for _ in 0..(TORSION_EVEN_POWER as i32 - f) {
        let (sp, sq) = (p, q);
        xdbl_e0(&mut p, &sp);
        xdbl_e0(&mut q, &sq);
    }

    pq2.p = p;
    pq2.q = q;
    difference_point(&mut pq2.pmq, &p, &q, curve);
}

/// Deterministically compute an E[2ᶠ] basis with Q above (0,0); returns a hint byte.
pub fn ec_curve_to_basis_2f_to_hint(pq2: &mut EcBasis, curve: &mut EcCurve, f: i32) -> u8 {
    ec_normalize_curve_and_a24(curve);

    if fp2_is_zero(&curve.a) != 0 {
        ec_basis_e0_2f(pq2, curve, f);
        return 0;
    }

    let hint_a = fp2_is_square(&curve.a) != 0;

    let mut p = EcPoint::default();
    let mut q = EcPoint::default();

    let hint = if hint_a {
        find_nqr_factor(&mut p.x, curve, 1)
    } else {
        find_na_x_coord(&mut p.x, curve, 1)
    }
    // Signing path: caller supplies an honest curve, so the search succeeds
    // with probability > 1 − 2⁻²⁵⁶. If it ever fails the curve is degenerate.
    .expect("basis search failed on honest curve");

    fp2_set_one(&mut p.z);
    fp2_add(&mut q.x, &curve.a, &p.x);
    let qx = q.x;
    fp2_neg(&mut q.x, &qx);
    fp2_set_one(&mut q.z);

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.q, &p, &q, curve);
    pq2.p = p;
    pq2.pmq = q;

    debug_assert!(hint < 128);
    (hint << 1) | (hint_a as u8)
}

/// Recompute the same E[2ᶠ] basis from a previously-returned hint.
/// Returns 1 on success, 0 if debug-mode validation fails.
pub fn ec_curve_to_basis_2f_from_hint(
    pq2: &mut EcBasis,
    curve: &mut EcCurve,
    f: i32,
    hint: u8,
) -> i32 {
    ec_normalize_curve_and_a24(curve);

    if fp2_is_zero(&curve.a) != 0 {
        ec_basis_e0_2f(pq2, curve, f);
        return 1;
    }

    let hint_a = (hint & 1) != 0;
    let hint_p = hint >> 1;

    let mut p = EcPoint::default();
    let mut q = EcPoint::default();

    if hint_p == 0 {
        let ok = if hint_a {
            find_nqr_factor(&mut p.x, curve, 128)
        } else {
            find_na_x_coord(&mut p.x, curve, 128)
        };
        if ok.is_none() {
            return 0;
        }
    } else if !hint_a {
        fp2_mul_small(&mut p.x, &curve.a, hint_p as u32);
    } else {
        fp_set_one(&mut p.x.re);
        fp_set_small(&mut p.x.im, hint_p as u64);
        fp2_inv(&mut p.x);
        fp2_mul_ip(&mut p.x, &curve.a);
        fp2_neg_ip(&mut p.x);
    }
    fp2_set_one(&mut p.z);

    #[cfg(debug_assertions)]
    {
        let mut passed = is_on_curve(&p.x, curve) != 0;
        passed &= fp2_is_square(&p.x) == 0;
        if !passed {
            return 0;
        }
    }

    fp2_add(&mut q.x, &curve.a, &p.x);
    let qx = q.x;
    fp2_neg(&mut q.x, &qx);
    fp2_set_one(&mut q.z);

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.q, &p, &q, curve);
    pq2.p = p;
    pq2.pmq = q;

    #[cfg(debug_assertions)]
    {
        if test_basis_order_twof(pq2, curve, f) == 0 {
            return 0;
        }
    }

    1
}
