//! Deterministic 2ᶠ-torsion basis generation (entangled bases).
//! Port of `ec/ref/lvlx/basis.c`.

use super::*;
use crate::precomp::{
    BASIS_E0_PX, BASIS_E0_QX, P_COFACTOR_FOR_2F, P_COFACTOR_FOR_2F_BITLENGTH, TORSION_EVEN_POWER,
};

/// Recover y from x on y² = x³ + Ax² + x. Returns y iff x is on-curve.
pub fn ec_recover_y(px: &Fp2, curve: &EcCurve) -> Option<Fp2> {
    let t0 = px.square();
    (t0 * curve.a + px + t0 * px).sqrt_verify()
}

/// Deterministic x(P-Q) from x(P), x(Q) per ePrint 2017/518 Prop. 3.
fn difference_point(pq: &mut EcPoint, p: &EcPoint, q: &EcPoint, curve: &EcCurve) {
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

/// Lift an x-only basis to a Jacobian pair, assuming (A:C)=(A:1) and P=(X:1).
pub fn lift_basis_normalized(
    p: &mut JacPoint,
    q: &mut JacPoint,
    b: &mut EcBasis,
    e: &EcCurve,
) -> bool {
    debug_assert!(b.p.z.is_one());
    debug_assert!(e.c.is_one());

    p.x = b.p.x;
    q.x = b.q.x;
    q.z = b.q.z;
    p.z = Fp2::ONE;
    let Some(py) = ec_recover_y(&p.x, e) else {
        return false;
    };
    p.y = py;

    let mut v1 = p.x * q.z;
    let mut v2 = q.x + v1;
    let v3 = (q.x - v1).square() * b.pmq.x;
    v1 = e.a.dbl() * q.z;
    v2 += v1;
    let v4 = p.x * q.x + q.z;
    v2 *= v4;
    v1 *= q.z;
    v2 -= v1;
    v2 *= b.pmq.z;
    q.y = v3 - v2;
    v1 = p.y.dbl() * q.z * b.pmq.z;
    q.x *= v1;
    q.z *= v1;

    v1 = q.z.square();
    q.y *= v1;
    q.x *= q.z;
    true
}

/// Lift an x-only basis to a Jacobian pair, normalizing first.
/// Returns `false` if any input projective coordinate is zero (point at
/// infinity), since the batched inverse would otherwise silently zero `e.a`.
pub fn lift_basis(p: &mut JacPoint, q: &mut JacPoint, b: &mut EcBasis, e: &mut EcCurve) -> bool {
    if (b.p.z.is_zero_ct() | e.c.is_zero_ct()) != 0 {
        return false;
    }
    let mut inverses = [b.p.z, e.c];
    fp2_batched_inv(&mut inverses);
    b.p.z = Fp2::ONE;
    e.c = Fp2::ONE;
    b.p.x *= inverses[0];
    e.a *= inverses[1];
    lift_basis_normalized(p, q, b, e)
}

/// On-curve test for an affine x (assumes C=1).
fn is_on_curve(x: &Fp2, curve: &EcCurve) -> bool {
    debug_assert!(curve.c.is_one());
    let mut t0 = x + curve.a;
    t0 *= x;
    t0 = t0.add_one();
    t0 *= x;
    t0.is_square()
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
    let mut z;
    let mut t0;
    let mut t1;

    let cap = n.saturating_add(BASIS_SEARCH_CAP);
    loop {
        let mut qr_b = true;
        while qr_b {
            if n >= cap {
                return None;
            }
            qr_b = (Fp::from_small((n as u64) * (n as u64) + 1).is_square_ct()) != 0;
            n += 1;
        }

        let b = Fp::from_small((n - 1) as u64);
        z = Fp2 { re: Fp::ONE, im: b };
        t0 = Fp2 {
            re: Fp::ZERO,
            im: b,
        };

        t1 = curve.a.square();
        t0 *= t1;
        t1 = z.square();
        t0 -= t1;
        if t0.is_square_ct() == 0 {
            break;
        }
    }

    *x = z;
    *x = x.inv();
    *x *= curve.a;
    x.neg_ip();

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
        *x = curve.a;
    } else {
        *x = curve.a.mul_small(n as u32);
    }
    let cap = n.saturating_add(BASIS_SEARCH_CAP);
    while !is_on_curve(x, curve) {
        if n >= cap {
            return None;
        }
        *x += curve.a;
        n += 1;
    }
    Some(if n < 128 { n as u8 } else { 0 })
}

/// Precomputed E₀ basis specialised for A=0.
fn ec_basis_e0_2f(pq2: &mut EcBasis, curve: &EcCurve, f: i32) {
    debug_assert!(curve.a.is_zero());
    let mut p = EcPoint::default();
    let mut q = EcPoint::default();
    p.x = BASIS_E0_PX;
    q.x = BASIS_E0_QX;
    p.z = Fp2::ONE;
    q.z = Fp2::ONE;

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
    curve.normalize_and_a24();

    if curve.a.is_zero() {
        ec_basis_e0_2f(pq2, curve, f);
        return 0;
    }

    let hint_a = curve.a.is_square();

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

    p.z = Fp2::ONE;
    q.x = curve.a + p.x;
    q.x.neg_ip();
    q.z = Fp2::ONE;

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.q, &p, &q, curve);
    pq2.p = p;
    pq2.pmq = q;

    debug_assert!(hint < 128);
    (hint << 1) | (hint_a as u8)
}

/// Recompute the same E[2ᶠ] basis from a previously-returned hint.
/// Returns `None` if hint validation fails or the search cap is reached.
pub fn ec_curve_to_basis_2f_from_hint(curve: &mut EcCurve, f: i32, hint: u8) -> Option<EcBasis> {
    curve.normalize_and_a24();

    let mut pq2 = EcBasis::default();
    if curve.a.is_zero() {
        ec_basis_e0_2f(&mut pq2, curve, f);
        return Some(pq2);
    }

    let hint_a = (hint & 1) != 0;
    let hint_p = hint >> 1;

    let mut p = EcPoint::default();
    let mut q = EcPoint::default();

    if hint_p == 0 {
        if hint_a {
            find_nqr_factor(&mut p.x, curve, 128)?;
        } else {
            find_na_x_coord(&mut p.x, curve, 128)?;
        }
    } else if !hint_a {
        p.x = curve.a.mul_small(hint_p as u32);
    } else {
        p.x.re = Fp::ONE;
        p.x.im = Fp::from_small(hint_p as u64);
        p.x = (p.x).inv();
        p.x *= curve.a;
        p.x.neg_ip();
    }
    p.z = Fp2::ONE;

    #[cfg(debug_assertions)]
    if !is_on_curve(&p.x, curve) || p.x.is_square() {
        return None;
    }

    q.x = curve.a + p.x;
    q.x.neg_ip();
    q.z = Fp2::ONE;

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.q, &p, &q, curve);
    pq2.p = p;
    pq2.pmq = q;

    #[cfg(debug_assertions)]
    if !test_basis_order_twof(&pq2, curve, f) {
        return None;
    }

    Some(pq2)
}
