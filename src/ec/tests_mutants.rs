//! Targeted unit tests added to kill `cargo-mutants` survivors.
//! These exercise predicate negative cases, special-case branches, and
//! parameter-shape variation that the original tests didn't cover.

#![cfg(test)]

use super::*;
use crate::precomp::TORSION_EVEN_POWER;

fn e0_normalized() -> EcCurve {
    let mut e = EcCurve::e0();
    ec_curve_normalize_a24(&mut e);
    e
}

#[test]
fn validity_predicates() {
    let e = e0_normalized();
    let mut e1 = e;
    let basis = ec_curve_to_basis_2f_from_hint(&mut e1, TORSION_EVEN_POWER as i32, 0).unwrap();
    let p = basis.p;

    assert_ne!(ec_has_zero_coordinate(&EcPoint::default()), 0);
    assert_ne!(ec_has_zero_coordinate(&EcPoint::IDENTITY), 0);
    assert_eq!(ec_has_zero_coordinate(&p), 0);

    assert_ne!(ec_is_equal(&p, &p), 0);
    let mut scaled = p;
    let three = Fp2::from_small(3);
    let (sx, sz) = (scaled.x, scaled.z);
    fp2_mul(&mut scaled.x, &sx, &three);
    fp2_mul(&mut scaled.z, &sz, &three);
    assert_ne!(ec_is_equal(&p, &scaled), 0);
    assert_eq!(ec_is_equal(&p, &basis.q), 0);
    assert_ne!(ec_is_equal(&EcPoint::default(), &EcPoint::default()), 0);

    assert!(ec_curve_verify_a(&Fp2::ZERO));
    let mut a = Fp2::from_small(2);
    assert!(!ec_curve_verify_a(&a));
    let s = a;
    fp2_neg(&mut a, &s);
    assert!(!ec_curve_verify_a(&a));
    assert!(ec_curve_init_from_a(&Fp2::from_small(6)).is_some());

    let copy = basis;
    assert_ne!(ec_is_equal(&copy.p, &basis.p), 0);
    assert_ne!(ec_is_equal(&copy.q, &basis.q), 0);
    assert_ne!(ec_is_equal(&copy.pmq, &basis.pmq), 0);
}

#[test]
fn torsion_predicates() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();

    assert_eq!(ec_is_two_torsion(&basis.p, &e), 0);
    assert_eq!(ec_is_four_torsion(&basis.p, &e), 0);

    let mut b4 = EcBasis::default();
    ec_dbl_iter_basis(&mut b4, TORSION_EVEN_POWER as i32 - 2, &basis, &mut e);
    assert_ne!(ec_is_four_torsion(&b4.p, &e), 0);
    assert_ne!(ec_is_basis_four_torsion(&b4, &e), 0);

    let mut b2 = EcBasis::default();
    ec_dbl_iter_basis(&mut b2, 1, &b4, &mut e);
    assert_ne!(ec_is_two_torsion(&b2.p, &e), 0);

    let mut bad = b4;
    bad.q = b4.p;
    assert_eq!(ec_is_basis_four_torsion(&bad, &e), 0);
    assert_eq!(ec_is_basis_four_torsion(&basis, &e), 0);
}

#[test]
fn biscalar_mul_kbits1() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();
    let mut b2 = EcBasis::default();
    ec_dbl_iter_basis(&mut b2, TORSION_EVEN_POWER as i32 - 1, &basis, &mut e);

    let r = ec_biscalar_mul(&[0u64], &[0u64], 1, &b2, &e).unwrap();
    assert_ne!(ec_is_zero(&r), 0);
    let r = ec_biscalar_mul(&[1u64], &[0u64], 1, &b2, &e).unwrap();
    assert_ne!(ec_is_equal(&r, &b2.p), 0);
    let r = ec_biscalar_mul(&[0u64], &[1u64], 1, &b2, &e).unwrap();
    assert_ne!(ec_is_equal(&r, &b2.q), 0);
    let r = ec_biscalar_mul(&[1u64], &[1u64], 1, &b2, &e).unwrap();
    assert_ne!(ec_is_equal(&r, &b2.pmq), 0);

    assert!(ec_biscalar_mul(&[0u64], &[0u64], 1, &basis, &e).is_none());
    let mut deg = b2;
    deg.pmq.z = Fp2::default();
    assert!(ec_biscalar_mul(&[0u64], &[0u64], 1, &deg, &e).is_none());
}

#[test]
fn dbl_iter_normalize_threshold() {
    let mut e1 = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e1, TORSION_EVEN_POWER as i32, 0).unwrap();

    let mut e2 = EcCurve::e0();
    assert!(!e2.is_a24_computed_and_normalized);

    let mut r1 = EcPoint::default();
    let mut r2 = EcPoint::default();
    ec_dbl_iter(&mut r1, 51, &basis.p, &mut e1);
    ec_dbl_iter(&mut r2, 50, &basis.p, &mut e2);
    let s = r2;
    ec_dbl_iter(&mut r2, 1, &s, &mut e2);
    assert_ne!(ec_is_equal(&r1, &r2), 0);

    let mut r3 = EcPoint::default();
    let mut e3 = EcCurve::e0();
    ec_mul(&mut r3, &[1u64 << 51], 52, &basis.p, &mut e3);
    assert_ne!(ec_is_equal(&r1, &r3), 0);
}

#[test]
fn jac_predicates_and_init() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();
    let mut jp = JacPoint::default();
    let mut jq = JacPoint::default();
    let mut bcopy = basis;
    lift_basis(&mut jp, &mut jq, &mut bcopy, &mut e);

    assert!(jac_is_equal(&jp, &jp));
    assert!(!jac_is_equal(&jp, &jq));
    assert!(!jac_is_equal(&jp, &-jp));

    let id = JacPoint::IDENTITY;
    let idx = EcPoint::from(id);
    assert_ne!(fp2_is_one(&idx.x), 0);
    assert_ne!(fp2_is_zero(&idx.z), 0);

    let mut pxz = EcPoint::default();
    jac_to_xz(&mut pxz, &jp);
    let mut z2 = Fp2::default();
    fp2_sqr(&mut z2, &jp.z);
    assert_ne!(fp2_is_equal(&pxz.x, &jp.x), 0);
    assert_ne!(fp2_is_equal(&pxz.z, &z2), 0);

    let mut t = Fp2::default();
    let mut ao3 = Fp2::default();
    let mut idw = JacPoint::default();
    jac_to_ws(&mut idw, &mut t, &mut ao3, &id, &e);
    let mut dbl_id = JacPoint::default();
    let mut u = Fp2::default();
    jac_dblw(&mut dbl_id, &mut u, &idw, &t);
    assert_ne!(fp2_is_zero(&dbl_id.z), 0);
}

#[test]
fn singular_isogeny_consistency() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();

    let mut b24 = EcPoint::default();
    let mut kps = EcKps2::default();
    xisog_2_singular(&mut kps, &mut b24, e.a24);

    let mut img = [basis.p];
    xeval_2_singular_inplace(&mut img, &kps);
    assert_eq!(ec_is_zero(&img[0]), 0);
    let mut img2 = [EcPoint::default()];
    xeval_2_singular(&mut img2, &[basis.p], &kps);
    assert_ne!(ec_is_equal(&img[0], &img2[0]), 0);
}

#[test]
#[cfg(debug_assertions)]
fn basis_from_bad_hint_rejected() {
    let mut e = EcCurve::e0();
    e.a = Fp2::from_small(6);
    ec_curve_normalize_a24(&mut e);
    let mut basis = EcBasis::default();
    let valid = ec_curve_to_basis_2f_to_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32);
    let bad = valid ^ 1;
    assert!(ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, bad).is_none());
}

#[test]
fn xdbladd_nonnormalized_a24() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();

    let mut a24_unnorm = EcPoint::default();
    ac_to_a24(&mut a24_unnorm, &e);
    let three = Fp2::from_small(3);
    let (ax, az) = (a24_unnorm.x, a24_unnorm.z);
    fp2_mul(&mut a24_unnorm.x, &ax, &three);
    fp2_mul(&mut a24_unnorm.z, &az, &three);

    let mut r1 = EcPoint::default();
    let mut s1 = EcPoint::default();
    let mut r2 = EcPoint::default();
    let mut s2 = EcPoint::default();
    xdbladd(
        &mut r1, &mut s1, &basis.p, &basis.q, &basis.pmq, &e.a24, true,
    );
    xdbladd(
        &mut r2,
        &mut s2,
        &basis.p,
        &basis.q,
        &basis.pmq,
        &a24_unnorm,
        false,
    );
    assert_ne!(ec_is_equal(&r1, &r2), 0);
    assert_ne!(ec_is_equal(&s1, &s2), 0);
}
