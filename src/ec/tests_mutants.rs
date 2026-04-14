// SPDX-License-Identifier: Apache-2.0
//! Targeted unit tests added to kill `cargo-mutants` survivors.
//! These exercise predicate negative cases, special-case branches, and
//! parameter-shape variation that the original tests didn't cover.

#![cfg(test)]

use super::*;
use crate::precomp::TORSION_EVEN_POWER;

fn e0_normalized() -> EcCurve {
    let mut e = EcCurve::default();
    ec_curve_init(&mut e);
    ec_curve_normalize_a24(&mut e);
    e
}

#[test]
fn validity_predicates() {
    let e = e0_normalized();
    let mut basis = EcBasis::default();
    let mut e1 = e;
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e1, TORSION_EVEN_POWER as i32, 0);
    let p = basis.p;

    let mut zero = EcPoint::default();
    assert_ne!(ec_has_zero_coordinate(&zero), 0);
    ec_point_init(&mut zero);
    assert_ne!(ec_has_zero_coordinate(&zero), 0);
    assert_eq!(ec_has_zero_coordinate(&p), 0);

    assert_ne!(ec_is_equal(&p, &p), 0);
    let mut scaled = p;
    let mut three = Fp2::default();
    fp2_set_small(&mut three, 3);
    let (sx, sz) = (scaled.x, scaled.z);
    fp2_mul(&mut scaled.x, &sx, &three);
    fp2_mul(&mut scaled.z, &sz, &three);
    assert_ne!(ec_is_equal(&p, &scaled), 0);
    assert_eq!(ec_is_equal(&p, &basis.q), 0);
    assert_ne!(ec_is_equal(&EcPoint::default(), &EcPoint::default()), 0);

    let mut a = Fp2::default();
    assert_eq!(ec_curve_verify_a(&a), 1);
    fp2_set_small(&mut a, 2);
    assert_eq!(ec_curve_verify_a(&a), 0);
    let s = a;
    fp2_neg(&mut a, &s);
    assert_eq!(ec_curve_verify_a(&a), 0);
    fp2_set_small(&mut a, 6);
    let mut e2 = EcCurve::default();
    assert_eq!(ec_curve_init_from_a(&mut e2, &a), 1);

    let mut copy = EcBasis::default();
    copy_basis(&mut copy, &basis);
    assert_ne!(ec_is_equal(&copy.p, &basis.p), 0);
    assert_ne!(ec_is_equal(&copy.q, &basis.q), 0);
    assert_ne!(ec_is_equal(&copy.pmq, &basis.pmq), 0);
}

#[test]
fn torsion_predicates() {
    let mut e = e0_normalized();
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, 0);

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
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, 0);
    let mut b2 = EcBasis::default();
    ec_dbl_iter_basis(&mut b2, TORSION_EVEN_POWER as i32 - 1, &basis, &mut e);

    let mut r = EcPoint::default();
    assert_eq!(ec_biscalar_mul(&mut r, &[0u64], &[0u64], 1, &b2, &e), 1);
    assert_ne!(ec_is_zero(&r), 0);
    assert_eq!(ec_biscalar_mul(&mut r, &[1u64], &[0u64], 1, &b2, &e), 1);
    assert_ne!(ec_is_equal(&r, &b2.p), 0);
    assert_eq!(ec_biscalar_mul(&mut r, &[0u64], &[1u64], 1, &b2, &e), 1);
    assert_ne!(ec_is_equal(&r, &b2.q), 0);
    assert_eq!(ec_biscalar_mul(&mut r, &[1u64], &[1u64], 1, &b2, &e), 1);
    assert_ne!(ec_is_equal(&r, &b2.pmq), 0);

    assert_eq!(ec_biscalar_mul(&mut r, &[0u64], &[0u64], 1, &basis, &e), 0);
    let mut deg = b2;
    deg.pmq.z = Fp2::default();
    assert_eq!(ec_biscalar_mul(&mut r, &[0u64], &[0u64], 1, &deg, &e), 0);
}

#[test]
fn dbl_iter_normalize_threshold() {
    let mut e1 = e0_normalized();
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e1, TORSION_EVEN_POWER as i32, 0);

    let mut e2 = EcCurve::default();
    ec_curve_init(&mut e2);
    assert!(!e2.is_a24_computed_and_normalized);

    let mut r1 = EcPoint::default();
    let mut r2 = EcPoint::default();
    ec_dbl_iter(&mut r1, 51, &basis.p, &mut e1);
    ec_dbl_iter(&mut r2, 50, &basis.p, &mut e2);
    let s = r2;
    ec_dbl_iter(&mut r2, 1, &s, &mut e2);
    assert_ne!(ec_is_equal(&r1, &r2), 0);

    let mut r3 = EcPoint::default();
    let mut e3 = EcCurve::default();
    ec_curve_init(&mut e3);
    ec_mul(&mut r3, &[1u64 << 51], 52, &basis.p, &mut e3);
    assert_ne!(ec_is_equal(&r1, &r3), 0);
}

#[test]
fn jac_predicates_and_init() {
    let mut e = e0_normalized();
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, 0);
    let mut jp = JacPoint::default();
    let mut jq = JacPoint::default();
    let mut bcopy = basis;
    lift_basis(&mut jp, &mut jq, &mut bcopy, &mut e);

    assert_ne!(jac_is_equal(&jp, &jp), 0);
    assert_eq!(jac_is_equal(&jp, &jq), 0);
    let mut neg_p = jp;
    let py = neg_p.y;
    fp2_neg(&mut neg_p.y, &py);
    assert_eq!(jac_is_equal(&jp, &neg_p), 0);

    let mut id = jp;
    jac_init(&mut id);
    let mut idx = EcPoint::default();
    jac_to_xz(&mut idx, &id);
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
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, 0);

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
    let mut e = EcCurve::default();
    ec_curve_init(&mut e);
    fp2_set_small(&mut e.a, 6);
    ec_curve_normalize_a24(&mut e);
    let mut basis = EcBasis::default();
    let valid = ec_curve_to_basis_2f_to_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32);
    let bad = valid ^ 1;
    let ok = ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, bad);
    assert_eq!(ok, 0);
}

#[test]
fn xdbladd_nonnormalized_a24() {
    let mut e = e0_normalized();
    let mut basis = EcBasis::default();
    ec_curve_to_basis_2f_from_hint(&mut basis, &mut e, TORSION_EVEN_POWER as i32, 0);

    let mut a24_unnorm = EcPoint::default();
    ac_to_a24(&mut a24_unnorm, &e);
    let mut three = Fp2::default();
    fp2_set_small(&mut three, 3);
    let (ax, az) = (a24_unnorm.x, a24_unnorm.z);
    fp2_mul(&mut a24_unnorm.x, &ax, &three);
    fp2_mul(&mut a24_unnorm.z, &az, &three);

    let mut r1 = EcPoint::default();
    let mut s1 = EcPoint::default();
    let mut r2 = EcPoint::default();
    let mut s2 = EcPoint::default();
    xdbladd(&mut r1, &mut s1, &basis.p, &basis.q, &basis.pmq, &e.a24, true);
    xdbladd(&mut r2, &mut s2, &basis.p, &basis.q, &basis.pmq, &a24_unnorm, false);
    assert_ne!(ec_is_equal(&r1, &r2), 0);
    assert_ne!(ec_is_equal(&s1, &s2), 0);
}
