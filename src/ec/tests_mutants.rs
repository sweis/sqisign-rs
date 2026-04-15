//! Targeted unit tests added to kill `cargo-mutants` survivors.
//! These exercise predicate negative cases, special-case branches, and
//! parameter-shape variation that the original tests didn't cover.

#![cfg(test)]

use super::*;
use crate::precomp::TORSION_EVEN_POWER;

fn e0_normalized() -> EcCurve {
    let mut e = EcCurve::e0();
    e.normalize_a24();
    e
}

#[test]
fn validity_predicates() {
    let e = e0_normalized();
    let mut e1 = e;
    let basis = ec_curve_to_basis_2f_from_hint(&mut e1, TORSION_EVEN_POWER as i32, 0).unwrap();
    let p = basis.p;

    assert_ne!(EcPoint::default().has_zero_coordinate_ct(), 0);
    assert_ne!(EcPoint::IDENTITY.has_zero_coordinate_ct(), 0);
    assert_eq!(p.has_zero_coordinate_ct(), 0);

    assert_eq!(p, p);
    let three = Fp2::from_small(3);
    let scaled = EcPoint {
        x: p.x * three,
        z: p.z * three,
    };
    assert_eq!(p, scaled);
    assert_ne!(p, basis.q);
    assert_eq!(EcPoint::default(), EcPoint::default());

    assert!(EcCurve::verify_a(&Fp2::ZERO));
    assert!(!EcCurve::verify_a(&Fp2::from_small(2)));
    assert!(!EcCurve::verify_a(&(-Fp2::from_small(2))));
    assert!(EcCurve::try_from_a(&Fp2::from_small(6)).is_some());

    let copy = basis;
    assert_eq!(copy.p, basis.p);
    assert_eq!(copy.q, basis.q);
    assert_eq!(copy.pmq, basis.pmq);
}

#[test]
fn torsion_predicates() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();

    assert_eq!(e.is_two_torsion_ct(&basis.p), 0);
    assert_eq!(e.is_four_torsion_ct(&basis.p), 0);

    let mut b4 = EcBasis::default();
    ec_dbl_iter_basis(&mut b4, TORSION_EVEN_POWER as i32 - 2, &basis, &mut e);
    assert_ne!(e.is_four_torsion_ct(&b4.p), 0);
    assert!(e.is_basis_four_torsion(&b4));

    let mut b2 = EcBasis::default();
    ec_dbl_iter_basis(&mut b2, 1, &b4, &mut e);
    assert_ne!(e.is_two_torsion_ct(&b2.p), 0);

    let mut bad = b4;
    bad.q = b4.p;
    assert!(!e.is_basis_four_torsion(&bad));
    assert!(!e.is_basis_four_torsion(&basis));
}

#[test]
fn biscalar_mul_kbits1() {
    let mut e = e0_normalized();
    let basis = ec_curve_to_basis_2f_from_hint(&mut e, TORSION_EVEN_POWER as i32, 0).unwrap();
    let mut b2 = EcBasis::default();
    ec_dbl_iter_basis(&mut b2, TORSION_EVEN_POWER as i32 - 1, &basis, &mut e);

    let r = ec_biscalar_mul(&[0u64], &[0u64], 1, &b2, &e).unwrap();
    assert_ne!(r.is_zero_ct(), 0);
    let r = ec_biscalar_mul(&[1u64], &[0u64], 1, &b2, &e).unwrap();
    assert_ne!(r.is_equal_ct(&b2.p), 0);
    let r = ec_biscalar_mul(&[0u64], &[1u64], 1, &b2, &e).unwrap();
    assert_ne!(r.is_equal_ct(&b2.q), 0);
    let r = ec_biscalar_mul(&[1u64], &[1u64], 1, &b2, &e).unwrap();
    assert_ne!(r.is_equal_ct(&b2.pmq), 0);

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
    assert_ne!(r1.is_equal_ct(&r2), 0);

    let mut r3 = EcPoint::default();
    let mut e3 = EcCurve::e0();
    ec_mul(&mut r3, &[1u64 << 51], 52, &basis.p, &mut e3);
    assert_ne!(r1.is_equal_ct(&r3), 0);
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
    assert_ne!(idx.x.is_one_ct(), 0);
    assert_ne!(idx.z.is_zero_ct(), 0);

    let mut pxz = EcPoint::default();
    jac_to_xz(&mut pxz, &jp);
    let z2 = jp.z.square();
    assert_ne!(pxz.x.is_equal_ct(&jp.x), 0);
    assert_ne!(pxz.z.is_equal_ct(&z2), 0);

    let mut t = Fp2::default();
    let mut ao3 = Fp2::default();
    let mut idw = JacPoint::default();
    jac_to_ws(&mut idw, &mut t, &mut ao3, &id, &e);
    let mut dbl_id = JacPoint::default();
    let mut u = Fp2::default();
    jac_dblw(&mut dbl_id, &mut u, &idw, &t);
    assert_ne!(dbl_id.z.is_zero_ct(), 0);
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
    assert_eq!(img[0].is_zero_ct(), 0);
    let mut img2 = [EcPoint::default()];
    xeval_2_singular(&mut img2, &[basis.p], &kps);
    assert_ne!(img[0].is_equal_ct(&img2[0]), 0);
}

#[test]
#[cfg(debug_assertions)]
fn basis_from_bad_hint_rejected() {
    let mut e = EcCurve::e0();
    e.a = Fp2::from_small(6);
    e.normalize_a24();
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
    a24_unnorm.x = ax * three;
    a24_unnorm.z = az * three;

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
    assert_ne!(r1.is_equal_ct(&r2), 0);
    assert_ne!(s1.is_equal_ct(&s2), 0);
}
