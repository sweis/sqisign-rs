// SPDX-License-Identifier: Apache-2.0
//! (2ⁿ, 2ⁿ)-isogeny chain in the theta model: gluing, generic steps, splitting.
//!
//! Ported from `src/hd/ref/lvlx/theta_isogenies.c`.

use super::*;
use crate::ec::{
    ec_curve_init, jac_to_xz_add_components, lift_basis, test_jac_order_twof, EcBasis,
};
use crate::gf::{fp2_batched_inv, fp2_sqrt};
use crate::precomp::{
    PrecompBasisChangeMatrix, CHI_EVAL, EVEN_INDEX, FP2_CONSTANTS, NORMALIZATION_TRANSFORMS,
    SPLITTING_TRANSFORMS,
};

// ===========================================================================
// Basis-change matrix helpers
// ===========================================================================

#[inline]
fn select_base_change_matrix(
    m: &mut BasisChangeMatrix,
    m1: &BasisChangeMatrix,
    m2: &PrecompBasisChangeMatrix,
    option: u32,
) {
    for i in 0..4 {
        for j in 0..4 {
            let a = m1.m[i][j];
            fp2_select(
                &mut m.m[i][j],
                &a,
                &FP2_CONSTANTS[m2.m[i][j] as usize],
                option,
            );
        }
    }
}

#[inline]
fn set_base_change_matrix_from_precomp(res: &mut BasisChangeMatrix, m: &PrecompBasisChangeMatrix) {
    for i in 0..4 {
        for j in 0..4 {
            res.m[i][j] = FP2_CONSTANTS[m.m[i][j] as usize];
        }
    }
}

#[inline]
fn choose_index_theta_point(res: &mut Fp2, ind: i32, t: &ThetaPoint) {
    let src = match ind & 3 {
        0 => &t.x,
        1 => &t.y,
        2 => &t.z,
        3 => &t.t,
        _ => unreachable!(),
    };
    fp2_copy(res, src);
}

/// res ← M · P, optionally skipping the t-column when `pt_not_zero == false`.
fn apply_isomorphism_general(
    res: &mut ThetaPoint,
    m: &BasisChangeMatrix,
    p: &ThetaPoint,
    pt_not_zero: bool,
) {
    let mut x1 = Fp2::default();
    let mut temp = ThetaPoint::default();

    fp2_mul(&mut temp.x, &p.x, &m.m[0][0]);
    fp2_mul(&mut x1, &p.y, &m.m[0][1]);
    fp2_add_ip(&mut temp.x, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[0][2]);
    fp2_add_ip(&mut temp.x, &x1);

    fp2_mul(&mut temp.y, &p.x, &m.m[1][0]);
    fp2_mul(&mut x1, &p.y, &m.m[1][1]);
    fp2_add_ip(&mut temp.y, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[1][2]);
    fp2_add_ip(&mut temp.y, &x1);

    fp2_mul(&mut temp.z, &p.x, &m.m[2][0]);
    fp2_mul(&mut x1, &p.y, &m.m[2][1]);
    fp2_add_ip(&mut temp.z, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[2][2]);
    fp2_add_ip(&mut temp.z, &x1);

    fp2_mul(&mut temp.t, &p.x, &m.m[3][0]);
    fp2_mul(&mut x1, &p.y, &m.m[3][1]);
    fp2_add_ip(&mut temp.t, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[3][2]);
    fp2_add_ip(&mut temp.t, &x1);

    if pt_not_zero {
        fp2_mul(&mut x1, &p.t, &m.m[0][3]);
        fp2_add_ip(&mut temp.x, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[1][3]);
        fp2_add_ip(&mut temp.y, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[2][3]);
        fp2_add_ip(&mut temp.z, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[3][3]);
        fp2_add_ip(&mut temp.t, &x1);
    }

    fp2_copy(&mut res.x, &temp.x);
    fp2_copy(&mut res.y, &temp.y);
    fp2_copy(&mut res.z, &temp.z);
    fp2_copy(&mut res.t, &temp.t);
}

#[inline]
fn apply_isomorphism(res: &mut ThetaPoint, m: &BasisChangeMatrix, p: &ThetaPoint) {
    apply_isomorphism_general(res, m, p, true);
}

/// res ← M1 · M2 (4×4 matrix product).
fn base_change_matrix_multiplication(
    res: &mut BasisChangeMatrix,
    m1: &BasisChangeMatrix,
    m2: &BasisChangeMatrix,
) {
    let mut tmp = BasisChangeMatrix::default();
    let mut sum = Fp2::default();
    let mut m_ik;
    for i in 0..4 {
        for j in 0..4 {
            fp2_set_zero(&mut sum);
            for k in 0..4 {
                m_ik = m1.m[i][k];
                let m_kj = m2.m[k][j];
                fp2_mul_ip(&mut m_ik, &m_kj);
                let ss = sum;
                fp2_add(&mut sum, &ss, &m_ik);
            }
            tmp.m[i][j] = sum;
        }
    }
    *res = tmp;
}

// ===========================================================================
// Gluing
// ===========================================================================

/// Map a couple point on E₁×E₂ to its theta point under the gluing's basis change.
fn base_change(out: &mut ThetaPoint, phi: &ThetaGluing, t: &ThetaCouplePoint) {
    let mut null_point = ThetaPoint::default();
    fp2_mul(&mut null_point.x, &t.p1.x, &t.p2.x);
    fp2_mul(&mut null_point.y, &t.p1.x, &t.p2.z);
    fp2_mul(&mut null_point.z, &t.p2.x, &t.p1.z);
    fp2_mul(&mut null_point.t, &t.p1.z, &t.p2.z);
    apply_isomorphism(out, &phi.m, &null_point);
}

fn action_by_translation_z_and_det(z_inv: &mut Fp2, det_inv: &mut Fp2, p4: &EcPoint, p2: &EcPoint) {
    fp2_copy(z_inv, &p4.z);
    let mut tmp = Fp2::default();
    fp2_mul(det_inv, &p4.x, &p2.z);
    fp2_mul(&mut tmp, &p4.z, &p2.x);
    fp2_sub_ip(det_inv, &tmp);
}

fn action_by_translation_compute_matrix(
    g: &mut TranslationMatrix,
    p4: &EcPoint,
    p2: &EcPoint,
    z_inv: &Fp2,
    det_inv: &Fp2,
) {
    let mut tmp = Fp2::default();

    fp2_mul(&mut tmp, &p4.x, z_inv);
    fp2_mul(&mut g.g10, &p4.x, &p2.x);
    fp2_mul_ip(&mut g.g10, det_inv);
    fp2_sub_ip(&mut g.g10, &tmp);

    fp2_mul(&mut g.g11, &p2.x, det_inv);
    fp2_mul_ip(&mut g.g11, &p4.z);

    let s = g.g11;
    fp2_neg(&mut g.g00, &s);

    fp2_mul(&mut g.g01, &p2.z, det_inv);
    fp2_mul_ip(&mut g.g01, &p4.z);
    fp2_neg_ip(&mut g.g01);
}

fn verify_two_torsion(
    k1_2: &ThetaCouplePoint,
    k2_2: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> bool {
    if (ec_is_zero(&k1_2.p1) | ec_is_zero(&k1_2.p2) | ec_is_zero(&k2_2.p1) | ec_is_zero(&k2_2.p2))
        != 0
    {
        return false;
    }
    if (ec_is_equal(&k1_2.p1, &k2_2.p1) | ec_is_equal(&k1_2.p2, &k2_2.p2)) != 0 {
        return false;
    }
    let mut o1 = ThetaCouplePoint::default();
    let mut o2 = ThetaCouplePoint::default();
    double_couple_point(&mut o1, k1_2, e12);
    double_couple_point(&mut o2, k2_2, e12);
    if (ec_is_zero(&o1.p1) & ec_is_zero(&o1.p2) & ec_is_zero(&o2.p1) & ec_is_zero(&o2.p2)) == 0 {
        return false;
    }
    true
}

fn action_by_translation(
    gi: &mut [TranslationMatrix; 4],
    k1_4: &ThetaCouplePoint,
    k2_4: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> bool {
    let mut k1_2 = ThetaCouplePoint::default();
    let mut k2_2 = ThetaCouplePoint::default();
    double_couple_point(&mut k1_2, k1_4, e12);
    double_couple_point(&mut k2_2, k2_4, e12);

    if !verify_two_torsion(&k1_2, &k2_2, e12) {
        return false;
    }

    let mut inv = [Fp2::default(); 8];
    {
        let (z, d) = inv.split_at_mut(4);
        action_by_translation_z_and_det(&mut z[0], &mut d[0], &k1_4.p1, &k1_2.p1);
        action_by_translation_z_and_det(&mut z[1], &mut d[1], &k1_4.p2, &k1_2.p2);
        action_by_translation_z_and_det(&mut z[2], &mut d[2], &k2_4.p1, &k2_2.p1);
        action_by_translation_z_and_det(&mut z[3], &mut d[3], &k2_4.p2, &k2_2.p2);
    }

    fp2_batched_inv(&mut inv);
    if fp2_is_zero(&inv[0]) != 0 {
        return false;
    }

    action_by_translation_compute_matrix(&mut gi[0], &k1_4.p1, &k1_2.p1, &inv[0], &inv[4]);
    action_by_translation_compute_matrix(&mut gi[1], &k1_4.p2, &k1_2.p2, &inv[1], &inv[5]);
    action_by_translation_compute_matrix(&mut gi[2], &k2_4.p1, &k2_2.p1, &inv[2], &inv[6]);
    action_by_translation_compute_matrix(&mut gi[3], &k2_4.p2, &k2_2.p2, &inv[3], &inv[7]);

    true
}

fn gluing_change_of_basis(
    m: &mut BasisChangeMatrix,
    k1_4: &ThetaCouplePoint,
    k2_4: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> bool {
    let mut gi = [TranslationMatrix::default(); 4];
    if !action_by_translation(&mut gi, k1_4, k2_4, e12) {
        return false;
    }

    let mut t001 = Fp2::default();
    let mut t101 = Fp2::default();
    let mut t002 = Fp2::default();
    let mut t102 = Fp2::default();
    let mut tmp = Fp2::default();

    fp2_mul(&mut t001, &gi[0].g00, &gi[2].g00);
    fp2_mul(&mut tmp, &gi[0].g01, &gi[2].g10);
    fp2_add_ip(&mut t001, &tmp);

    fp2_mul(&mut t101, &gi[0].g10, &gi[2].g00);
    fp2_mul(&mut tmp, &gi[0].g11, &gi[2].g10);
    fp2_add_ip(&mut t101, &tmp);

    fp2_mul(&mut t002, &gi[1].g00, &gi[3].g00);
    fp2_mul(&mut tmp, &gi[1].g01, &gi[3].g10);
    fp2_add_ip(&mut t002, &tmp);

    fp2_mul(&mut t102, &gi[1].g10, &gi[3].g00);
    fp2_mul(&mut tmp, &gi[1].g11, &gi[3].g10);
    fp2_add_ip(&mut t102, &tmp);

    // Row 0: trace.
    fp2_set_one(&mut m.m[0][0]);
    fp2_mul(&mut tmp, &t001, &t002);
    fp2_add_ip(&mut m.m[0][0], &tmp);
    fp2_mul(&mut tmp, &gi[2].g00, &gi[3].g00);
    fp2_add_ip(&mut m.m[0][0], &tmp);
    fp2_mul(&mut tmp, &gi[0].g00, &gi[1].g00);
    fp2_add_ip(&mut m.m[0][0], &tmp);

    fp2_mul(&mut m.m[0][1], &t001, &t102);
    fp2_mul(&mut tmp, &gi[2].g00, &gi[3].g10);
    fp2_add_ip(&mut m.m[0][1], &tmp);
    fp2_mul(&mut tmp, &gi[0].g00, &gi[1].g10);
    fp2_add_ip(&mut m.m[0][1], &tmp);

    fp2_mul(&mut m.m[0][2], &t101, &t002);
    fp2_mul(&mut tmp, &gi[2].g10, &gi[3].g00);
    fp2_add_ip(&mut m.m[0][2], &tmp);
    fp2_mul(&mut tmp, &gi[0].g10, &gi[1].g00);
    fp2_add_ip(&mut m.m[0][2], &tmp);

    fp2_mul(&mut m.m[0][3], &t101, &t102);
    fp2_mul(&mut tmp, &gi[2].g10, &gi[3].g10);
    fp2_add_ip(&mut m.m[0][3], &tmp);
    fp2_mul(&mut tmp, &gi[0].g10, &gi[1].g10);
    fp2_add_ip(&mut m.m[0][3], &tmp);

    // Row 1: action of (0, K2_4.P2).
    let r0 = m.m[0];
    fp2_mul(&mut tmp, &gi[3].g01, &r0[1]);
    fp2_mul(&mut m.m[1][0], &gi[3].g00, &r0[0]);
    fp2_add_ip(&mut m.m[1][0], &tmp);

    fp2_mul(&mut tmp, &gi[3].g11, &r0[1]);
    fp2_mul(&mut m.m[1][1], &gi[3].g10, &r0[0]);
    fp2_add_ip(&mut m.m[1][1], &tmp);

    fp2_mul(&mut tmp, &gi[3].g01, &r0[3]);
    fp2_mul(&mut m.m[1][2], &gi[3].g00, &r0[2]);
    fp2_add_ip(&mut m.m[1][2], &tmp);

    fp2_mul(&mut tmp, &gi[3].g11, &r0[3]);
    fp2_mul(&mut m.m[1][3], &gi[3].g10, &r0[2]);
    fp2_add_ip(&mut m.m[1][3], &tmp);

    // Row 2: action of (K1_4.P1, 0).
    fp2_mul(&mut tmp, &gi[0].g01, &r0[2]);
    fp2_mul(&mut m.m[2][0], &gi[0].g00, &r0[0]);
    fp2_add_ip(&mut m.m[2][0], &tmp);

    fp2_mul(&mut tmp, &gi[0].g01, &r0[3]);
    fp2_mul(&mut m.m[2][1], &gi[0].g00, &r0[1]);
    fp2_add_ip(&mut m.m[2][1], &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &r0[2]);
    fp2_mul(&mut m.m[2][2], &gi[0].g10, &r0[0]);
    fp2_add_ip(&mut m.m[2][2], &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &r0[3]);
    fp2_mul(&mut m.m[2][3], &gi[0].g10, &r0[1]);
    fp2_add_ip(&mut m.m[2][3], &tmp);

    // Row 3: action of (K1_4.P1, K2_4.P2).
    let r1 = m.m[1];
    fp2_mul(&mut tmp, &gi[0].g01, &r1[2]);
    fp2_mul(&mut m.m[3][0], &gi[0].g00, &r1[0]);
    fp2_add_ip(&mut m.m[3][0], &tmp);

    fp2_mul(&mut tmp, &gi[0].g01, &r1[3]);
    fp2_mul(&mut m.m[3][1], &gi[0].g00, &r1[1]);
    fp2_add_ip(&mut m.m[3][1], &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &r1[2]);
    fp2_mul(&mut m.m[3][2], &gi[0].g10, &r1[0]);
    fp2_add_ip(&mut m.m[3][2], &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &r1[3]);
    fp2_mul(&mut m.m[3][3], &gi[0].g10, &r1[1]);
    fp2_add_ip(&mut m.m[3][3], &tmp);

    true
}

fn gluing_compute(
    out: &mut ThetaGluing,
    e12: &ThetaCoupleCurve,
    xy_k1_8: &ThetaCoupleJacPoint,
    xy_k2_8: &ThetaCoupleJacPoint,
    verify: bool,
) -> bool {
    // On the verify path the kernel is derived from attacker-controlled
    // signature data, so a wrong-order point means "reject", not "panic".
    // The C reference does `#ifndef NDEBUG assert(...)` here (debug crash).
    if verify {
        if test_jac_order_twof(&xy_k1_8.p1, &e12.e1, 3) == 0
            || test_jac_order_twof(&xy_k2_8.p1, &e12.e1, 3) == 0
            || test_jac_order_twof(&xy_k1_8.p2, &e12.e2, 3) == 0
            || test_jac_order_twof(&xy_k2_8.p2, &e12.e2, 3) == 0
        {
            return false;
        }
    } else {
        debug_assert!(test_jac_order_twof(&xy_k1_8.p1, &e12.e1, 3) != 0);
        debug_assert!(test_jac_order_twof(&xy_k2_8.p1, &e12.e1, 3) != 0);
        debug_assert!(test_jac_order_twof(&xy_k1_8.p2, &e12.e2, 3) != 0);
        debug_assert!(test_jac_order_twof(&xy_k2_8.p2, &e12.e2, 3) != 0);
    }

    out.xy_k1_8 = *xy_k1_8;
    out.domain = *e12;

    let mut xy_k1_4 = ThetaCoupleJacPoint::default();
    let mut xy_k2_4 = ThetaCoupleJacPoint::default();
    double_couple_jac_point(&mut xy_k1_4, xy_k1_8, e12);
    double_couple_jac_point(&mut xy_k2_4, xy_k2_8, e12);

    let mut k1_8 = ThetaCouplePoint::default();
    let mut k2_8 = ThetaCouplePoint::default();
    let mut k1_4 = ThetaCouplePoint::default();
    let mut k2_4 = ThetaCouplePoint::default();
    couple_jac_to_xz(&mut k1_8, xy_k1_8);
    couple_jac_to_xz(&mut k2_8, xy_k2_8);
    couple_jac_to_xz(&mut k1_4, &xy_k1_4);
    couple_jac_to_xz(&mut k2_4, &xy_k2_4);

    if !gluing_change_of_basis(&mut out.m, &k1_4, &k2_4, e12) {
        return false;
    }

    let mut tt1 = ThetaPoint::default();
    let mut tt2 = ThetaPoint::default();
    base_change(&mut tt1, out, &k1_8);
    base_change(&mut tt2, out, &k2_8);

    to_squared_theta_ip(&mut tt1);
    to_squared_theta_ip(&mut tt2);

    if (fp2_is_zero(&tt1.t) & fp2_is_zero(&tt2.t)) == 0 {
        return false;
    }
    if (fp2_is_zero(&tt1.x)
        | fp2_is_zero(&tt2.x)
        | fp2_is_zero(&tt1.y)
        | fp2_is_zero(&tt2.z)
        | fp2_is_zero(&tt1.z))
        != 0
    {
        return false;
    }

    fp2_mul(&mut out.codomain.x, &tt1.x, &tt2.x);
    fp2_mul(&mut out.codomain.y, &tt1.y, &tt2.x);
    fp2_mul(&mut out.codomain.z, &tt1.x, &tt2.z);
    fp2_set_zero(&mut out.codomain.t);

    fp2_mul(&mut out.precomputation.x, &tt1.y, &tt2.z);
    fp2_copy(&mut out.precomputation.y, &out.codomain.z);
    fp2_copy(&mut out.precomputation.z, &out.codomain.y);
    fp2_set_zero(&mut out.precomputation.t);

    fp2_mul(&mut out.image_k1_8.x, &tt1.x, &out.precomputation.x);
    fp2_mul(&mut out.image_k1_8.y, &tt1.z, &out.precomputation.z);

    if verify {
        let mut t1 = Fp2::default();
        let mut t2 = Fp2::default();
        fp2_mul(&mut t1, &tt1.y, &out.precomputation.y);
        if fp2_is_equal(&out.image_k1_8.x, &t1) == 0 {
            return false;
        }
        fp2_mul(&mut t1, &tt2.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt2.z, &out.precomputation.z);
        if fp2_is_equal(&t2, &t1) == 0 {
            return false;
        }
    }

    hadamard_ip(&mut out.codomain);
    true
}

fn gluing_eval_point(image: &mut ThetaPoint, p: &ThetaCoupleJacPoint, phi: &ThetaGluing) {
    let mut t1 = ThetaPoint::default();
    let mut t2 = ThetaPoint::default();
    let mut ac1 = AddComponents::default();
    let mut ac2 = AddComponents::default();

    jac_to_xz_add_components(&mut ac1, &p.p1, &phi.xy_k1_8.p1, &phi.domain.e1);
    jac_to_xz_add_components(&mut ac2, &p.p2, &phi.xy_k1_8.p2, &phi.domain.e2);

    fp2_mul(&mut t1.x, &ac1.u, &ac2.u);
    fp2_mul(&mut t2.t, &ac1.v, &ac2.v);
    fp2_add_ip(&mut t1.x, &t2.t);
    fp2_mul(&mut t1.y, &ac1.u, &ac2.w);
    fp2_mul(&mut t1.z, &ac1.w, &ac2.u);
    fp2_mul(&mut t1.t, &ac1.w, &ac2.w);
    fp2_add(&mut t2.x, &ac1.u, &ac1.v);
    fp2_add(&mut t2.y, &ac2.u, &ac2.v);
    let (sx, sy) = (t2.x, t2.y);
    fp2_mul(&mut t2.x, &sx, &sy);
    let (s2x, s1x) = (t2.x, t1.x);
    fp2_sub(&mut t2.x, &s2x, &s1x);
    fp2_mul(&mut t2.y, &ac1.v, &ac2.w);
    fp2_mul(&mut t2.z, &ac1.w, &ac2.v);
    fp2_set_zero(&mut t2.t);

    let s = t1;
    apply_isomorphism_general(&mut t1, &phi.m, &s, true);
    let s = t2;
    apply_isomorphism_general(&mut t2, &phi.m, &s, false);
    pointwise_square_ip(&mut t1);
    pointwise_square_ip(&mut t2);

    let (a, b) = (t1.x, t2.x);
    fp2_sub(&mut t1.x, &a, &b);
    let (a, b) = (t1.y, t2.y);
    fp2_sub(&mut t1.y, &a, &b);
    let (a, b) = (t1.z, t2.z);
    fp2_sub(&mut t1.z, &a, &b);
    let (a, b) = (t1.t, t2.t);
    fp2_sub(&mut t1.t, &a, &b);
    hadamard_ip(&mut t1);

    fp2_mul(&mut image.x, &t1.x, &phi.image_k1_8.y);
    fp2_mul(&mut image.y, &t1.y, &phi.image_k1_8.y);
    fp2_mul(&mut image.z, &t1.z, &phi.image_k1_8.x);
    fp2_mul(&mut image.t, &t1.t, &phi.image_k1_8.x);

    hadamard_ip(image);
}

fn gluing_eval_point_special_case(
    image: &mut ThetaPoint,
    p: &ThetaCouplePoint,
    phi: &ThetaGluing,
) -> bool {
    let mut t = ThetaPoint::default();
    base_change(&mut t, phi, p);
    to_squared_theta_ip(&mut t);

    if fp2_is_zero(&t.t) == 0 {
        return false;
    }

    fp2_mul(&mut image.x, &t.x, &phi.precomputation.x);
    fp2_mul(&mut image.y, &t.y, &phi.precomputation.y);
    fp2_mul(&mut image.z, &t.z, &phi.precomputation.z);
    fp2_set_zero(&mut image.t);

    hadamard_ip(image);
    true
}

fn gluing_eval_basis(
    image1: &mut ThetaPoint,
    image2: &mut ThetaPoint,
    xy_t1: &ThetaCoupleJacPoint,
    xy_t2: &ThetaCoupleJacPoint,
    phi: &ThetaGluing,
) {
    gluing_eval_point(image1, xy_t1, phi);
    gluing_eval_point(image2, xy_t2, phi);
}

// ===========================================================================
// Generic (2,2)-isogeny step
// ===========================================================================

fn theta_isogeny_compute(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_8: &ThetaPoint,
    t2_8: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
    verify: bool,
) -> bool {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_8;
    out.t2_8 = *t2_8;
    out.codomain.precomputation = false;

    let mut tt1 = ThetaPoint::default();
    let mut tt2 = ThetaPoint::default();

    if hadamard_bool_1 {
        hadamard(&mut tt1, t1_8);
        to_squared_theta_ip(&mut tt1);
        hadamard(&mut tt2, t2_8);
        to_squared_theta_ip(&mut tt2);
    } else {
        to_squared_theta(&mut tt1, t1_8);
        to_squared_theta(&mut tt2, t2_8);
    }

    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();

    if (fp2_is_zero(&tt2.x)
        | fp2_is_zero(&tt2.y)
        | fp2_is_zero(&tt2.z)
        | fp2_is_zero(&tt2.t)
        | fp2_is_zero(&tt1.x)
        | fp2_is_zero(&tt1.y))
        != 0
    {
        return false;
    }

    fp2_mul(&mut t1, &tt1.x, &tt2.y);
    fp2_mul(&mut t2, &tt1.y, &tt2.x);
    fp2_mul(&mut out.codomain.null_point.x, &tt2.x, &t1);
    fp2_mul(&mut out.codomain.null_point.y, &tt2.y, &t2);
    fp2_mul(&mut out.codomain.null_point.z, &tt2.z, &t1);
    fp2_mul(&mut out.codomain.null_point.t, &tt2.t, &t2);
    let mut t3 = Fp2::default();
    fp2_mul(&mut t3, &tt2.z, &tt2.t);
    fp2_mul(&mut out.precomputation.x, &t3, &tt1.y);
    fp2_mul(&mut out.precomputation.y, &t3, &tt1.x);
    fp2_copy(&mut out.precomputation.z, &out.codomain.null_point.t);
    fp2_copy(&mut out.precomputation.t, &out.codomain.null_point.z);

    if verify {
        fp2_mul(&mut t1, &tt1.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt1.y, &out.precomputation.y);
        if fp2_is_equal(&t1, &t2) == 0 {
            return false;
        }
        fp2_mul(&mut t1, &tt1.z, &out.precomputation.z);
        fp2_mul(&mut t2, &tt1.t, &out.precomputation.t);
        if fp2_is_equal(&t1, &t2) == 0 {
            return false;
        }
        fp2_mul(&mut t1, &tt2.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt2.z, &out.precomputation.z);
        if fp2_is_equal(&t1, &t2) == 0 {
            return false;
        }
        fp2_mul(&mut t1, &tt2.y, &out.precomputation.y);
        fp2_mul(&mut t2, &tt2.t, &out.precomputation.t);
        if fp2_is_equal(&t1, &t2) == 0 {
            return false;
        }
    }

    if hadamard_bool_2 {
        hadamard_ip(&mut out.codomain.null_point);
    }
    true
}

fn theta_isogeny_compute_4(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_4: &ThetaPoint,
    t2_4: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_4;
    out.t2_8 = *t2_4;
    out.codomain.precomputation = false;

    let mut tt1 = ThetaPoint::default();
    let mut tt2 = ThetaPoint::default();

    if hadamard_bool_1 {
        hadamard(&mut tt1, t1_4);
        to_squared_theta_ip(&mut tt1);
        hadamard(&mut tt2, &a.null_point);
        to_squared_theta_ip(&mut tt2);
    } else {
        to_squared_theta(&mut tt1, t1_4);
        to_squared_theta(&mut tt2, &a.null_point);
    }

    let mut sqaabb = Fp2::default();
    let mut sqaacc = Fp2::default();
    fp2_mul(&mut sqaabb, &tt2.x, &tt2.y);
    fp2_mul(&mut sqaacc, &tt2.x, &tt2.z);
    fp2_sqrt(&mut sqaabb);
    fp2_sqrt(&mut sqaacc);

    fp2_mul(&mut out.codomain.null_point.y, &sqaabb, &sqaacc);
    let s = out.codomain.null_point.y;
    fp2_mul(&mut out.precomputation.t, &s, &tt1.z);
    fp2_mul_ip(&mut out.codomain.null_point.y, &tt1.x);

    fp2_mul(&mut out.codomain.null_point.t, &tt1.z, &sqaabb);
    fp2_mul_ip(&mut out.codomain.null_point.t, &tt2.x);

    fp2_mul(&mut out.codomain.null_point.x, &tt1.x, &tt2.x);
    let s = out.codomain.null_point.x;
    fp2_mul(&mut out.codomain.null_point.z, &s, &tt2.z);
    fp2_mul_ip(&mut out.codomain.null_point.x, &sqaacc);

    fp2_mul(&mut out.precomputation.x, &tt1.x, &tt2.t);
    let s = out.precomputation.x;
    fp2_mul(&mut out.precomputation.z, &s, &tt2.y);
    fp2_mul_ip(&mut out.precomputation.x, &tt2.z);
    let s = out.precomputation.x;
    fp2_mul(&mut out.precomputation.y, &s, &sqaabb);
    fp2_mul_ip(&mut out.precomputation.x, &tt2.y);
    fp2_mul_ip(&mut out.precomputation.z, &sqaacc);
    fp2_mul_ip(&mut out.precomputation.t, &tt2.y);

    if hadamard_bool_2 {
        hadamard_ip(&mut out.codomain.null_point);
    }
}

fn theta_isogeny_compute_2(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_2: &ThetaPoint,
    t2_2: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_2;
    out.t2_8 = *t2_2;
    out.codomain.precomputation = false;

    let mut tt2 = ThetaPoint::default();
    if hadamard_bool_1 {
        hadamard(&mut tt2, &a.null_point);
        to_squared_theta_ip(&mut tt2);
    } else {
        to_squared_theta(&mut tt2, &a.null_point);
    }

    fp2_copy(&mut out.codomain.null_point.x, &tt2.x);
    fp2_mul(&mut out.codomain.null_point.y, &tt2.x, &tt2.y);
    fp2_mul(&mut out.codomain.null_point.z, &tt2.x, &tt2.z);
    fp2_mul(&mut out.codomain.null_point.t, &tt2.x, &tt2.t);
    fp2_sqrt(&mut out.codomain.null_point.y);
    fp2_sqrt(&mut out.codomain.null_point.z);
    fp2_sqrt(&mut out.codomain.null_point.t);

    fp2_mul(&mut out.precomputation.x, &tt2.z, &tt2.t);
    let s = out.precomputation.x;
    fp2_mul(&mut out.precomputation.y, &s, &out.codomain.null_point.y);
    fp2_mul_ip(&mut out.precomputation.x, &tt2.y);
    fp2_mul(
        &mut out.precomputation.z,
        &tt2.t,
        &out.codomain.null_point.z,
    );
    fp2_mul_ip(&mut out.precomputation.z, &tt2.y);
    fp2_mul(
        &mut out.precomputation.t,
        &tt2.z,
        &out.codomain.null_point.t,
    );
    fp2_mul_ip(&mut out.precomputation.t, &tt2.y);

    if hadamard_bool_2 {
        hadamard_ip(&mut out.codomain.null_point);
    }
}

fn theta_isogeny_eval(out: &mut ThetaPoint, phi: &ThetaIsogeny, p: &ThetaPoint) {
    if phi.hadamard_bool_1 {
        hadamard(out, p);
        to_squared_theta_ip(out);
    } else {
        to_squared_theta(out, p);
    }
    fp2_mul_ip(&mut out.x, &phi.precomputation.x);
    fp2_mul_ip(&mut out.y, &phi.precomputation.y);
    fp2_mul_ip(&mut out.z, &phi.precomputation.z);
    fp2_mul_ip(&mut out.t, &phi.precomputation.t);

    if phi.hadamard_bool_2 {
        hadamard_ip(out);
    }
}

// ===========================================================================
// Splitting
// ===========================================================================

#[cfg(feature = "sign")]
fn sample_random_index() -> u8 {
    use crate::common::ctrdrbg::randombytes;
    let mut seed_arr = [0u8; 4];
    let mut seed: u32;
    loop {
        randombytes(&mut seed_arr);
        seed = u32::from_le_bytes(seed_arr);
        if seed < 4_294_967_292 {
            break;
        }
    }
    let secret_index = seed - (((seed as u64) * 2_863_311_531u64 >> 34) as u32) * 6;
    debug_assert_eq!(secret_index, seed % 6);
    secret_index as u8
}

fn splitting_compute(
    out: &mut ThetaSplitting,
    a: &ThetaStructure,
    zero_index: i32,
    randomize: bool,
) -> bool {
    let mut count: u32 = 0;
    let mut u_cst = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();

    out.m = BasisChangeMatrix::default();

    for i in 0..10 {
        fp2_set_zero(&mut u_cst);
        let chi_row = &CHI_EVAL[EVEN_INDEX[i][0] as usize];
        for (t, &chi) in chi_row.iter().enumerate() {
            choose_index_theta_point(&mut t2, t as i32, &a.null_point);
            choose_index_theta_point(&mut t1, (t as i32) ^ EVEN_INDEX[i][1], &a.null_point);
            fp2_mul_ip(&mut t1, &t2);

            let ctl = (chi >> 1) as u32;
            debug_assert!(ctl == 0 || ctl == 0xFFFF_FFFF);

            let s = t1;
            fp2_neg(&mut t2, &s);
            let s = t1;
            fp2_select(&mut t1, &s, &t2, ctl);

            fp2_add_ip(&mut u_cst, &t1);
        }

        let ctl = fp2_is_zero(&u_cst);
        count = count.wrapping_sub(ctl);
        let m_prev = out.m;
        select_base_change_matrix(&mut out.m, &m_prev, &SPLITTING_TRANSFORMS[i], ctl);
        if zero_index != -1 && i as i32 == zero_index && ctl == 0 {
            return false;
        }
    }

    #[cfg(feature = "sign")]
    if randomize {
        let secret_index = sample_random_index();
        let mut m_random = BasisChangeMatrix::default();
        set_base_change_matrix_from_precomp(&mut m_random, &NORMALIZATION_TRANSFORMS[0]);
        for i in 1u8..6 {
            let mut mask = i as i32 - secret_index as i32;
            mask = (mask | mask.wrapping_neg()) >> 31;
            let prev = m_random;
            select_base_change_matrix(
                &mut m_random,
                &prev,
                &NORMALIZATION_TRANSFORMS[i as usize],
                !(mask as u32),
            );
        }
        let prev = out.m;
        base_change_matrix_multiplication(&mut out.m, &m_random, &prev);
    }
    #[cfg(not(feature = "sign"))]
    {
        let _ = randomize;
        debug_assert!(!randomize);
    }

    apply_isomorphism(&mut out.b.null_point, &out.m, &a.null_point);

    count == 1
}

fn theta_product_structure_to_elliptic_product(
    e12: &mut ThetaCoupleCurve,
    a: &ThetaStructure,
) -> bool {
    let mut xx = Fp2::default();
    let mut yy = Fp2::default();

    if is_product_theta_point(&a.null_point) == 0 {
        return false;
    }

    ec_curve_init(&mut e12.e1);
    ec_curve_init(&mut e12.e2);

    if (fp2_is_zero(&a.null_point.x) | fp2_is_zero(&a.null_point.y) | fp2_is_zero(&a.null_point.z))
        != 0
    {
        return false;
    }

    fp2_sqr(&mut xx, &a.null_point.x);
    fp2_sqr(&mut yy, &a.null_point.y);
    fp2_sqr_ip(&mut xx);
    fp2_sqr_ip(&mut yy);

    fp2_add(&mut e12.e2.a, &xx, &yy);
    fp2_sub(&mut e12.e2.c, &xx, &yy);
    fp2_dbl_ip(&mut e12.e2.a);
    fp2_neg_ip(&mut e12.e2.a);

    fp2_sqr(&mut xx, &a.null_point.x);
    fp2_sqr(&mut yy, &a.null_point.z);
    fp2_sqr_ip(&mut xx);
    fp2_sqr_ip(&mut yy);

    fp2_add(&mut e12.e1.a, &xx, &yy);
    fp2_sub(&mut e12.e1.c, &xx, &yy);
    fp2_dbl_ip(&mut e12.e1.a);
    fp2_neg_ip(&mut e12.e1.a);

    if (fp2_is_zero(&e12.e1.c) | fp2_is_zero(&e12.e2.c)) != 0 {
        return false;
    }
    true
}

fn theta_point_to_montgomery_point(
    p12: &mut ThetaCouplePoint,
    p: &ThetaPoint,
    a: &ThetaStructure,
) -> bool {
    let mut temp = Fp2::default();

    if is_product_theta_point(p) == 0 {
        return false;
    }

    let (mut x, mut z) = (&p.x, &p.y);
    if (fp2_is_zero(x) & fp2_is_zero(z)) != 0 {
        x = &p.z;
        z = &p.t;
    }
    if (fp2_is_zero(x) & fp2_is_zero(z)) != 0 {
        return false;
    }
    fp2_mul(&mut p12.p2.x, &a.null_point.y, x);
    fp2_mul(&mut temp, &a.null_point.x, z);
    let s = p12.p2.x;
    fp2_sub(&mut p12.p2.z, &temp, &s);
    fp2_add_ip(&mut p12.p2.x, &temp);

    let (mut x, mut z) = (&p.x, &p.z);
    if (fp2_is_zero(x) & fp2_is_zero(z)) != 0 {
        x = &p.y;
        z = &p.t;
    }
    fp2_mul(&mut p12.p1.x, &a.null_point.z, x);
    fp2_mul(&mut temp, &a.null_point.x, z);
    let s = p12.p1.x;
    fp2_sub(&mut p12.p1.z, &temp, &s);
    fp2_add_ip(&mut p12.p1.x, &temp);
    true
}

// ===========================================================================
// Chain driver
// ===========================================================================

#[allow(clippy::too_many_arguments)]
fn theta_chain_compute_impl(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    verify: bool,
    randomize: bool,
) -> i32 {
    // The strategy initialisation at `todo[0] = (n - 2 + extra)` requires this;
    // verify-path callers already guarantee it via the pow_dim2_deg_resp gate.
    debug_assert!(n >= 2 || extra_torsion);
    let num_p = p12.len();

    let mut theta = ThetaStructure::default();

    let mut xy_t1 = ThetaCoupleJacPoint::default();
    let mut xy_t2 = ThetaCoupleJacPoint::default();

    let mut bas1 = EcBasis {
        p: ker.t1.p1,
        q: ker.t2.p1,
        pmq: ker.t1m2.p1,
    };
    let mut bas2 = EcBasis {
        p: ker.t1.p2,
        q: ker.t2.p2,
        pmq: ker.t1m2.p2,
    };
    if lift_basis(&mut xy_t1.p1, &mut xy_t2.p1, &mut bas1, &mut e12.e1) == 0 {
        return 0;
    }
    if lift_basis(&mut xy_t1.p2, &mut xy_t2.p2, &mut bas2, &mut e12.e2) == 0 {
        return 0;
    }

    let extra: u32 = if extra_torsion { HD_EXTRA_TORSION } else { 0 };

    debug_assert!(extra == 0 || extra == 2);
    if verify {
        if test_point_order_twof(&bas2.p, &e12.e2, (n + extra) as i32) == 0
            || test_jac_order_twof(&xy_t2.p2, &e12.e2, (n + extra) as i32) == 0
        {
            return 0;
        }
    } else {
        debug_assert!(test_point_order_twof(&bas2.p, &e12.e2, (n + extra) as i32) != 0);
        debug_assert!(test_jac_order_twof(&xy_t2.p2, &e12.e2, (n + extra) as i32) != 0);
    }

    let mut pts: Vec<ThetaPoint> = vec![ThetaPoint::default(); num_p.max(1)];

    let mut space: usize = 1;
    let mut i: u32 = 1;
    while i < n {
        space += 1;
        i *= 2;
    }

    let mut todo: Vec<u16> = vec![0; space];
    todo[0] = (n - 2 + extra) as u16;

    let mut current: i32 = 0;

    let mut jac_q1: Vec<ThetaCoupleJacPoint> = vec![ThetaCoupleJacPoint::default(); space];
    let mut jac_q2: Vec<ThetaCoupleJacPoint> = vec![ThetaCoupleJacPoint::default(); space];
    jac_q1[0] = xy_t1;
    jac_q2[0] = xy_t2;
    while todo[current as usize] != 1 {
        debug_assert!(todo[current as usize] >= 2);
        current += 1;
        debug_assert!((current as usize) < space);
        let prev = todo[(current - 1) as usize];
        let num_dbls: u32 = if prev >= 16 {
            (prev / 2) as u32
        } else {
            (prev - 1) as u32
        };
        debug_assert!(num_dbls > 0 && num_dbls < prev as u32);
        let src1 = jac_q1[(current - 1) as usize];
        double_couple_jac_point_iter(&mut jac_q1[current as usize], num_dbls, &src1, e12);
        let src2 = jac_q2[(current - 1) as usize];
        double_couple_jac_point_iter(&mut jac_q2[current as usize], num_dbls, &src2, e12);
        todo[current as usize] = prev - num_dbls as u16;
    }

    let mut theta_q1: Vec<ThetaPoint> = vec![ThetaPoint::default(); space];
    let mut theta_q2: Vec<ThetaPoint> = vec![ThetaPoint::default(); space];

    // Gluing step.
    let mut first_step = ThetaGluing::default();
    {
        debug_assert_eq!(todo[current as usize], 1);
        if !gluing_compute(
            &mut first_step,
            e12,
            &jac_q1[current as usize],
            &jac_q2[current as usize],
            verify,
        ) {
            return 0;
        }

        for j in 0..num_p {
            debug_assert!(ec_is_zero(&p12[j].p1) != 0 || ec_is_zero(&p12[j].p2) != 0);
            if !gluing_eval_point_special_case(&mut pts[j], &p12[j], &first_step) {
                return 0;
            }
        }

        for j in 0..current as usize {
            gluing_eval_basis(
                &mut theta_q1[j],
                &mut theta_q2[j],
                &jac_q1[j],
                &jac_q2[j],
                &first_step,
            );
            todo[j] -= 1;
        }

        current -= 1;
    }

    theta.null_point = first_step.codomain;
    theta.precomputation = false;
    theta_precomputation(&mut theta);

    let mut step = ThetaIsogeny::default();

    // Remaining generic steps.
    let mut i: u32 = 1;
    while current >= 0 && todo[current as usize] != 0 {
        debug_assert!((current as usize) < space);
        while todo[current as usize] != 1 {
            debug_assert!(todo[current as usize] >= 2);
            current += 1;
            debug_assert!((current as usize) < space);
            let prev = todo[(current - 1) as usize];
            let num_dbls = (prev / 2) as i32;
            debug_assert!(num_dbls > 0 && (num_dbls as u16) < prev);
            let src1 = theta_q1[(current - 1) as usize];
            double_iter(&mut theta_q1[current as usize], &mut theta, &src1, num_dbls);
            let src2 = theta_q2[(current - 1) as usize];
            double_iter(&mut theta_q2[current as usize], &mut theta, &src2, num_dbls);
            todo[current as usize] = prev - num_dbls as u16;
        }

        let ret = if i == n - 2 {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                false,
                false,
                verify,
            )
        } else if i == n - 1 {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                true,
                false,
                false,
            )
        } else {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                false,
                true,
                verify,
            )
        };
        if !ret {
            return 0;
        }

        for pt in pts.iter_mut().take(num_p) {
            let s = *pt;
            theta_isogeny_eval(pt, &step, &s);
        }

        theta = step.codomain;

        debug_assert_eq!(todo[current as usize], 1);
        for j in 0..current as usize {
            let s1 = theta_q1[j];
            theta_isogeny_eval(&mut theta_q1[j], &step, &s1);
            let s2 = theta_q2[j];
            theta_isogeny_eval(&mut theta_q2[j], &step, &s2);
            debug_assert!(todo[j] != 0);
            todo[j] -= 1;
        }

        current -= 1;
        i += 1;
    }

    debug_assert_eq!(current, -1);

    if !extra_torsion {
        if n >= 3 {
            let s1 = theta_q1[0];
            theta_isogeny_eval(&mut theta_q1[0], &step, &s1);
            let s2 = theta_q2[0];
            theta_isogeny_eval(&mut theta_q2[0], &step, &s2);
        }

        let (q1, q2) = (theta_q1[0], theta_q2[0]);
        theta_isogeny_compute_4(&mut step, &theta, &q1, &q2, false, false);
        for pt in pts.iter_mut().take(num_p) {
            let s = *pt;
            theta_isogeny_eval(pt, &step, &s);
        }
        theta = step.codomain;
        let s1 = theta_q1[0];
        theta_isogeny_eval(&mut theta_q1[0], &step, &s1);
        let s2 = theta_q2[0];
        theta_isogeny_eval(&mut theta_q2[0], &step, &s2);

        let (q1, q2) = (theta_q1[0], theta_q2[0]);
        theta_isogeny_compute_2(&mut step, &theta, &q1, &q2, true, false);
        for pt in pts.iter_mut().take(num_p) {
            let s = *pt;
            theta_isogeny_eval(pt, &step, &s);
        }
        theta = step.codomain;
    }

    let mut last_step = ThetaSplitting::default();
    let is_split = splitting_compute(
        &mut last_step,
        &theta,
        if extra_torsion { 8 } else { -1 },
        randomize,
    );
    if !is_split {
        return 0;
    }

    if !theta_product_structure_to_elliptic_product(e34, &last_step.b) {
        return 0;
    }

    for j in 0..num_p {
        let s = pts[j];
        apply_isomorphism(&mut pts[j], &last_step.m, &s);
        if !theta_point_to_montgomery_point(&mut p12[j], &pts[j], &last_step.b) {
            return 0;
        }
    }

    1
}

/// (2ⁿ,2ⁿ)-isogeny chain (signing path: no extra checks, no randomization).
pub fn theta_chain_compute_and_eval(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, false, false)
}

/// (2ⁿ,2ⁿ)-isogeny chain with isotropy checks (verification path).
pub fn theta_chain_compute_and_eval_verify(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, true, false)
}

/// (2ⁿ,2ⁿ)-isogeny chain with randomized codomain model (signing path).
#[cfg(feature = "sign")]
pub fn theta_chain_compute_and_eval_randomized(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, false, true)
}
