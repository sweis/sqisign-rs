// SPDX-License-Identifier: Apache-2.0
//! (2,2)-isogenies between products of elliptic curves via the theta model.
//!
//! Ported from `src/hd/ref/` in the C reference implementation.

use crate::ec::{
    copy_point, ec_dbl, ec_is_equal, ec_is_zero, jac_add, jac_dbl, jac_dblw, jac_from_ws,
    jac_to_ws, jac_to_xz, test_point_order_twof, AddComponents, EcBasis, EcCurve, EcPoint,
    JacPoint,
};
use crate::gf::{
    fp2_add, fp2_add_ip, fp2_copy, fp2_dbl_ip, fp2_is_equal, fp2_is_zero, fp2_mul, fp2_mul_ip,
    fp2_neg, fp2_neg_ip, fp2_rsub_ip, fp2_select, fp2_set_one, fp2_set_zero, fp2_sqr, fp2_sqr_ip,
    fp2_sub, fp2_sub_ip, Fp2,
};

mod theta_isogenies;
pub use theta_isogenies::*;

/// Extra two-power torsion above the kernel needed by the chain algorithm.
pub const HD_EXTRA_TORSION: u32 = 2;

// ===========================================================================
// Data structures (hd.h)
// ===========================================================================

/// A pair of x-only points (P₁, P₂) ∈ E₁ × E₂.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaCouplePoint {
    pub p1: EcPoint,
    pub p2: EcPoint,
}

/// Kernel data for a (2,2)-isogeny: T₁, T₂ and T₁−T₂ on E₁ × E₂.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaKernelCouplePoints {
    pub t1: ThetaCouplePoint,
    pub t2: ThetaCouplePoint,
    pub t1m2: ThetaCouplePoint,
}

/// A pair of Jacobian points on E₁ × E₂.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaCoupleJacPoint {
    pub p1: JacPoint,
    pub p2: JacPoint,
}

/// An ordered product of two Montgomery curves.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaCoupleCurve {
    pub e1: EcCurve,
    pub e2: EcCurve,
}

/// A product E₁ × E₂ together with a 2ⁿ-torsion basis on each factor.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaCoupleCurveWithBasis {
    pub e1: EcCurve,
    pub e2: EcCurve,
    pub b1: EcBasis,
    pub b2: EcBasis,
}

/// A level-2 theta point with four projective coordinates.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaPoint {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
    pub t: Fp2,
}

/// A theta point that has only two distinct components.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaPointCompact {
    pub x: Fp2,
    pub y: Fp2,
}

/// A theta null point with cached doubling/isogeny precomputation.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaStructure {
    pub null_point: ThetaPoint,
    pub precomputation: bool,
    pub xyz0_d: Fp2,
    pub yzt0_d: Fp2,
    pub xzt0_d: Fp2,
    pub xyt0_d: Fp2,
    pub xyz0: Fp2,
    pub yzt0: Fp2,
    pub xzt0: Fp2,
    pub xyt0: Fp2,
}

/// 2×2 matrix used for action by translation during gluing.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct TranslationMatrix {
    pub g00: Fp2,
    pub g01: Fp2,
    pub g10: Fp2,
    pub g11: Fp2,
}

/// 4×4 change-of-basis matrix for theta points.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct BasisChangeMatrix {
    pub m: [[Fp2; 4]; 4],
}

/// State for the gluing isogeny E₁ × E₂ → A.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaGluing {
    pub domain: ThetaCoupleCurve,
    pub xy_k1_8: ThetaCoupleJacPoint,
    pub image_k1_8: ThetaPointCompact,
    pub m: BasisChangeMatrix,
    pub precomputation: ThetaPoint,
    pub codomain: ThetaPoint,
}

/// State for a generic (2,2)-isogeny step in the theta model.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaIsogeny {
    pub t1_8: ThetaPoint,
    pub t2_8: ThetaPoint,
    pub hadamard_bool_1: bool,
    pub hadamard_bool_2: bool,
    pub domain: ThetaStructure,
    pub precomputation: ThetaPoint,
    pub codomain: ThetaStructure,
}

/// State for the splitting isomorphism back to an elliptic product.
#[derive(Clone, Copy, Default, Debug)]
#[repr(C)]
pub struct ThetaSplitting {
    pub m: BasisChangeMatrix,
    pub b: ThetaStructure,
}

// ===========================================================================
// Couple-point arithmetic on E₁ × E₂ (hd.c)
// ===========================================================================

/// out ← [2]·in on E₁ × E₂.
pub fn double_couple_point(
    out: &mut ThetaCouplePoint,
    in_: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) {
    ec_dbl(&mut out.p1, &in_.p1, &e12.e1);
    ec_dbl(&mut out.p2, &in_.p2, &e12.e2);
}

/// out ← [2ⁿ]·in on E₁ × E₂.
pub fn double_couple_point_iter(
    out: &mut ThetaCouplePoint,
    n: u32,
    in_: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) {
    if n == 0 {
        *out = *in_;
    } else {
        double_couple_point(out, in_, e12);
        for _ in 0..n - 1 {
            let s = *out;
            double_couple_point(out, &s, e12);
        }
    }
}

/// out ← T₁ + T₂ on E₁ × E₂ (Jacobian).
pub fn add_couple_jac_points(
    out: &mut ThetaCoupleJacPoint,
    t1: &ThetaCoupleJacPoint,
    t2: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    jac_add(&mut out.p1, &t1.p1, &t2.p1, &e12.e1);
    jac_add(&mut out.p2, &t1.p2, &t2.p2, &e12.e2);
}

/// out ← [2]·in on E₁ × E₂ (Jacobian).
pub fn double_couple_jac_point(
    out: &mut ThetaCoupleJacPoint,
    in_: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    jac_dbl(&mut out.p1, &in_.p1, &e12.e1);
    jac_dbl(&mut out.p2, &in_.p2, &e12.e2);
}

/// out ← [2ⁿ]·in on E₁ × E₂ (Jacobian).
///
/// For n ≥ 2, converts to short-Weierstrass form, uses the cheaper `DBLW`
/// repeatedly, then converts back.
pub fn double_couple_jac_point_iter(
    out: &mut ThetaCoupleJacPoint,
    n: u32,
    in_: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    if n == 0 {
        *out = *in_;
    } else if n == 1 {
        double_couple_jac_point(out, in_, e12);
    } else {
        let mut a1 = Fp2::default();
        let mut a2 = Fp2::default();
        let mut t1 = Fp2::default();
        let mut t2 = Fp2::default();

        jac_to_ws(&mut out.p1, &mut t1, &mut a1, &in_.p1, &e12.e1);
        jac_to_ws(&mut out.p2, &mut t2, &mut a2, &in_.p2, &e12.e2);

        let (sp, st) = (out.p1, t1);
        jac_dblw(&mut out.p1, &mut t1, &sp, &st);
        let (sp, st) = (out.p2, t2);
        jac_dblw(&mut out.p2, &mut t2, &sp, &st);
        for _ in 0..n - 1 {
            let (sp, st) = (out.p1, t1);
            jac_dblw(&mut out.p1, &mut t1, &sp, &st);
            let (sp, st) = (out.p2, t2);
            jac_dblw(&mut out.p2, &mut t2, &sp, &st);
        }

        let s = out.p1;
        jac_from_ws(&mut out.p1, &s, &a1, &e12.e1);
        let s = out.p2;
        jac_from_ws(&mut out.p2, &s, &a2, &e12.e2);
    }
}

/// Forget the Y-coordinates of a Jacobian couple point.
pub fn couple_jac_to_xz(p: &mut ThetaCouplePoint, xyp: &ThetaCoupleJacPoint) {
    jac_to_xz(&mut p.p1, &xyp.p1);
    jac_to_xz(&mut p.p2, &xyp.p2);
}

/// Pack two bases (P,Q,P−Q) on E₁, E₂ into a kernel triple on E₁ × E₂.
pub fn copy_bases_to_kernel(ker: &mut ThetaKernelCouplePoints, b1: &EcBasis, b2: &EcBasis) {
    copy_point(&mut ker.t1.p1, &b1.p);
    copy_point(&mut ker.t2.p1, &b1.q);
    copy_point(&mut ker.t1m2.p1, &b1.pmq);

    copy_point(&mut ker.t1.p2, &b2.p);
    copy_point(&mut ker.t2.p2, &b2.q);
    copy_point(&mut ker.t1m2.p2, &b2.pmq);
}

/// Debug helper: tests both components of a couple point have order exactly 2ᵗ.
pub fn test_couple_point_order_twof(
    t: &ThetaCouplePoint,
    e: &ThetaCoupleCurve,
    tw: i32,
) -> u32 {
    test_point_order_twof(&t.p1, &e.e1, tw) & test_point_order_twof(&t.p2, &e.e2, tw)
}

// ===========================================================================
// Theta structure primitives (theta_structure.{h,c})
// ===========================================================================

/// Hadamard transform: (x,y,z,t) ↦ (x+y+z+t, x−y+z−t, x+y−z−t, x−y−z+t).
#[inline]
pub fn hadamard(out: &mut ThetaPoint, in_: &ThetaPoint) {
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    let mut t4 = Fp2::default();

    fp2_add(&mut t1, &in_.x, &in_.y);
    fp2_sub(&mut t2, &in_.x, &in_.y);
    fp2_add(&mut t3, &in_.z, &in_.t);
    fp2_sub(&mut t4, &in_.z, &in_.t);

    fp2_add(&mut out.x, &t1, &t3);
    fp2_add(&mut out.y, &t2, &t4);
    fp2_sub(&mut out.z, &t1, &t3);
    fp2_sub(&mut out.t, &t2, &t4);
}

/// Coordinate-wise squaring.
#[inline]
pub fn pointwise_square(out: &mut ThetaPoint, in_: &ThetaPoint) {
    fp2_sqr(&mut out.x, &in_.x);
    fp2_sqr(&mut out.y, &in_.y);
    fp2_sqr(&mut out.z, &in_.z);
    fp2_sqr(&mut out.t, &in_.t);
}

/// Square coordinates then apply the Hadamard transform.
#[inline]
pub fn to_squared_theta(out: &mut ThetaPoint, in_: &ThetaPoint) {
    pointwise_square(out, in_);
    hadamard_ip(out);
}

/// In-place Hadamard. The transform reads all inputs into locals before
/// writing, so it is alias-safe.
#[inline]
pub fn hadamard_ip(r: &mut ThetaPoint) {
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    let mut t4 = Fp2::default();
    fp2_add(&mut t1, &r.x, &r.y);
    fp2_sub(&mut t2, &r.x, &r.y);
    fp2_add(&mut t3, &r.z, &r.t);
    fp2_sub(&mut t4, &r.z, &r.t);
    fp2_add(&mut r.x, &t1, &t3);
    fp2_add(&mut r.y, &t2, &t4);
    fp2_sub(&mut r.z, &t1, &t3);
    fp2_sub(&mut r.t, &t2, &t4);
}

/// In-place coordinate-wise squaring.
#[inline]
pub fn pointwise_square_ip(r: &mut ThetaPoint) {
    fp2_sqr_ip(&mut r.x);
    fp2_sqr_ip(&mut r.y);
    fp2_sqr_ip(&mut r.z);
    fp2_sqr_ip(&mut r.t);
}

/// In-place `to_squared_theta`.
#[inline]
pub fn to_squared_theta_ip(r: &mut ThetaPoint) {
    pointwise_square_ip(r);
    hadamard_ip(r);
}

/// Fill the eight cached products on a theta structure.
pub fn theta_precomputation(a: &mut ThetaStructure) {
    if a.precomputation {
        return;
    }

    let mut a_dual = ThetaPoint::default();
    to_squared_theta(&mut a_dual, &a.null_point);

    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_mul(&mut t1, &a_dual.x, &a_dual.y);
    fp2_mul(&mut t2, &a_dual.z, &a_dual.t);
    fp2_mul(&mut a.xyz0_d, &t1, &a_dual.z);
    fp2_mul(&mut a.xyt0_d, &t1, &a_dual.t);
    fp2_mul(&mut a.yzt0_d, &t2, &a_dual.y);
    fp2_mul(&mut a.xzt0_d, &t2, &a_dual.x);

    fp2_mul(&mut t1, &a.null_point.x, &a.null_point.y);
    fp2_mul(&mut t2, &a.null_point.z, &a.null_point.t);
    fp2_mul(&mut a.xyz0, &t1, &a.null_point.z);
    fp2_mul(&mut a.xyt0, &t1, &a.null_point.t);
    fp2_mul(&mut a.yzt0, &t2, &a.null_point.y);
    fp2_mul(&mut a.xzt0, &t2, &a.null_point.x);

    a.precomputation = true;
}

/// out ← [2]·in on the theta structure A.
pub fn double_point(out: &mut ThetaPoint, a: &mut ThetaStructure, in_: &ThetaPoint) {
    to_squared_theta(out, in_);
    fp2_sqr_ip(&mut out.x);
    fp2_sqr_ip(&mut out.y);
    fp2_sqr_ip(&mut out.z);
    fp2_sqr_ip(&mut out.t);

    if !a.precomputation {
        theta_precomputation(a);
    }
    fp2_mul_ip(&mut out.x, &a.yzt0_d);
    fp2_mul_ip(&mut out.y, &a.xzt0_d);
    fp2_mul_ip(&mut out.z, &a.xyt0_d);
    fp2_mul_ip(&mut out.t, &a.xyz0_d);

    hadamard_ip(out);

    fp2_mul_ip(&mut out.x, &a.yzt0);
    fp2_mul_ip(&mut out.y, &a.xzt0);
    fp2_mul_ip(&mut out.z, &a.xyt0);
    fp2_mul_ip(&mut out.t, &a.xyz0);
}

/// out ← [2ᵉˣᵖ]·in on the theta structure A.
pub fn double_iter(out: &mut ThetaPoint, a: &mut ThetaStructure, in_: &ThetaPoint, exp: i32) {
    if exp == 0 {
        *out = *in_;
    } else {
        double_point(out, a, in_);
        for _ in 1..exp {
            let s = *out;
            double_point(out, a, &s);
        }
    }
}

/// Returns 0xFFFFFFFF if x·t = y·z (i.e. P lies on a product theta structure).
pub fn is_product_theta_point(p: &ThetaPoint) -> u32 {
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    fp2_mul(&mut t1, &p.x, &p.t);
    fp2_mul(&mut t2, &p.y, &p.z);
    fp2_is_equal(&t1, &t2)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::{fp2_set_one, fp2_set_small};

    fn tp(x: u64, y: u64, z: u64, t: u64) -> ThetaPoint {
        let mut p = ThetaPoint::default();
        fp2_set_small(&mut p.x, x);
        fp2_set_small(&mut p.y, y);
        fp2_set_small(&mut p.z, z);
        fp2_set_small(&mut p.t, t);
        p
    }

    #[test]
    fn hadamard_involution_up_to_scale() {
        // H² = 4·id for the 4×4 Hadamard matrix used here.
        let p = tp(3, 7, 11, 19);
        let mut hp = ThetaPoint::default();
        hadamard(&mut hp, &p);
        let mut hhp = ThetaPoint::default();
        hadamard(&mut hhp, &hp);
        let expect = tp(12, 28, 44, 76);
        assert_eq!(fp2_is_equal(&hhp.x, &expect.x), 0xFFFFFFFF);
        assert_eq!(fp2_is_equal(&hhp.y, &expect.y), 0xFFFFFFFF);
        assert_eq!(fp2_is_equal(&hhp.z, &expect.z), 0xFFFFFFFF);
        assert_eq!(fp2_is_equal(&hhp.t, &expect.t), 0xFFFFFFFF);
    }

    #[test]
    fn product_theta_point_detection() {
        // (a, b, c, bc/a) satisfies x·t = y·z. Use a=1 to avoid division.
        let p = tp(1, 5, 7, 35);
        assert_eq!(is_product_theta_point(&p), 0xFFFFFFFF);
        let q = tp(1, 5, 7, 36);
        assert_eq!(is_product_theta_point(&q), 0);
        let z = ThetaPoint::default();
        assert_eq!(is_product_theta_point(&z), 0xFFFFFFFF);
    }

    #[test]
    fn to_squared_theta_explicit() {
        // (1,1,1,1) → squares (1,1,1,1) → Hadamard (4,0,0,0).
        let p = tp(1, 1, 1, 1);
        let mut out = ThetaPoint::default();
        to_squared_theta(&mut out, &p);
        let mut four = Fp2::default();
        fp2_set_small(&mut four, 4);
        assert_eq!(fp2_is_equal(&out.x, &four), 0xFFFFFFFF);
        assert_eq!(fp2_is_zero(&out.y), 0xFFFFFFFF);
        assert_eq!(fp2_is_zero(&out.z), 0xFFFFFFFF);
        assert_eq!(fp2_is_zero(&out.t), 0xFFFFFFFF);
    }

    #[test]
    fn copy_bases_to_kernel_layout() {
        let mut one = Fp2::default();
        fp2_set_one(&mut one);
        let p = EcPoint { x: one, z: one };
        let mut b1 = EcBasis::default();
        b1.p = p; b1.q = p; b1.pmq = p;
        let b2 = EcBasis::default();
        let mut ker = ThetaKernelCouplePoints::default();
        copy_bases_to_kernel(&mut ker, &b1, &b2);
        assert_eq!(fp2_is_equal(&ker.t1.p1.x, &one), 0xFFFFFFFF);
        assert_eq!(fp2_is_zero(&ker.t1.p2.x), 0xFFFFFFFF);
    }
}
