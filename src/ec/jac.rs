//! Jacobian-coordinate point arithmetic on Montgomery curves.
//! Port of `ec/ref/lvlx/ec_jac.c`.

use super::{AddComponents, EcCurve, EcPoint, JacPoint};
use crate::gf::*;

/// Initialize as the identity (0:1:0).
pub fn jac_init(p: &mut JacPoint) {
    fp2_set_zero(&mut p.x);
    fp2_set_one(&mut p.y);
    fp2_set_zero(&mut p.z);
}

/// Projective equality on Jacobian (X:Y:Z) ↔ (X/Z², Y/Z³).
pub fn jac_is_equal(p: &JacPoint, q: &JacPoint) -> u32 {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();

    fp2_sqr(&mut t0, &q.z);
    fp2_mul(&mut t2, &p.x, &t0);
    fp2_sqr(&mut t1, &p.z);
    fp2_mul(&mut t3, &q.x, &t1);
    fp2_sub_ip(&mut t2, &t3);

    fp2_mul_ip(&mut t0, &q.z);
    fp2_mul_ip(&mut t0, &p.y);
    fp2_mul_ip(&mut t1, &p.z);
    fp2_mul_ip(&mut t1, &q.y);
    fp2_sub_ip(&mut t0, &t1);

    fp2_is_zero(&t0) & fp2_is_zero(&t2)
}

/// Drop Y to convert (X:Y:Z) → (X:Z²), fixing up (0:1:0) → (1:0).
pub fn jac_to_xz(p: &mut EcPoint, xy: &JacPoint) {
    fp2_copy(&mut p.x, &xy.x);
    fp2_copy(&mut p.z, &xy.z);
    fp2_sqr_ip(&mut p.z);

    let mut one = Fp2::default();
    fp2_set_one(&mut one);
    let c1 = fp2_is_zero(&p.x);
    let c2 = fp2_is_zero(&p.z);
    let px = p.x;
    fp2_select(&mut p.x, &px, &one, c1 & c2);
}

/// Convert Montgomery-Jacobian → short-Weierstrass modified-Jacobian (X:Y:Z:T=a·Z⁴).
pub fn jac_to_ws(q: &mut JacPoint, t: &mut Fp2, ao3: &mut Fp2, p: &JacPoint, curve: &EcCurve) {
    let mut a = Fp2::default();
    let mut one = Fp::default();
    fp_set_one(&mut one);
    if fp2_is_zero(&curve.a) == 0 {
        fp_div3(&mut ao3.re, &curve.a.re);
        fp_div3(&mut ao3.im, &curve.a.im);
        fp2_sqr(t, &p.z);
        let st = *t;
        fp2_mul(&mut q.x, ao3, &st);
        let qx = q.x;
        fp2_add(&mut q.x, &qx, &p.x);
        let st = *t;
        fp2_sqr(t, &st);
        fp2_mul(&mut a, ao3, &curve.a);
        let are = a.re;
        fp_sub(&mut a.re, &one, &are);
        let aim = a.im;
        fp_neg(&mut a.im, &aim);
        let st = *t;
        fp2_mul(t, &st, &a);
    } else {
        fp2_copy(&mut q.x, &p.x);
        fp2_sqr(t, &p.z);
        let st = *t;
        fp2_sqr(t, &st);
    }
    fp2_copy(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

/// Inverse of [`jac_to_ws`].
pub fn jac_from_ws(q: &mut JacPoint, p: &JacPoint, ao3: &Fp2, curve: &EcCurve) {
    if fp2_is_zero(&curve.a) == 0 {
        let mut t = Fp2::default();
        fp2_sqr(&mut t, &p.z);
        fp2_mul_ip(&mut t, ao3);
        fp2_sub(&mut q.x, &p.x, &t);
    }
    fp2_copy(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

#[inline]
pub fn jac_neg(q: &mut JacPoint, p: &JacPoint) {
    fp2_copy(&mut q.x, &p.x);
    fp2_neg(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

/// Jacobian doubling on a Montgomery curve. Cost 6M + 6S.
pub fn jac_dbl(q: &mut JacPoint, p: &JacPoint, ac: &EcCurve) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();

    let flag = fp2_is_zero(&p.x) & fp2_is_zero(&p.z);

    fp2_sqr(&mut t0, &p.x);
    fp2_add(&mut t1, &t0, &t0);
    fp2_add_ip(&mut t0, &t1);
    fp2_sqr(&mut t1, &p.z);
    fp2_mul(&mut t2, &p.x, &ac.a);
    fp2_dbl_ip(&mut t2);
    fp2_add_ip(&mut t2, &t1);
    fp2_mul_ip(&mut t2, &t1);
    fp2_add_ip(&mut t2, &t0);
    fp2_mul(&mut q.z, &p.y, &p.z);
    let qz = q.z;
    fp2_add(&mut q.z, &qz, &qz);
    fp2_sqr(&mut t0, &q.z);
    fp2_mul_ip(&mut t0, &ac.a);
    fp2_sqr(&mut t1, &p.y);
    fp2_dbl_ip(&mut t1);
    fp2_add(&mut t3, &p.x, &p.x);
    fp2_mul_ip(&mut t3, &t1);
    fp2_sqr(&mut q.x, &t2);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t0);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t3);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t3);
    let qx = q.x;
    fp2_sub(&mut q.y, &t3, &qx);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &t2);
    fp2_sqr_ip(&mut t1);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &t1);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &t1);

    // C passes `-flag` (uint32) which yields 1, violating the 0/-1 contract;
    // harmless there because the doubling formula already gives O for O. We
    // pass the mask directly so the select is well-defined.
    let qx = q.x;
    fp2_select(&mut q.x, &qx, &p.x, flag);
    let qz = q.z;
    fp2_select(&mut q.z, &qz, &p.z, flag);
}

/// Weierstrass modified-Jacobian doubling. Cost 3M + 5S.
pub fn jac_dblw(q: &mut JacPoint, u: &mut Fp2, p: &JacPoint, t: &Fp2) {
    let flag = fp2_is_zero(&p.x) & fp2_is_zero(&p.z);

    let mut xx = Fp2::default();
    let mut c = Fp2::default();
    let mut cc = Fp2::default();
    let mut r = Fp2::default();
    let mut s = Fp2::default();
    let mut m = Fp2::default();

    fp2_sqr(&mut xx, &p.x);
    fp2_sqr(&mut c, &p.y);
    let sc = c;
    fp2_add(&mut c, &sc, &sc);
    fp2_sqr(&mut cc, &c);
    fp2_add(&mut r, &cc, &cc);
    fp2_add(&mut s, &p.x, &c);
    let ss = s;
    fp2_sqr(&mut s, &ss);
    let ss = s;
    fp2_sub(&mut s, &ss, &xx);
    let ss = s;
    fp2_sub(&mut s, &ss, &cc);
    fp2_add(&mut m, &xx, &xx);
    let sm = m;
    fp2_add(&mut m, &sm, &xx);
    let sm = m;
    fp2_add(&mut m, &sm, t);
    fp2_sqr(&mut q.x, &m);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &s);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &s);
    fp2_mul(&mut q.z, &p.y, &p.z);
    let qz = q.z;
    fp2_add(&mut q.z, &qz, &qz);
    let qx = q.x;
    fp2_sub(&mut q.y, &s, &qx);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &m);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &r);
    fp2_mul(u, t, &r);
    let su = *u;
    fp2_add(u, &su, &su);

    // C passes `-flag` (uint32) which yields 1, violating the 0/-1 contract;
    // harmless there because the doubling formula already gives O for O. We
    // pass the mask directly so the select is well-defined.
    let qx = q.x;
    fp2_select(&mut q.x, &qx, &p.x, flag);
    let qz = q.z;
    fp2_select(&mut q.z, &qz, &p.z, flag);
}

#[inline]
pub fn select_jac_point(q: &mut JacPoint, p1: &JacPoint, p2: &JacPoint, option: u64) {
    let ctl = option as u32;
    fp2_select(&mut q.x, &p1.x, &p2.x, ctl);
    fp2_select(&mut q.y, &p1.y, &p2.y, ctl);
    fp2_select(&mut q.z, &p1.z, &p2.z, ctl);
}

/// Complete Jacobian addition on a Montgomery curve. Cost 17M + 6S + 13a.
pub fn jac_add(r: &mut JacPoint, p: &JacPoint, q: &JacPoint, ac: &EcCurve) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    let mut u1 = Fp2::default();
    let mut u2 = Fp2::default();
    let mut v1 = Fp2::default();
    let mut dx = Fp2::default();
    let mut dy = Fp2::default();

    let ctl1 = fp2_is_zero(&p.z);
    let ctl2 = fp2_is_zero(&q.z);

    fp2_sqr(&mut t0, &p.z);
    fp2_sqr(&mut t1, &q.z);

    fp2_mul(&mut v1, &t1, &q.z);
    fp2_mul(&mut t2, &t0, &p.z);
    fp2_mul_ip(&mut v1, &p.y);
    fp2_mul_ip(&mut t2, &q.y);
    fp2_sub(&mut dy, &t2, &v1);
    fp2_mul(&mut u2, &t0, &q.x);
    fp2_mul(&mut u1, &t1, &p.x);
    fp2_sub(&mut dx, &u2, &u1);

    fp2_add(&mut t1, &p.y, &p.y);
    fp2_add(&mut t2, &ac.a, &ac.a);
    fp2_mul_ip(&mut t2, &p.x);
    fp2_add_ip(&mut t2, &t0);
    fp2_mul_ip(&mut t2, &t0);
    fp2_sqr(&mut t0, &p.x);
    fp2_add_ip(&mut t2, &t0);
    fp2_add_ip(&mut t2, &t0);
    fp2_add_ip(&mut t2, &t0);
    fp2_mul_ip(&mut t2, &q.z);

    let ctl = fp2_is_zero(&dx) & fp2_is_zero(&dy);
    let sdx = dx;
    fp2_select(&mut dx, &sdx, &t1, ctl);
    let sdy = dy;
    fp2_select(&mut dy, &sdy, &t2, ctl);

    fp2_mul(&mut t0, &p.z, &q.z);
    fp2_sqr(&mut t1, &t0);
    fp2_sqr(&mut t2, &dx);
    fp2_sqr(&mut t3, &dy);

    fp2_mul(&mut r.x, &ac.a, &t1);
    let rx = r.x;
    fp2_add(&mut r.x, &rx, &u1);
    let rx = r.x;
    fp2_add(&mut r.x, &rx, &u2);
    let rx = r.x;
    fp2_mul(&mut r.x, &rx, &t2);
    let rx = r.x;
    fp2_sub(&mut r.x, &t3, &rx);

    let rx = r.x;
    fp2_mul(&mut r.y, &u1, &t2);
    let ry = r.y;
    fp2_sub(&mut r.y, &ry, &rx);
    let ry = r.y;
    fp2_mul(&mut r.y, &ry, &dy);
    fp2_mul(&mut t3, &t2, &dx);
    fp2_mul_ip(&mut t3, &v1);
    let ry = r.y;
    fp2_sub(&mut r.y, &ry, &t3);

    fp2_mul(&mut r.z, &dx, &t0);

    let sr = *r;
    select_jac_point(r, &sr, q, ctl1 as u64);
    let sr = *r;
    select_jac_point(r, &sr, p, ctl2 as u64);
}

/// Compute (u,v,w) such that x(P+Q)=(u-v:w) and x(P-Q)=(u+v:w).
pub fn jac_to_xz_add_components(out: &mut AddComponents, p: &JacPoint, q: &JacPoint, ac: &EcCurve) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    let mut t4 = Fp2::default();
    let mut t5 = Fp2::default();
    let mut t6 = Fp2::default();

    fp2_sqr(&mut t0, &p.z);
    fp2_sqr(&mut t1, &q.z);
    fp2_mul(&mut t2, &p.x, &t1);
    fp2_mul(&mut t3, &t0, &q.x);
    fp2_mul(&mut t4, &p.y, &q.z);
    fp2_mul_ip(&mut t4, &t1);
    fp2_mul(&mut t5, &p.z, &q.y);
    fp2_mul_ip(&mut t5, &t0);
    fp2_mul_ip(&mut t0, &t1);
    fp2_mul(&mut t6, &t4, &t5);
    fp2_add(&mut out.v, &t6, &t6);
    fp2_sqr_ip(&mut t4);
    fp2_sqr_ip(&mut t5);
    fp2_add_ip(&mut t4, &t5);
    fp2_add(&mut t5, &t2, &t3);
    fp2_add(&mut t6, &t3, &t3);
    fp2_rsub_ip(&mut t6, &t5);
    fp2_sqr_ip(&mut t6);
    fp2_mul(&mut t1, &ac.a, &t0);
    fp2_add_ip(&mut t1, &t5);
    fp2_mul_ip(&mut t1, &t6);
    fp2_sub(&mut out.u, &t4, &t1);
    fp2_mul(&mut out.w, &t6, &t0);
}
