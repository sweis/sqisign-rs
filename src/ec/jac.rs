//! Jacobian-coordinate point arithmetic on Montgomery curves.
//! Port of `ec/ref/lvlx/ec_jac.c`.

use super::{AddComponents, EcCurve, EcPoint, JacPoint};
use crate::gf::*;
use core::ops::Neg;

impl JacPoint {
    /// The Jacobian identity (0:1:0).
    pub const IDENTITY: Self = Self {
        x: Fp2::ZERO,
        y: Fp2::ONE,
        z: Fp2::ZERO,
    };
}

impl Neg for JacPoint {
    type Output = Self;
    #[inline]
    fn neg(mut self) -> Self {
        self.y.neg_ip();
        self
    }
}

impl From<JacPoint> for EcPoint {
    /// Drop Y to convert (X:Y:Z) → (X:Z²), fixing up (0:1:0) → (1:0).
    fn from(j: JacPoint) -> Self {
        let mut p = EcPoint::default();
        jac_to_xz(&mut p, &j);
        p
    }
}

/// Projective equality on Jacobian (X:Y:Z) ↔ (X/Z², Y/Z³).
pub fn jac_is_equal(p: &JacPoint, q: &JacPoint) -> bool {
    let mut t0 = q.z.square();
    let mut t2 = p.x * t0;
    let mut t1 = p.z.square();
    let t3 = q.x * t1;
    t2 -= t3;

    t0 *= q.z;
    t0 *= p.y;
    t1 *= p.z;
    t1 *= q.y;
    t0 -= t1;

    (t0.is_zero_ct() & t2.is_zero_ct()) != 0
}

/// Drop Y to convert (X:Y:Z) → (X:Z²), fixing up (0:1:0) → (1:0).
pub fn jac_to_xz(p: &mut EcPoint, xy: &JacPoint) {
    p.x = xy.x;
    p.z = xy.z;
    p.z.square_ip();

    let c1 = p.x.is_zero_ct();
    let c2 = p.z.is_zero_ct();
    p.x = Fp2::select(&p.x, &Fp2::ONE, c1 & c2);
}

/// Convert Montgomery-Jacobian → short-Weierstrass modified-Jacobian (X:Y:Z:T=a·Z⁴).
pub fn jac_to_ws(q: &mut JacPoint, t: &mut Fp2, ao3: &mut Fp2, p: &JacPoint, curve: &EcCurve) {
    if curve.a.is_zero() {
        q.x = p.x;
        *t = p.z.square().square();
    } else {
        ao3.re = curve.a.re.div3();
        ao3.im = curve.a.im.div3();
        *t = p.z.square();
        q.x = *ao3 * *t + p.x;
        t.square_ip();
        let a = Fp2::ONE - *ao3 * curve.a;
        *t *= a;
    }
    q.y = p.y;
    q.z = p.z;
}

/// Inverse of [`jac_to_ws`].
pub fn jac_from_ws(q: &mut JacPoint, p: &JacPoint, ao3: &Fp2, curve: &EcCurve) {
    if curve.a.is_zero_ct() == 0 {
        let mut t = p.z.square();
        t *= ao3;
        q.x = p.x - t;
    }
    q.y = p.y;
    q.z = p.z;
}

/// Jacobian doubling on a Montgomery curve. Cost 6M + 6S.
pub fn jac_dbl(q: &mut JacPoint, p: &JacPoint, ac: &EcCurve) {
    let flag = p.x.is_zero_ct() & p.z.is_zero_ct();

    let mut t0 = p.x.square();
    let mut t1 = t0 + t0;
    t0 += t1;
    t1 = p.z.square();
    let mut t2 = p.x * ac.a;
    t2.dbl_ip();
    t2 += t1;
    t2 *= t1;
    t2 += t0;
    q.z = p.y * p.z;
    q.z.dbl_ip();
    t0 = q.z.square();
    t0 *= ac.a;
    t1 = p.y.square();
    t1.dbl_ip();
    let mut t3 = p.x + p.x;
    t3 *= t1;
    q.x = t2.square();
    q.x -= t0;
    q.x -= t3;
    q.x -= t3;
    let qx = q.x;
    q.y = t3 - qx;
    q.y *= t2;
    t1.square_ip();
    q.y -= t1;
    q.y -= t1;

    // C passes `-flag` (uint32) which yields 1, violating the 0/-1 contract;
    // harmless there because the doubling formula already gives O for O. We
    // pass the mask directly so the select is well-defined.
    q.x = Fp2::select(&q.x, &p.x, flag);
    q.z = Fp2::select(&q.z, &p.z, flag);
}

/// Weierstrass modified-Jacobian doubling. Cost 3M + 5S.
pub fn jac_dblw(q: &mut JacPoint, u: &mut Fp2, p: &JacPoint, t: &Fp2) {
    let flag = p.x.is_zero_ct() & p.z.is_zero_ct();
    let _cc = Fp2::default();
    let xx = p.x.square();
    let c = p.y.square().dbl();
    let cc = c.square();
    let r = cc.dbl();
    let s = (p.x + c).square() - xx - cc;
    let m = xx.dbl() + xx + t;
    q.x = m.square() - s - s;
    q.z = (p.y * p.z).dbl();
    q.y = (s - q.x) * m - r;
    *u = (t * r).dbl();

    // C passes `-flag` (uint32) which yields 1, violating the 0/-1 contract;
    // harmless there because the doubling formula already gives O for O. We
    // pass the mask directly so the select is well-defined.
    q.x = Fp2::select(&q.x, &p.x, flag);
    q.z = Fp2::select(&q.z, &p.z, flag);
}

#[inline]
pub fn select_jac_point(q: &mut JacPoint, p1: &JacPoint, p2: &JacPoint, option: u64) {
    let ctl = option as u32;
    q.x = Fp2::select(&p1.x, &p2.x, ctl);
    q.y = Fp2::select(&p1.y, &p2.y, ctl);
    q.z = Fp2::select(&p1.z, &p2.z, ctl);
}

/// Complete Jacobian addition on a Montgomery curve. Cost 17M + 6S + 13a.
pub fn jac_add(r: &mut JacPoint, p: &JacPoint, q: &JacPoint, ac: &EcCurve) {
    let ctl1 = p.z.is_zero_ct();
    let ctl2 = q.z.is_zero_ct();

    let mut t0 = p.z.square();
    let mut t1 = q.z.square();

    let mut v1 = t1 * q.z;
    let mut t2 = t0 * p.z;
    v1 *= p.y;
    t2 *= q.y;
    let mut dy = t2 - v1;
    let u2 = t0 * q.x;
    let u1 = t1 * p.x;
    let mut dx = u2 - u1;

    t1 = p.y + p.y;
    t2 = ac.a + ac.a;
    t2 *= p.x;
    t2 += t0;
    t2 *= t0;
    t0 = p.x.square();
    t2 += t0;
    t2 += t0;
    t2 += t0;
    t2 *= q.z;

    let ctl = dx.is_zero_ct() & dy.is_zero_ct();
    dx = Fp2::select(&dx, &t1, ctl);
    dy = Fp2::select(&dy, &t2, ctl);

    t0 = p.z * q.z;
    t1 = t0.square();
    t2 = dx.square();
    let mut t3 = dy.square();

    r.x = ac.a * t1;
    r.x += u1;
    r.x += u2;
    r.x *= t2;
    r.x = t3 - r.x;

    let rx = r.x;
    r.y = u1 * t2;
    r.y -= rx;
    r.y *= dy;
    t3 = t2 * dx;
    t3 *= v1;
    r.y -= t3;

    r.z = dx * t0;

    let sr = *r;
    select_jac_point(r, &sr, q, ctl1 as u64);
    let sr = *r;
    select_jac_point(r, &sr, p, ctl2 as u64);
}

/// Compute (u,v,w) such that x(P+Q)=(u-v:w) and x(P-Q)=(u+v:w).
pub fn jac_to_xz_add_components(out: &mut AddComponents, p: &JacPoint, q: &JacPoint, ac: &EcCurve) {
    let mut t0 = p.z.square();
    let mut t1 = q.z.square();
    let t2 = p.x * t1;
    let t3 = t0 * q.x;
    let mut t4 = p.y * q.z;
    t4 *= t1;
    let mut t5 = p.z * q.y;
    t5 *= t0;
    t0 *= t1;
    let mut t6 = t4 * t5;
    out.v = t6 + t6;
    t4.square_ip();
    t5.square_ip();
    t4 += t5;
    t5 = t2 + t3;
    t6 = t3 + t3;
    t6 = t5 - t6;
    t6.square_ip();
    t1 = ac.a * t0;
    t1 += t5;
    t1 *= t6;
    out.u = t4 - t1;
    out.w = t6 * t0;
}
