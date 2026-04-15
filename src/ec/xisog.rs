//! Vélu-style 2- and 4-isogeny construction and evaluation.
//! Port of `ec/ref/lvlx/{xisog.c, xeval.c}`.

use super::{EcKps2, EcKps4, EcPoint};
use crate::gf::*;

/// Degree-2 isogeny with kernel ⟨P⟩, P ≠ (0,0). Writes new A24 into `b`.
pub fn xisog_2(kps: &mut EcKps2, b: &mut EcPoint, p: &EcPoint) {
    b.x = p.x.square();
    b.z = p.z.square();
    let (bx, bz) = (b.x, b.z);
    b.x = bz - bx;
    kps.k.x = p.x + p.z;
    kps.k.z = p.x - p.z;
}

/// Degree-2 isogeny with kernel ⟨(0,0)⟩.
pub fn xisog_2_singular(kps: &mut EcKps2, b24: &mut EcPoint, mut a24: EcPoint) {
    let four = Fp2::from_small(4);
    let mut t0 = a24.x + a24.x;
    t0 -= a24.z;
    t0.dbl_ip();
    a24.z = (a24.z).inv();
    t0 *= a24.z;
    kps.k.x = t0;
    b24.x = t0 + t0;
    t0.square_ip();
    t0 -= four;
    t0 = t0.sqrt();
    kps.k.z = -t0;
    b24.z = t0 + t0;
    let (bx, bz) = (b24.x, b24.z);
    b24.x = bx + bz;
    b24.z += bz;
}

/// Degree-4 isogeny with kernel ⟨P⟩, [2]P ≠ (0,0). Writes new A24 into `b`.
pub fn xisog_4(kps: &mut EcKps4, b: &mut EcPoint, p: &EcPoint) {
    let mut k0x = p.x.square();
    let k0z = p.z.square();
    let mut k1x = k0z + k0x;
    let k1z = k0z - k0x;
    b.x = k1x * k1z;
    b.z = k0z.square();

    let k2x = p.x + p.z;
    k1x = p.x - p.z;
    let s = k0z;
    k0x = s + s;
    k0x.dbl_ip();

    kps.k[0].x = k0x;
    kps.k[0].z = k0z;
    kps.k[1].x = k1x;
    kps.k[1].z = k1z;
    kps.k[2].x = k2x;
}

/// Evaluate a degree-2 isogeny on a slice of points.
pub fn xeval_2(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    for (r, q) in r.iter_mut().zip(q) {
        let t2 = kps.k.x * (q.x - q.z);
        let t1 = kps.k.z * (q.x + q.z);
        r.x = q.x * (t2 + t1);
        r.z = q.z * (t2 - t1);
    }
}

/// Evaluate the singular degree-2 isogeny on a slice of points.
pub fn xeval_2_singular(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    for (r, q) in r.iter_mut().zip(q) {
        let t1 = (kps.k.x * q.z + q.x) * q.x;
        r.x = q.z.square() + t1;
        r.z = q.x * q.z * kps.k.z;
    }
}

/// Evaluate a degree-4 isogeny on a slice of points.
pub fn xeval_4(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps4) {
    for i in 0..q.len() {
        let mut t0 = q[i].x + q[i].z;
        let t1 = q[i].x - q[i].z;
        r[i].x = t0 * kps.k[1].x;
        r[i].z = t1 * kps.k[2].x;
        t0 *= t1;
        t0 *= kps.k[0].x;
        let (rx, rz) = (r[i].x, r[i].z);
        let t1 = (rx + rz).square();
        r[i].z = (rx - rz).square();
        r[i].x = (t0 + t1) * t1;
        r[i].z *= t0 - r[i].z;
    }
}

/// In-place wrapper for `xeval_2` (handles C aliasing pattern `xeval_2(R, R, ...)`).
/// Each output element depends only on the same-index input, so a per-element
/// snapshot suffices — no heap allocation.
#[inline]
pub fn xeval_2_inplace(r: &mut [EcPoint], kps: &EcKps2) {
    for j in 0..r.len() {
        let q = r[j];
        xeval_2(
            core::slice::from_mut(&mut r[j]),
            core::slice::from_ref(&q),
            kps,
        );
    }
}

/// In-place wrapper for `xeval_2_singular`.
#[inline]
pub fn xeval_2_singular_inplace(r: &mut [EcPoint], kps: &EcKps2) {
    for j in 0..r.len() {
        let q = r[j];
        xeval_2_singular(
            core::slice::from_mut(&mut r[j]),
            core::slice::from_ref(&q),
            kps,
        );
    }
}

/// In-place wrapper for `xeval_4`.
#[inline]
pub fn xeval_4_inplace(r: &mut [EcPoint], kps: &EcKps4) {
    for j in 0..r.len() {
        let q = r[j];
        xeval_4(
            core::slice::from_mut(&mut r[j]),
            core::slice::from_ref(&q),
            kps,
        );
    }
}
