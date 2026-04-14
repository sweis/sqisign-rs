// SPDX-License-Identifier: Apache-2.0
//! Vélu-style 2- and 4-isogeny construction and evaluation.
//! Port of `ec/ref/lvlx/{xisog.c, xeval.c}`.

use super::{EcKps2, EcKps4, EcPoint};
use crate::gf::*;

/// Degree-2 isogeny with kernel ⟨P⟩, P ≠ (0,0). Writes new A24 into `b`.
pub fn xisog_2(kps: &mut EcKps2, b: &mut EcPoint, p: &EcPoint) {
    fp2_sqr(&mut b.x, &p.x);
    fp2_sqr(&mut b.z, &p.z);
    let (bx, bz) = (b.x, b.z);
    fp2_sub(&mut b.x, &bz, &bx);
    fp2_add(&mut kps.k.x, &p.x, &p.z);
    fp2_sub(&mut kps.k.z, &p.x, &p.z);
}

/// Degree-2 isogeny with kernel ⟨(0,0)⟩.
pub fn xisog_2_singular(kps: &mut EcKps2, b24: &mut EcPoint, mut a24: EcPoint) {
    let mut t0 = Fp2::default();
    let mut four = Fp2::default();
    fp2_set_small(&mut four, 4);
    fp2_add(&mut t0, &a24.x, &a24.x);
    let s = t0;
    fp2_sub(&mut t0, &s, &a24.z);
    let s = t0;
    fp2_add(&mut t0, &s, &s);
    fp2_inv(&mut a24.z);
    let s = t0;
    fp2_mul(&mut t0, &s, &a24.z);
    fp2_copy(&mut kps.k.x, &t0);
    fp2_add(&mut b24.x, &t0, &t0);
    let s = t0;
    fp2_sqr(&mut t0, &s);
    let s = t0;
    fp2_sub(&mut t0, &s, &four);
    fp2_sqrt(&mut t0);
    fp2_neg(&mut kps.k.z, &t0);
    fp2_add(&mut b24.z, &t0, &t0);
    let (bx, bz) = (b24.x, b24.z);
    fp2_add(&mut b24.x, &bx, &bz);
    let bz = b24.z;
    fp2_add(&mut b24.z, &bz, &bz);
}

/// Degree-4 isogeny with kernel ⟨P⟩, [2]P ≠ (0,0). Writes new A24 into `b`.
pub fn xisog_4(kps: &mut EcKps4, b: &mut EcPoint, p: &EcPoint) {
    let mut k0x = Fp2::default();
    let mut k0z = Fp2::default();
    let mut k1x = Fp2::default();
    let mut k1z = Fp2::default();
    let mut k2x = Fp2::default();

    fp2_sqr(&mut k0x, &p.x);
    fp2_sqr(&mut k0z, &p.z);
    fp2_add(&mut k1x, &k0z, &k0x);
    fp2_sub(&mut k1z, &k0z, &k0x);
    fp2_mul(&mut b.x, &k1x, &k1z);
    fp2_sqr(&mut b.z, &k0z);

    fp2_add(&mut k2x, &p.x, &p.z);
    fp2_sub(&mut k1x, &p.x, &p.z);
    let s = k0z;
    fp2_add(&mut k0x, &s, &s);
    let s = k0x;
    fp2_add(&mut k0x, &s, &s);

    kps.k[0].x = k0x;
    kps.k[0].z = k0z;
    kps.k[1].x = k1x;
    kps.k[1].z = k1z;
    kps.k[2].x = k2x;
}

/// Evaluate a degree-2 isogeny on a slice of points.
pub fn xeval_2(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    for j in 0..q.len() {
        fp2_add(&mut t0, &q[j].x, &q[j].z);
        fp2_sub(&mut t1, &q[j].x, &q[j].z);
        fp2_mul(&mut t2, &kps.k.x, &t1);
        fp2_mul(&mut t1, &kps.k.z, &t0);
        fp2_add(&mut t0, &t2, &t1);
        let s = t1;
        fp2_sub(&mut t1, &t2, &s);
        fp2_mul(&mut r[j].x, &q[j].x, &t0);
        fp2_mul(&mut r[j].z, &q[j].z, &t1);
    }
}

/// Evaluate the singular degree-2 isogeny on a slice of points.
pub fn xeval_2_singular(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    for i in 0..q.len() {
        fp2_mul(&mut t0, &q[i].x, &q[i].z);
        fp2_mul(&mut t1, &kps.k.x, &q[i].z);
        let s = t1;
        fp2_add(&mut t1, &s, &q[i].x);
        let s = t1;
        fp2_mul(&mut t1, &s, &q[i].x);
        fp2_sqr(&mut r[i].x, &q[i].z);
        let rx = r[i].x;
        fp2_add(&mut r[i].x, &rx, &t1);
        fp2_mul(&mut r[i].z, &t0, &kps.k.z);
    }
}

/// Evaluate a degree-4 isogeny on a slice of points.
pub fn xeval_4(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps4) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    for i in 0..q.len() {
        fp2_add(&mut t0, &q[i].x, &q[i].z);
        fp2_sub(&mut t1, &q[i].x, &q[i].z);
        fp2_mul(&mut r[i].x, &t0, &kps.k[1].x);
        fp2_mul(&mut r[i].z, &t1, &kps.k[2].x);
        let s = t0;
        fp2_mul(&mut t0, &s, &t1);
        let s = t0;
        fp2_mul(&mut t0, &s, &kps.k[0].x);
        let (rx, rz) = (r[i].x, r[i].z);
        fp2_add(&mut t1, &rx, &rz);
        fp2_sub(&mut r[i].z, &rx, &rz);
        let s = t1;
        fp2_sqr(&mut t1, &s);
        let rz = r[i].z;
        fp2_sqr(&mut r[i].z, &rz);
        fp2_add(&mut r[i].x, &t0, &t1);
        let rz = r[i].z;
        let s = t0;
        fp2_sub(&mut t0, &s, &rz);
        let rx = r[i].x;
        fp2_mul(&mut r[i].x, &rx, &t1);
        let rz = r[i].z;
        fp2_mul(&mut r[i].z, &rz, &t0);
    }
}

/// In-place wrapper for `xeval_2` (handles C aliasing pattern `xeval_2(R, R, ...)`).
#[inline]
pub fn xeval_2_inplace(r: &mut [EcPoint], kps: &EcKps2) {
    let q: Vec<EcPoint> = r.to_vec();
    xeval_2(r, &q, kps);
}

/// In-place wrapper for `xeval_2_singular`.
#[inline]
pub fn xeval_2_singular_inplace(r: &mut [EcPoint], kps: &EcKps2) {
    let q: Vec<EcPoint> = r.to_vec();
    xeval_2_singular(r, &q, kps);
}

/// In-place wrapper for `xeval_4`.
#[inline]
pub fn xeval_4_inplace(r: &mut [EcPoint], kps: &EcKps4) {
    let q: Vec<EcPoint> = r.to_vec();
    xeval_4(r, &q, kps);
}
