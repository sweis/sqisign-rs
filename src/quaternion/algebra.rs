// SPDX-License-Identifier: Apache-2.0
//! Quaternion algebra element arithmetic, ported from `algebra.c`.

use super::dim4::*;
use super::intbig::*;
use super::types::*;

/// Coordinate-only quaternion product in basis (1, i, j, ij) with
/// i² = -1, j² = -p, k = ij.
pub fn quat_alg_coord_mul(res: &mut IbzVec4, a: &IbzVec4, b: &IbzVec4, alg: &QuatAlg) {
    let mut prod = Ibz::new();
    let mut sum = ibz_vec_4_init();

    // 1 component: -p(a2 b2 + a3 b3) + a0 b0 - a1 b1
    ibz_mul(&mut prod, &a[2], &b[2]);
    let t = sum[0].clone();
    ibz_sub(&mut sum[0], &t, &prod);
    ibz_mul(&mut prod, &a[3], &b[3]);
    let t = sum[0].clone();
    ibz_sub(&mut sum[0], &t, &prod);
    let t = sum[0].clone();
    ibz_mul(&mut sum[0], &t, &alg.p);
    ibz_mul(&mut prod, &a[0], &b[0]);
    let t = sum[0].clone();
    ibz_add(&mut sum[0], &t, &prod);
    ibz_mul(&mut prod, &a[1], &b[1]);
    let t = sum[0].clone();
    ibz_sub(&mut sum[0], &t, &prod);

    // i component: p(a2 b3 - a3 b2) + a0 b1 + a1 b0
    ibz_mul(&mut prod, &a[2], &b[3]);
    let t = sum[1].clone();
    ibz_add(&mut sum[1], &t, &prod);
    ibz_mul(&mut prod, &a[3], &b[2]);
    let t = sum[1].clone();
    ibz_sub(&mut sum[1], &t, &prod);
    let t = sum[1].clone();
    ibz_mul(&mut sum[1], &t, &alg.p);
    ibz_mul(&mut prod, &a[0], &b[1]);
    let t = sum[1].clone();
    ibz_add(&mut sum[1], &t, &prod);
    ibz_mul(&mut prod, &a[1], &b[0]);
    let t = sum[1].clone();
    ibz_add(&mut sum[1], &t, &prod);

    // j component: a0 b2 + a2 b0 - a1 b3 + a3 b1
    ibz_mul(&mut prod, &a[0], &b[2]);
    let t = sum[2].clone();
    ibz_add(&mut sum[2], &t, &prod);
    ibz_mul(&mut prod, &a[2], &b[0]);
    let t = sum[2].clone();
    ibz_add(&mut sum[2], &t, &prod);
    ibz_mul(&mut prod, &a[1], &b[3]);
    let t = sum[2].clone();
    ibz_sub(&mut sum[2], &t, &prod);
    ibz_mul(&mut prod, &a[3], &b[1]);
    let t = sum[2].clone();
    ibz_add(&mut sum[2], &t, &prod);

    // ij component: a0 b3 + a3 b0 - a2 b1 + a1 b2
    ibz_mul(&mut prod, &a[0], &b[3]);
    let t = sum[3].clone();
    ibz_add(&mut sum[3], &t, &prod);
    ibz_mul(&mut prod, &a[3], &b[0]);
    let t = sum[3].clone();
    ibz_add(&mut sum[3], &t, &prod);
    ibz_mul(&mut prod, &a[2], &b[1]);
    let t = sum[3].clone();
    ibz_sub(&mut sum[3], &t, &prod);
    ibz_mul(&mut prod, &a[1], &b[2]);
    let t = sum[3].clone();
    ibz_add(&mut sum[3], &t, &prod);

    ibz_vec_4_copy(res, &sum);
}

/// Put `a` and `b` on a common denominator.
pub fn quat_alg_equal_denom(
    res_a: &mut QuatAlgElem,
    res_b: &mut QuatAlgElem,
    a: &QuatAlgElem,
    b: &QuatAlgElem,
) {
    let mut gcd = Ibz::new();
    let mut r = Ibz::new();
    ibz_gcd(&mut gcd, &a.denom, &b.denom);
    ibz_div(&mut res_a.denom, &mut r, &a.denom, &gcd);
    ibz_div(&mut res_b.denom, &mut r, &b.denom, &gcd);
    for i in 0..4 {
        ibz_mul(&mut res_a.coord[i], &a.coord[i], &res_b.denom);
        ibz_mul(&mut res_b.coord[i], &b.coord[i], &res_a.denom);
    }
    let t = res_a.denom.clone();
    ibz_mul(&mut res_a.denom, &t, &res_b.denom);
    ibz_mul(&mut res_b.denom, &res_a.denom, &gcd);
    let t = res_a.denom.clone();
    ibz_mul(&mut res_a.denom, &t, &gcd);
}

pub fn quat_alg_add(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem) {
    let mut ra = QuatAlgElem::default();
    let mut rb = QuatAlgElem::default();
    quat_alg_equal_denom(&mut ra, &mut rb, a, b);
    ibz_copy(&mut res.denom, &ra.denom);
    ibz_vec_4_add(&mut res.coord, &ra.coord, &rb.coord);
}

pub fn quat_alg_sub(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem) {
    let mut ra = QuatAlgElem::default();
    let mut rb = QuatAlgElem::default();
    quat_alg_equal_denom(&mut ra, &mut rb, a, b);
    ibz_copy(&mut res.denom, &ra.denom);
    ibz_vec_4_sub(&mut res.coord, &ra.coord, &rb.coord);
}

pub fn quat_alg_mul(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem, alg: &QuatAlg) {
    ibz_mul(&mut res.denom, &a.denom, &b.denom);
    quat_alg_coord_mul(&mut res.coord, &a.coord, &b.coord, alg);
}

pub fn quat_alg_conj(conj: &mut QuatAlgElem, x: &QuatAlgElem) {
    ibz_copy(&mut conj.denom, &x.denom);
    ibz_copy(&mut conj.coord[0], &x.coord[0]);
    ibz_neg(&mut conj.coord[1], &x.coord[1]);
    ibz_neg(&mut conj.coord[2], &x.coord[2]);
    ibz_neg(&mut conj.coord[3], &x.coord[3]);
}

/// Reduced norm `Nrd(a) = a · ā` as a fraction `num/denom` (positive, gcd-reduced).
pub fn quat_alg_norm(res_num: &mut Ibz, res_denom: &mut Ibz, a: &QuatAlgElem, alg: &QuatAlg) {
    let mut r = Ibz::new();
    let mut g = Ibz::new();
    let mut norm = QuatAlgElem::default();
    quat_alg_conj(&mut norm, a);
    let conj = norm.clone();
    quat_alg_mul(&mut norm, a, &conj, alg);
    ibz_gcd(&mut g, &norm.coord[0], &norm.denom);
    ibz_div(res_num, &mut r, &norm.coord[0], &g);
    ibz_div(res_denom, &mut r, &norm.denom, &g);
    let t = res_denom.clone();
    ibz_abs(res_denom, &t);
    let t = res_num.clone();
    ibz_abs(res_num, &t);
    debug_assert!(ibz_cmp(res_denom, ibz_const_zero()) > 0);
}

pub fn quat_alg_scalar(elem: &mut QuatAlgElem, num: &Ibz, denom: &Ibz) {
    ibz_copy(&mut elem.denom, denom);
    ibz_copy(&mut elem.coord[0], num);
    ibz_set(&mut elem.coord[1], 0);
    ibz_set(&mut elem.coord[2], 0);
    ibz_set(&mut elem.coord[3], 0);
}

/// Normalize so `gcd(denom, content(coord)) = 1` and `denom > 0`.
pub fn quat_alg_normalize(x: &mut QuatAlgElem) {
    let mut gcd = Ibz::new();
    let mut r = Ibz::new();
    ibz_vec_4_content(&mut gcd, &x.coord);
    let g = gcd.clone();
    ibz_gcd(&mut gcd, &g, &x.denom);
    let d = x.denom.clone();
    ibz_div(&mut x.denom, &mut r, &d, &gcd);
    let c = x.coord.clone();
    ibz_vec_4_scalar_div(&mut x.coord, &gcd, &c);
    let sign = Ibz::from(2 * (ibz_cmp(ibz_const_zero(), &x.denom) < 0) as i32 - 1);
    let c = x.coord.clone();
    ibz_vec_4_scalar_mul(&mut x.coord, &sign, &c);
    let d = x.denom.clone();
    ibz_mul(&mut x.denom, &sign, &d);
}

pub fn quat_alg_elem_is_zero(x: &QuatAlgElem) -> i32 {
    ibz_vec_4_is_zero(&x.coord)
}

pub fn quat_alg_elem_equal(a: &QuatAlgElem, b: &QuatAlgElem) -> i32 {
    let mut diff = QuatAlgElem::default();
    quat_alg_sub(&mut diff, a, b);
    quat_alg_elem_is_zero(&diff)
}

pub fn quat_alg_elem_set(elem: &mut QuatAlgElem, denom: i32, c0: i32, c1: i32, c2: i32, c3: i32) {
    ibz_vec_4_set(&mut elem.coord, c0, c1, c2, c3);
    ibz_set(&mut elem.denom, denom);
}

pub fn quat_alg_elem_copy(dst: &mut QuatAlgElem, src: &QuatAlgElem) {
    ibz_copy(&mut dst.denom, &src.denom);
    ibz_vec_4_copy(&mut dst.coord, &src.coord);
}

pub fn quat_alg_elem_copy_ibz(
    elem: &mut QuatAlgElem,
    denom: &Ibz,
    c0: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
    c3: &Ibz,
) {
    ibz_vec_4_copy_ibz(&mut elem.coord, c0, c1, c2, c3);
    ibz_copy(&mut elem.denom, denom);
}

pub fn quat_alg_elem_mul_by_scalar(res: &mut QuatAlgElem, scalar: &Ibz, elem: &QuatAlgElem) {
    for i in 0..4 {
        ibz_mul(&mut res.coord[i], &elem.coord[i], scalar);
    }
    ibz_copy(&mut res.denom, &elem.denom);
}

// `quat_alg_make_primitive` depends on `quat_lattice_contains` (lattice.c).
// Ported in the lattice module.

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }
    fn elem(d: i32, c: [i32; 4]) -> QuatAlgElem {
        let mut e = QuatAlgElem::default();
        quat_alg_elem_set(&mut e, d, c[0], c[1], c[2], c[3]);
        e
    }

    #[test]
    fn coord_mul_basic() {
        let alg = QuatAlg::from_ui(7);
        // i * j = k
        let i = [z(0), z(1), z(0), z(0)];
        let j = [z(0), z(0), z(1), z(0)];
        let mut r = ibz_vec_4_init();
        quat_alg_coord_mul(&mut r, &i, &j, &alg);
        assert_eq!(r, [z(0), z(0), z(0), z(1)]);
        // j * i = -k
        quat_alg_coord_mul(&mut r, &j, &i, &alg);
        assert_eq!(r, [z(0), z(0), z(0), z(-1)]);
        // i² = -1
        quat_alg_coord_mul(&mut r, &i, &i, &alg);
        assert_eq!(r, [z(-1), z(0), z(0), z(0)]);
        // j² = -p
        quat_alg_coord_mul(&mut r, &j, &j, &alg);
        assert_eq!(r, [z(-7), z(0), z(0), z(0)]);
        // k² = (ij)(ij) = i(ji)j = i(-ij)j = -i²j² = -(-1)(-p) = -p
        let k = [z(0), z(0), z(0), z(1)];
        quat_alg_coord_mul(&mut r, &k, &k, &alg);
        assert_eq!(r, [z(-7), z(0), z(0), z(0)]);
    }

    #[test]
    fn norm_is_multiplicative() {
        let alg = QuatAlg::from_ui(11);
        let a = elem(1, [1, 2, 3, 4]);
        let b = elem(1, [5, -1, 2, 0]);
        let (mut na_n, mut na_d) = (Ibz::new(), Ibz::new());
        let (mut nb_n, mut nb_d) = (Ibz::new(), Ibz::new());
        quat_alg_norm(&mut na_n, &mut na_d, &a, &alg);
        quat_alg_norm(&mut nb_n, &mut nb_d, &b, &alg);
        let mut ab = QuatAlgElem::default();
        quat_alg_mul(&mut ab, &a, &b, &alg);
        let (mut nab_n, mut nab_d) = (Ibz::new(), Ibz::new());
        quat_alg_norm(&mut nab_n, &mut nab_d, &ab, &alg);
        assert_eq!(na_d, 1);
        assert_eq!(nb_d, 1);
        assert_eq!(nab_d, 1);
        assert_eq!(nab_n, na_n.clone() * nb_n.clone());
        // Direct: Nrd(a) = a0² + a1² + p(a2² + a3²) for denom 1.
        let expect = 1 + 4 + 11 * (9 + 16);
        assert_eq!(na_n, expect);
    }

    #[test]
    fn conj_and_norm_real() {
        let alg = QuatAlg::from_ui(3);
        let x = elem(1, [2, 1, 1, 1]);
        let mut xb = QuatAlgElem::default();
        quat_alg_conj(&mut xb, &x);
        let mut xxb = QuatAlgElem::default();
        quat_alg_mul(&mut xxb, &x, &xb, &alg);
        // x·x̄ is real (i,j,k components zero).
        assert_eq!(xxb.coord[1], 0);
        assert_eq!(xxb.coord[2], 0);
        assert_eq!(xxb.coord[3], 0);
    }

    #[test]
    fn add_sub_with_denoms() {
        let a = elem(2, [1, 0, 0, 0]); // 1/2
        let b = elem(3, [1, 0, 0, 0]); // 1/3
        let mut s = QuatAlgElem::default();
        quat_alg_add(&mut s, &a, &b);
        quat_alg_normalize(&mut s);
        assert_eq!(s.denom, 6);
        assert_eq!(s.coord[0], 5);
    }

    #[test]
    fn normalize_sign_and_gcd() {
        let mut x = elem(-4, [2, 6, 10, 14]);
        quat_alg_normalize(&mut x);
        assert_eq!(x.denom, 2);
        assert_eq!(x.coord, [z(-1), z(-3), z(-5), z(-7)]);
    }

    #[test]
    fn equal() {
        let a = elem(2, [4, 6, 0, 0]);
        let b = elem(1, [2, 3, 0, 0]);
        assert_eq!(quat_alg_elem_equal(&a, &b), 1);
    }
}
