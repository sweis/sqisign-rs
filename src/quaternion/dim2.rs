// SPDX-License-Identifier: Apache-2.0
//! 2-vector and 2×2 matrix operations over `Ibz`, ported from `dim2.c`.

use super::intbig::*;
use super::types::*;

pub fn ibz_vec_2_set(vec: &mut IbzVec2, a0: i32, a1: i32) {
    ibz_set(&mut vec[0], a0);
    ibz_set(&mut vec[1], a1);
}

pub fn ibz_mat_2x2_set(mat: &mut IbzMat2x2, a00: i32, a01: i32, a10: i32, a11: i32) {
    ibz_set(&mut mat[0][0], a00);
    ibz_set(&mut mat[0][1], a01);
    ibz_set(&mut mat[1][0], a10);
    ibz_set(&mut mat[1][1], a11);
}

pub fn ibz_mat_2x2_copy(copy: &mut IbzMat2x2, src: &IbzMat2x2) {
    for i in 0..2 {
        for j in 0..2 {
            ibz_copy(&mut copy[i][j], &src[i][j]);
        }
    }
}

pub fn ibz_mat_2x2_add(sum: &mut IbzMat2x2, a: &IbzMat2x2, b: &IbzMat2x2) {
    for i in 0..2 {
        for j in 0..2 {
            ibz_add(&mut sum[i][j], &a[i][j], &b[i][j]);
        }
    }
}

/// `det = a11*a22 - a12*a21`.
pub fn ibz_mat_2x2_det_from_ibz(det: &mut Ibz, a11: &Ibz, a12: &Ibz, a21: &Ibz, a22: &Ibz) {
    let mut prod = Ibz::new();
    ibz_mul(&mut prod, a12, a21);
    ibz_mul(det, a11, a22);
    let d = det.clone();
    ibz_sub(det, &d, &prod);
}

/// `res = mat * vec`.
pub fn ibz_mat_2x2_eval(res: &mut IbzVec2, mat: &IbzMat2x2, vec: &IbzVec2) {
    let mut prod = Ibz::new();
    let mut mv = ibz_vec_2_init();
    ibz_mul(&mut mv[0], &mat[0][0], &vec[0]);
    ibz_mul(&mut prod, &mat[0][1], &vec[1]);
    let t = mv[0].clone();
    ibz_add(&mut mv[0], &t, &prod);
    ibz_mul(&mut mv[1], &mat[1][0], &vec[0]);
    ibz_mul(&mut prod, &mat[1][1], &vec[1]);
    let t = mv[1].clone();
    ibz_add(&mut mv[1], &t, &prod);
    ibz_copy(&mut res[0], &mv[0]);
    ibz_copy(&mut res[1], &mv[1]);
}

/// `prod = (mat_a * mat_b) mod m` (entrywise reduction during accumulation).
pub fn ibz_2x2_mul_mod(prod: &mut IbzMat2x2, mat_a: &IbzMat2x2, mat_b: &IbzMat2x2, m: &Ibz) {
    let mut mul = Ibz::new();
    let mut sums = ibz_mat_2x2_init();
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                ibz_mul(&mut mul, &mat_a[i][k], &mat_b[k][j]);
                let t = sums[i][j].clone();
                ibz_add(&mut sums[i][j], &t, &mul);
                let t = sums[i][j].clone();
                ibz_mod(&mut sums[i][j], &t, m);
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            ibz_copy(&mut prod[i][j], &sums[i][j]);
        }
    }
}

/// Modular inverse of a 2×2 matrix. Returns 1 if `det` is invertible mod `m`,
/// 0 otherwise (and writes the zero matrix).
pub fn ibz_mat_2x2_inv_mod(inv: &mut IbzMat2x2, mat: &IbzMat2x2, m: &Ibz) -> i32 {
    let mut det = Ibz::new();
    let mut prod = Ibz::new();
    ibz_mul(&mut det, &mat[0][0], &mat[1][1]);
    let t = det.clone();
    ibz_mod(&mut det, &t, m);
    ibz_mul(&mut prod, &mat[0][1], &mat[1][0]);
    let t = det.clone();
    ibz_sub(&mut det, &t, &prod);
    let t = det.clone();
    ibz_mod(&mut det, &t, m);
    let res = {
        let d = det.clone();
        ibz_invmod(&mut det, &d, m)
    };
    // Zero matrix if non-invertible determinant.
    ibz_set(&mut prod, res);
    let t = det.clone();
    ibz_mul(&mut det, &t, &prod);
    // Adjugate.
    ibz_copy(&mut prod, &mat[0][0]);
    ibz_copy(&mut inv[0][0], &mat[1][1]);
    ibz_copy(&mut inv[1][1], &prod);
    ibz_neg(&mut inv[1][0], &mat[1][0]);
    ibz_neg(&mut inv[0][1], &mat[0][1]);
    for i in 0..2 {
        for j in 0..2 {
            let t = inv[i][j].clone();
            ibz_mul(&mut inv[i][j], &t, &det);
            let t = inv[i][j].clone();
            ibz_mod(&mut inv[i][j], &t, m);
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }

    #[test]
    fn det2x2() {
        let mut d = Ibz::new();
        ibz_mat_2x2_det_from_ibz(&mut d, &z(3), &z(2), &z(5), &z(7));
        assert_eq!(d, 11);
    }

    #[test]
    fn eval_and_mul_mod() {
        let mut m = ibz_mat_2x2_init();
        ibz_mat_2x2_set(&mut m, 1, 2, 3, 4);
        let v: IbzVec2 = [z(5), z(6)];
        let mut r = ibz_vec_2_init();
        ibz_mat_2x2_eval(&mut r, &m, &v);
        assert_eq!(r[0], 17);
        assert_eq!(r[1], 39);

        let modn = z(7);
        let mut p = ibz_mat_2x2_init();
        ibz_2x2_mul_mod(&mut p, &m, &m, &modn);
        // [[1,2],[3,4]]² = [[7,10],[15,22]] → mod 7 = [[0,3],[1,1]]
        assert_eq!(p[0][0], 0);
        assert_eq!(p[0][1], 3);
        assert_eq!(p[1][0], 1);
        assert_eq!(p[1][1], 1);
    }

    #[test]
    fn inv_mod_roundtrip() {
        let mut m = ibz_mat_2x2_init();
        ibz_mat_2x2_set(&mut m, 1, 2, 3, 4); // det = -2, invertible mod 7
        let modn = z(7);
        let mut inv = ibz_mat_2x2_init();
        assert_eq!(ibz_mat_2x2_inv_mod(&mut inv, &m, &modn), 1);
        let mut prod = ibz_mat_2x2_init();
        ibz_2x2_mul_mod(&mut prod, &m, &inv, &modn);
        assert_eq!(prod[0][0], 1);
        assert_eq!(prod[0][1], 0);
        assert_eq!(prod[1][0], 0);
        assert_eq!(prod[1][1], 1);
        // Non-invertible: det = 0 mod 5
        ibz_mat_2x2_set(&mut m, 1, 2, 2, 4);
        assert_eq!(ibz_mat_2x2_inv_mod(&mut inv, &m, &z(5)), 0);
    }
}
