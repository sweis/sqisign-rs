// SPDX-License-Identifier: Apache-2.0
//! 4-vector and 4×4 matrix operations over `Ibz`, ported from `dim4.c`.

use super::dim2::ibz_mat_2x2_det_from_ibz;
use super::intbig::*;
use super::types::*;

// ---------------------------------------------------------------------------
// IbzVec4

pub fn ibz_vec_4_set(vec: &mut IbzVec4, c0: i32, c1: i32, c2: i32, c3: i32) {
    ibz_set(&mut vec[0], c0);
    ibz_set(&mut vec[1], c1);
    ibz_set(&mut vec[2], c2);
    ibz_set(&mut vec[3], c3);
}

pub fn ibz_vec_4_copy(dst: &mut IbzVec4, src: &IbzVec4) {
    for i in 0..4 {
        ibz_copy(&mut dst[i], &src[i]);
    }
}

pub fn ibz_vec_4_copy_ibz(res: &mut IbzVec4, c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) {
    ibz_copy(&mut res[0], c0);
    ibz_copy(&mut res[1], c1);
    ibz_copy(&mut res[2], c2);
    ibz_copy(&mut res[3], c3);
}

pub fn ibz_vec_4_content(content: &mut Ibz, v: &IbzVec4) {
    ibz_gcd(content, &v[0], &v[1]);
    let t = content.clone();
    ibz_gcd(content, &v[2], &t);
    let t = content.clone();
    ibz_gcd(content, &v[3], &t);
}

pub fn ibz_vec_4_negate(neg: &mut IbzVec4, vec: &IbzVec4) {
    for i in 0..4 {
        ibz_neg(&mut neg[i], &vec[i]);
    }
}

pub fn ibz_vec_4_add(res: &mut IbzVec4, a: &IbzVec4, b: &IbzVec4) {
    for i in 0..4 {
        ibz_add(&mut res[i], &a[i], &b[i]);
    }
}

pub fn ibz_vec_4_sub(res: &mut IbzVec4, a: &IbzVec4, b: &IbzVec4) {
    for i in 0..4 {
        ibz_sub(&mut res[i], &a[i], &b[i]);
    }
}

pub fn ibz_vec_4_is_zero(x: &IbzVec4) -> i32 {
    x.iter().all(|v| ibz_is_zero(v) != 0) as i32
}

pub fn ibz_vec_4_linear_combination(
    lc: &mut IbzVec4,
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
) {
    let mut prod = Ibz::new();
    let mut sums = ibz_vec_4_init();
    for i in 0..4 {
        ibz_mul(&mut sums[i], coeff_a, &vec_a[i]);
        ibz_mul(&mut prod, coeff_b, &vec_b[i]);
        let t = sums[i].clone();
        ibz_add(&mut sums[i], &t, &prod);
    }
    ibz_vec_4_copy(lc, &sums);
}

pub fn ibz_vec_4_scalar_mul(prod: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4) {
    for i in 0..4 {
        ibz_mul(&mut prod[i], &vec[i], scalar);
    }
}

pub fn ibz_vec_4_scalar_div(quot: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4) -> i32 {
    let mut r = Ibz::new();
    let mut ok = 1;
    for i in 0..4 {
        ibz_div(&mut quot[i], &mut r, &vec[i], scalar);
        ok &= ibz_is_zero(&r);
    }
    ok
}

// ---------------------------------------------------------------------------
// IbzMat4x4

pub fn ibz_mat_4x4_mul(res: &mut IbzMat4x4, a: &IbzMat4x4, b: &IbzMat4x4) {
    let mut prod = Ibz::new();
    let mut mat = ibz_mat_4x4_init();
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut mat[i][j], 0);
            for k in 0..4 {
                ibz_mul(&mut prod, &a[i][k], &b[k][j]);
                let t = mat[i][j].clone();
                ibz_add(&mut mat[i][j], &t, &prod);
            }
        }
    }
    ibz_mat_4x4_copy(res, &mat);
}

pub fn ibz_mat_4x4_copy(dst: &mut IbzMat4x4, src: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_copy(&mut dst[i][j], &src[i][j]);
        }
    }
}

pub fn ibz_mat_4x4_negate(neg: &mut IbzMat4x4, mat: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_neg(&mut neg[i][j], &mat[i][j]);
        }
    }
}

pub fn ibz_mat_4x4_transpose(transposed: &mut IbzMat4x4, mat: &IbzMat4x4) {
    let mut work = ibz_mat_4x4_init();
    for i in 0..4 {
        for j in 0..4 {
            ibz_copy(&mut work[i][j], &mat[j][i]);
        }
    }
    ibz_mat_4x4_copy(transposed, &work);
}

pub fn ibz_mat_4x4_zero(m: &mut IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut m[i][j], 0);
        }
    }
}

pub fn ibz_mat_4x4_identity(id: &mut IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut id[i][j], (i == j) as i32);
        }
    }
}

pub fn ibz_mat_4x4_is_identity(mat: &IbzMat4x4) -> i32 {
    let mut res = true;
    for i in 0..4 {
        for j in 0..4 {
            res = res && ((ibz_is_one(&mat[i][j]) != 0) == (i == j));
        }
    }
    res as i32
}

pub fn ibz_mat_4x4_equal(a: &IbzMat4x4, b: &IbzMat4x4) -> i32 {
    for i in 0..4 {
        for j in 0..4 {
            if ibz_cmp(&a[i][j], &b[i][j]) != 0 {
                return 0;
            }
        }
    }
    1
}

pub fn ibz_mat_4x4_scalar_mul(prod: &mut IbzMat4x4, scalar: &Ibz, mat: &IbzMat4x4) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_mul(&mut prod[i][j], &mat[i][j], scalar);
        }
    }
}

pub fn ibz_mat_4x4_gcd(gcd: &mut Ibz, mat: &IbzMat4x4) {
    let mut d = mat[0][0].clone();
    for i in 0..4 {
        for j in 0..4 {
            let t = d.clone();
            ibz_gcd(&mut d, &t, &mat[i][j]);
        }
    }
    ibz_copy(gcd, &d);
}

pub fn ibz_mat_4x4_scalar_div(quot: &mut IbzMat4x4, scalar: &Ibz, mat: &IbzMat4x4) -> i32 {
    let mut r = Ibz::new();
    let mut ok = 1;
    for i in 0..4 {
        for j in 0..4 {
            ibz_div(&mut quot[i][j], &mut r, &mat[i][j], scalar);
            ok &= ibz_is_zero(&r);
        }
    }
    ok
}

// 4×4 inversion via Laplace expansion (geometrictools method).

fn coeff_pmp(out: &mut Ibz, a1: &Ibz, a2: &Ibz, b1: &Ibz, b2: &Ibz, c1: &Ibz, c2: &Ibz) {
    let mut prod = Ibz::new();
    let mut sum = Ibz::new();
    ibz_mul(&mut sum, a1, a2);
    ibz_mul(&mut prod, b1, b2);
    let t = sum.clone();
    ibz_sub(&mut sum, &t, &prod);
    ibz_mul(&mut prod, c1, c2);
    ibz_add(out, &sum, &prod);
}

fn coeff_mpm(out: &mut Ibz, a1: &Ibz, a2: &Ibz, b1: &Ibz, b2: &Ibz, c1: &Ibz, c2: &Ibz) {
    let mut prod = Ibz::new();
    let mut sum = Ibz::new();
    ibz_mul(&mut sum, b1, b2);
    ibz_mul(&mut prod, a1, a2);
    let t = sum.clone();
    ibz_sub(&mut sum, &t, &prod);
    ibz_mul(&mut prod, c1, c2);
    ibz_sub(out, &sum, &prod);
}

pub fn ibz_inv_dim4_make_coeff_pmp(
    out: &mut Ibz,
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) {
    coeff_pmp(out, a1, a2, b1, b2, c1, c2);
}
pub fn ibz_inv_dim4_make_coeff_mpm(
    out: &mut Ibz,
    a1: &Ibz,
    a2: &Ibz,
    b1: &Ibz,
    b2: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
) {
    coeff_mpm(out, a1, a2, b1, b2, c1, c2);
}

/// Computes the adjugate of `mat` into `inv` (so `mat * inv = det * I`) and
/// writes the determinant. Returns 1 if `det != 0`.
pub fn ibz_mat_4x4_inv_with_det_as_denom(
    inv: Option<&mut IbzMat4x4>,
    det: &mut Ibz,
    mat: &IbzMat4x4,
) -> i32 {
    let mut prod = Ibz::new();
    let mut work_det = Ibz::new();
    let mut work = ibz_mat_4x4_init();
    let mut s: [Ibz; 6] = core::array::from_fn(|_| Ibz::new());
    let mut c: [Ibz; 6] = core::array::from_fn(|_| Ibz::new());

    // 2×2 minors.
    for i in 0..3 {
        ibz_mat_2x2_det_from_ibz(
            &mut s[i],
            &mat[0][0],
            &mat[0][i + 1],
            &mat[1][0],
            &mat[1][i + 1],
        );
        ibz_mat_2x2_det_from_ibz(
            &mut c[i],
            &mat[2][0],
            &mat[2][i + 1],
            &mat[3][0],
            &mat[3][i + 1],
        );
    }
    for i in 0..2 {
        ibz_mat_2x2_det_from_ibz(
            &mut s[3 + i],
            &mat[0][1],
            &mat[0][2 + i],
            &mat[1][1],
            &mat[1][2 + i],
        );
        ibz_mat_2x2_det_from_ibz(
            &mut c[3 + i],
            &mat[2][1],
            &mat[2][2 + i],
            &mat[3][1],
            &mat[3][2 + i],
        );
    }
    ibz_mat_2x2_det_from_ibz(&mut s[5], &mat[0][2], &mat[0][3], &mat[1][2], &mat[1][3]);
    ibz_mat_2x2_det_from_ibz(&mut c[5], &mat[2][2], &mat[2][3], &mat[3][2], &mat[3][3]);

    // det = Σ ±s[i]*c[5-i].
    ibz_set(&mut work_det, 0);
    for i in 0..6 {
        ibz_mul(&mut prod, &s[i], &c[5 - i]);
        let t = work_det.clone();
        if i != 1 && i != 4 {
            ibz_add(&mut work_det, &t, &prod);
        } else {
            ibz_sub(&mut work_det, &t, &prod);
        }
    }

    // Transposed adjugate.
    let b = |x: bool| x as usize;
    for j in 0..4usize {
        for k in 0..2usize {
            let row = 1 - k;
            let cidx0 = 6 - j - b(j == 0);
            let cidx1 = 4 - j - b(j == 1);
            let cidx2 = 3 - j - b(j == 1) - b(j == 2);
            let col0 = b(j == 0);
            let col1 = 2 - b(j > 1);
            let col2 = 3 - b(j == 3);
            if (k + j + 1) % 2 == 1 {
                coeff_pmp(
                    &mut work[j][k],
                    &mat[row][col0],
                    &c[cidx0],
                    &mat[row][col1],
                    &c[cidx1],
                    &mat[row][col2],
                    &c[cidx2],
                );
            } else {
                coeff_mpm(
                    &mut work[j][k],
                    &mat[row][col0],
                    &c[cidx0],
                    &mat[row][col1],
                    &c[cidx1],
                    &mat[row][col2],
                    &c[cidx2],
                );
            }
        }
        for k in 2..4usize {
            let row = 3 - b(k == 3);
            let sidx0 = 6 - j - b(j == 0);
            let sidx1 = 4 - j - b(j == 1);
            let sidx2 = 3 - j - b(j == 1) - b(j == 2);
            let col0 = b(j == 0);
            let col1 = 2 - b(j > 1);
            let col2 = 3 - b(j == 3);
            if (k + j + 1) % 2 == 1 {
                coeff_pmp(
                    &mut work[j][k],
                    &mat[row][col0],
                    &s[sidx0],
                    &mat[row][col1],
                    &s[sidx1],
                    &mat[row][col2],
                    &s[sidx2],
                );
            } else {
                coeff_mpm(
                    &mut work[j][k],
                    &mat[row][col0],
                    &s[sidx0],
                    &mat[row][col1],
                    &s[sidx1],
                    &mat[row][col2],
                    &s[sidx2],
                );
            }
        }
    }

    if let Some(inv) = inv {
        ibz_set(&mut prod, (!work_det.is_zero()) as i32);
        ibz_mat_4x4_scalar_mul(inv, &prod, &work);
    }
    ibz_copy(det, &work_det);
    (!det.is_zero()) as i32
}

/// `res = mat * vec`.
pub fn ibz_mat_4x4_eval(res: &mut IbzVec4, mat: &IbzMat4x4, vec: &IbzVec4) {
    let mut prod = Ibz::new();
    let mut sum = ibz_vec_4_init();
    for i in 0..4 {
        for j in 0..4 {
            ibz_mul(&mut prod, &mat[i][j], &vec[j]);
            let t = sum[i].clone();
            ibz_add(&mut sum[i], &t, &prod);
        }
    }
    ibz_vec_4_copy(res, &sum);
}

/// `res = vecᵀ * mat`.
pub fn ibz_mat_4x4_eval_t(res: &mut IbzVec4, vec: &IbzVec4, mat: &IbzMat4x4) {
    let mut prod = Ibz::new();
    let mut sum = ibz_vec_4_init();
    for i in 0..4 {
        for j in 0..4 {
            ibz_mul(&mut prod, &mat[j][i], &vec[j]);
            let t = sum[i].clone();
            ibz_add(&mut sum[i], &t, &prod);
        }
    }
    ibz_vec_4_copy(res, &sum);
}

/// Quadratic form evaluation: `res = coordᵀ · qf · coord`.
pub fn quat_qf_eval(res: &mut Ibz, qf: &IbzMat4x4, coord: &IbzVec4) {
    let mut prod = Ibz::new();
    let mut sum = ibz_vec_4_init();
    ibz_mat_4x4_eval(&mut sum, qf, coord);
    for i in 0..4 {
        ibz_mul(&mut prod, &sum[i], &coord[i]);
        if i > 0 {
            let t = sum[0].clone();
            ibz_add(&mut sum[0], &t, &prod);
        } else {
            ibz_copy(&mut sum[0], &prod);
        }
    }
    ibz_copy(res, &sum[0]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }

    fn mat_from(rows: [[i64; 4]; 4]) -> IbzMat4x4 {
        core::array::from_fn(|i| core::array::from_fn(|j| z(rows[i][j])))
    }

    #[test]
    fn vec4_ops() {
        let a: IbzVec4 = [z(1), z(2), z(3), z(4)];
        let b: IbzVec4 = [z(5), z(6), z(7), z(8)];
        let mut r = ibz_vec_4_init();
        ibz_vec_4_add(&mut r, &a, &b);
        assert_eq!(r[3], 12);
        let mut g = Ibz::new();
        ibz_vec_4_content(&mut g, &[z(6), z(10), z(14), z(22)]);
        assert_eq!(g, 2);
        ibz_vec_4_linear_combination(&mut r, &z(2), &a, &z(3), &b);
        assert_eq!(r[0], 17);
        assert_eq!(r[3], 32);
    }

    #[test]
    fn mat4x4_mul_identity() {
        let mut id = ibz_mat_4x4_init();
        ibz_mat_4x4_identity(&mut id);
        assert_eq!(ibz_mat_4x4_is_identity(&id), 1);
        let m = mat_from([
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
        ]);
        let mut r = ibz_mat_4x4_init();
        ibz_mat_4x4_mul(&mut r, &id, &m);
        assert_eq!(ibz_mat_4x4_equal(&r, &m), 1);
        ibz_mat_4x4_mul(&mut r, &m, &id);
        assert_eq!(ibz_mat_4x4_equal(&r, &m), 1);
    }

    #[test]
    fn mat4x4_inv_with_det() {
        // A matrix with known det.
        let m = mat_from([[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 5, 0], [0, 0, 0, 7]]);
        let mut inv = ibz_mat_4x4_init();
        let mut det = Ibz::new();
        assert_eq!(
            ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), &mut det, &m),
            1
        );
        assert_eq!(det, 210);
        // m * inv should equal det * I.
        let mut prod = ibz_mat_4x4_init();
        ibz_mat_4x4_mul(&mut prod, &m, &inv);
        for i in 0..4 {
            for j in 0..4 {
                let expect = if i == j { 210 } else { 0 };
                assert_eq!(prod[i][j], expect, "[{i}][{j}]");
            }
        }
        // A non-diagonal example.
        let m2 = mat_from([[1, 2, 0, 1], [0, 1, 3, 0], [2, 0, 1, 1], [1, 1, 0, 2]]);
        let mut det2 = Ibz::new();
        let mut inv2 = ibz_mat_4x4_init();
        assert_eq!(
            ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv2), &mut det2, &m2),
            1
        );
        let mut prod2 = ibz_mat_4x4_init();
        ibz_mat_4x4_mul(&mut prod2, &m2, &inv2);
        for i in 0..4 {
            for j in 0..4 {
                let expect = if i == j { det2.clone() } else { z(0) };
                assert_eq!(prod2[i][j], expect, "m2[{i}][{j}]");
            }
        }
        // Singular: det = 0.
        let sing = mat_from([[1, 2, 3, 4], [2, 4, 6, 8], [0, 1, 0, 1], [1, 0, 1, 0]]);
        let mut det3 = Ibz::new();
        assert_eq!(ibz_mat_4x4_inv_with_det_as_denom(None, &mut det3, &sing), 0);
        assert_eq!(det3, 0);
    }

    #[test]
    fn qf_eval() {
        // Identity QF: vᵀ·I·v = Σ vᵢ².
        let mut id = ibz_mat_4x4_init();
        ibz_mat_4x4_identity(&mut id);
        let v: IbzVec4 = [z(1), z(2), z(3), z(4)];
        let mut r = Ibz::new();
        quat_qf_eval(&mut r, &id, &v);
        assert_eq!(r, 30);
    }
}
