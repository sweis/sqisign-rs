// SPDX-License-Identifier: Apache-2.0
//! L²-LLL reduction (Nguyen–Stehlé) on 4-dimensional lattices.
//! Ported from `lll/{l2.c, lll_applications.c, lll_verification.c}`.

use super::dim4::*;
use super::dpe::Dpe;
use super::intbig::*;
use super::lattice::*;
use super::rationals::*;
use super::types::*;

pub const DELTABAR: f64 = 0.995;
pub const DELTA_NUM: i32 = 99;
pub const DELTA_DENOM: i32 = 100;
pub const ETABAR: f64 = 0.505;
pub const EPSILON_NUM: i32 = 1;
pub const EPSILON_DENOM: i32 = 100;

/// Symmetric-matrix accessor (lower-triangular storage).
#[inline]
fn sym_get(g: &IbzMat4x4, i: usize, j: usize) -> &Ibz {
    if i < j { &g[j][i] } else { &g[i][j] }
}
#[inline]
fn sym_get_mut(g: &mut IbzMat4x4, i: usize, j: usize) -> &mut Ibz {
    if i < j { &mut g[j][i] } else { &mut g[i][j] }
}

/// In-place L² reduction. `g` is the (symmetric) Gram matrix of `basis`;
/// both are updated.
pub fn quat_lll_core(g: &mut IbzMat4x4, basis: &mut IbzMat4x4) {
    let mut delta_bar = Dpe::ZERO;
    delta_bar.set_d(DELTABAR);

    let mut r = [[Dpe::ZERO; 4]; 4];
    let mut u = [[Dpe::ZERO; 4]; 4];
    let mut lovasz = [Dpe::ZERO; 4];

    let mut xf = Dpe::ZERO;
    let mut tmpf = Dpe::ZERO;
    let mut x = Ibz::new();
    let mut tmpi = Ibz::new();

    r[0][0].set_z(&g[0][0]);
    let mut kappa: usize = 1;
    while kappa < 4 {
        // Size-reduce b_κ.
        let mut done = false;
        while !done {
            // Recompute the κ-th row of the Cholesky factorization.
            for j in 0..=kappa {
                r[kappa][j].set_z(&g[kappa][j]);
                for k in 0..j {
                    let rk = r[kappa][k];
                    let uk = u[j][k];
                    tmpf.mul(&rk, &uk);
                    let cur = r[kappa][j];
                    r[kappa][j].sub(&cur, &tmpf);
                }
                if j < kappa {
                    let num = r[kappa][j];
                    let den = r[j][j];
                    u[kappa][j].div(&num, &den);
                }
            }

            done = true;
            for i in (0..kappa).rev() {
                if u[kappa][i].cmp_d(ETABAR) > 0 || u[kappa][i].cmp_d(-ETABAR) < 0 {
                    done = false;
                    xf.set(&u[kappa][i]);
                    let xf0 = xf;
                    xf.round(&xf0);
                    xf.get_z(&mut x);
                    // b_κ ← b_κ − X·b_i
                    for j in 0..4 {
                        ibz_mul(&mut tmpi, &x, &basis[j][i]);
                        let t = basis[j][kappa].clone();
                        ibz_sub(&mut basis[j][kappa], &t, &tmpi);
                    }
                    // Update lower half of Gram matrix.
                    ibz_mul(&mut tmpi, &x, &g[kappa][i]);
                    let t = g[kappa][kappa].clone();
                    ibz_sub(&mut g[kappa][kappa], &t, &tmpi);
                    for j in 0..4 {
                        ibz_mul(&mut tmpi, &x, sym_get(g, i, j));
                        let t = sym_get(g, kappa, j).clone();
                        ibz_sub(sym_get_mut(g, kappa, j), &t, &tmpi);
                    }
                    // Update u[κ][0..i).
                    for j in 0..i {
                        let uij = u[i][j];
                        tmpf.mul(&xf, &uij);
                        let cur = u[kappa][j];
                        u[kappa][j].sub(&cur, &tmpf);
                    }
                }
            }
        }

        // Lovász conditions.
        lovasz[0].set_z(&g[kappa][kappa]);
        for i in 1..kappa {
            let uk = u[kappa][i - 1];
            let rk = r[kappa][i - 1];
            tmpf.mul(&uk, &rk);
            let prev = lovasz[i - 1];
            lovasz[i].sub(&prev, &tmpf);
        }
        let mut swap = kappa;
        while swap > 0 {
            let rss = r[swap - 1][swap - 1];
            tmpf.mul(&delta_bar, &rss);
            if tmpf.cmp(&lovasz[swap - 1]) < 0 {
                break;
            }
            swap -= 1;
        }

        if kappa != swap {
            for jc in (swap + 1..=kappa).rev() {
                for ir in 0..4 {
                    basis[ir].swap(jc, jc - 1);
                    if ir == jc - 1 {
                        // Swap diagonal entries g[ir][ir] and g[jc][jc].
                        let t = g[ir][ir].clone();
                        let s = g[jc][jc].clone();
                        g[ir][ir] = s;
                        g[jc][jc] = t;
                    } else if ir != jc {
                        let a = sym_get(g, ir, jc).clone();
                        let b = sym_get(g, ir, jc - 1).clone();
                        *sym_get_mut(g, ir, jc) = b;
                        *sym_get_mut(g, ir, jc - 1) = a;
                    }
                }
            }
            for i in 0..swap {
                u[swap][i] = u[kappa][i];
                r[swap][i] = r[kappa][i];
            }
            r[swap][swap] = lovasz[swap];
            kappa = swap;
        }

        kappa += 1;
    }

    #[cfg(debug_assertions)]
    {
        let mut deltabar_dpe = Dpe::ZERO;
        deltabar_dpe.set_d(DELTABAR);
        for i in 0..4 {
            for j in 0..i {
                let mut a = Dpe::ZERO;
                a.abs(&u[i][j]);
                debug_assert!(a.cmp_d(ETABAR) <= 0, "size-reduce fail i={i} j={j}");
            }
        }
        for i in 1..4 {
            let mut t = Dpe::ZERO;
            t.mul(&u[i][i - 1], &u[i][i - 1]);
            let mut s = Dpe::ZERO;
            s.sub(&deltabar_dpe, &t);
            let mut p = Dpe::ZERO;
            p.mul(&s, &r[i - 1][i - 1]);
            debug_assert!(p.cmp(&r[i][i]) <= 0, "Lovász fail i={i}");
        }
    }

    // Mirror lower half to upper half.
    for i in 0..4 {
        for j in (i + 1)..4 {
            let t = g[j][i].clone();
            ibz_copy(&mut g[i][j], &t);
        }
    }
}

pub fn quat_lattice_lll(red: &mut IbzMat4x4, lattice: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut g = ibz_mat_4x4_init();
    quat_lattice_gram(&mut g, lattice, alg);
    ibz_mat_4x4_copy(red, &lattice.basis);
    quat_lll_core(&mut g, red);
    0
}

// ---------------------------------------------------------------------------
// LLL verification (rationals; debug/test only).

pub fn quat_lll_set_ibq_parameters(delta: &mut Ibq, eta: &mut Ibq) {
    let mut num = Ibz::new();
    let mut denom = Ibz::new();
    ibq_set(delta, ibz_const_one(), ibz_const_two());
    ibz_set(&mut num, EPSILON_NUM);
    ibz_set(&mut denom, EPSILON_DENOM);
    ibq_set(eta, &num, &denom);
    let e = eta.clone();
    let d = delta.clone();
    ibq_add(eta, &e, &d);
    ibz_set(&mut num, DELTA_NUM);
    ibz_set(&mut denom, DELTA_DENOM);
    ibq_set(delta, &num, &denom);
}

pub fn ibq_vec_4_copy_ibz(vec: &mut IbqVec4, c0: &Ibz, c1: &Ibz, c2: &Ibz, c3: &Ibz) {
    let one = Ibz::from(1);
    ibq_set(&mut vec[0], c0, &one);
    ibq_set(&mut vec[1], c1, &one);
    ibq_set(&mut vec[2], c2, &one);
    ibq_set(&mut vec[3], c3, &one);
}

pub fn quat_lll_bilinear(b: &mut Ibq, vec0: &IbqVec4, vec1: &IbqVec4, q: &Ibz) {
    let one = Ibz::from(1);
    let mut sum = ibq_init();
    let mut prod = ibq_init();
    let mut norm_q = ibq_init();
    ibq_set(&mut norm_q, q, &one);

    ibq_mul(&mut sum, &vec0[0], &vec1[0]);
    ibq_mul(&mut prod, &vec0[1], &vec1[1]);
    let s = sum.clone();
    ibq_add(&mut sum, &s, &prod);
    ibq_mul(&mut prod, &vec0[2], &vec1[2]);
    let p = prod.clone();
    ibq_mul(&mut prod, &p, &norm_q);
    let s = sum.clone();
    ibq_add(&mut sum, &s, &prod);
    ibq_mul(&mut prod, &vec0[3], &vec1[3]);
    let p = prod.clone();
    ibq_mul(&mut prod, &p, &norm_q);
    ibq_add(b, &sum, &prod);
}

pub fn quat_lll_gram_schmidt_transposed_with_ibq(
    orthogonalised_transposed: &mut IbqMat4x4,
    mat: &IbzMat4x4,
    q: &Ibz,
) {
    let mut work = ibq_mat_4x4_init();
    let mut vec = ibq_vec_4_init();
    let mut norm = ibq_init();
    let mut b = ibq_init();
    let mut coeff = ibq_init();
    let mut prod = ibq_init();
    for i in 0..4 {
        ibq_vec_4_copy_ibz(&mut work[i], &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
    }
    for i in 0..4 {
        let wi = work[i].clone();
        quat_lll_bilinear(&mut norm, &wi, &wi, q);
        let n = norm.clone();
        ibq_inv(&mut norm, &n);
        for j in (i + 1)..4 {
            ibq_vec_4_copy_ibz(&mut vec, &mat[0][j], &mat[1][j], &mat[2][j], &mat[3][j]);
            let wi = work[i].clone();
            quat_lll_bilinear(&mut b, &wi, &vec, q);
            ibq_mul(&mut coeff, &norm, &b);
            for k in 0..4 {
                ibq_mul(&mut prod, &coeff, &work[i][k]);
                let wjk = work[j][k].clone();
                ibq_sub(&mut work[j][k], &wjk, &prod);
            }
        }
    }
    for i in 0..4 {
        for j in 0..4 {
            ibq_copy(&mut orthogonalised_transposed[i][j], &work[i][j]);
        }
    }
}

pub fn quat_lll_verify(mat: &IbzMat4x4, delta: &Ibq, eta: &Ibq, alg: &QuatAlg) -> i32 {
    let mut res = true;
    let mut ot = ibq_mat_4x4_init();
    let mut tmp_vec = ibq_vec_4_init();
    let mut div = ibq_init();
    let mut tmp = ibq_init();
    let mut mu = ibq_init();
    let mut norm = ibq_init();
    let mut b = ibq_init();

    quat_lll_gram_schmidt_transposed_with_ibq(&mut ot, mat, &alg.p);
    for i in 0..4 {
        for j in 0..i {
            ibq_vec_4_copy_ibz(&mut tmp_vec, &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
            quat_lll_bilinear(&mut b, &ot[j], &tmp_vec, &alg.p);
            quat_lll_bilinear(&mut norm, &ot[j], &ot[j], &alg.p);
            ibq_inv(&mut tmp, &norm);
            ibq_mul(&mut mu, &b, &tmp);
            let m = mu.clone();
            ibq_abs(&mut mu, &m);
            res = res && ibq_cmp(&mu, eta) <= 0;
        }
    }
    for i in 1..4 {
        ibq_vec_4_copy_ibz(&mut tmp_vec, &mat[0][i], &mat[1][i], &mat[2][i], &mat[3][i]);
        quat_lll_bilinear(&mut b, &ot[i - 1], &tmp_vec, &alg.p);
        quat_lll_bilinear(&mut norm, &ot[i - 1], &ot[i - 1], &alg.p);
        ibq_inv(&mut tmp, &norm);
        ibq_mul(&mut mu, &b, &tmp);
        ibq_mul(&mut tmp, &mu, &mu);
        let t = tmp.clone();
        ibq_sub(&mut mu, delta, &t);
        quat_lll_bilinear(&mut tmp, &ot[i], &ot[i], &alg.p);
        ibq_mul(&mut div, &norm, &mu);
        res = res && ibq_cmp(&tmp, &div) >= 0;
    }
    res as i32
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
    fn lll_reduces_and_verifies() {
        let alg = QuatAlg::from_ui(103);
        // A non-reduced basis: identity scaled by 5 plus an off-diagonal mess.
        let mut lat = QuatLattice::default();
        lat.basis = mat_from([[5, 11, 0, 0], [0, 5, 7, 0], [0, 0, 5, 9], [0, 0, 0, 5]]);
        let mut red = ibz_mat_4x4_init();
        quat_lattice_lll(&mut red, &lat, &alg);
        // Verify it's LLL-reduced for δ=99/100, η=51/100.
        let mut delta = ibq_init();
        let mut eta = ibq_init();
        quat_lll_set_ibq_parameters(&mut delta, &mut eta);
        assert_eq!(quat_lll_verify(&red, &delta, &eta, &alg), 1);
        // det preserved.
        let mut d0 = Ibz::new();
        let mut d1 = Ibz::new();
        ibz_mat_4x4_inv_with_det_as_denom(None, &mut d0, &lat.basis);
        ibz_mat_4x4_inv_with_det_as_denom(None, &mut d1, &red);
        assert_eq!(d0.clone().abs(), d1.clone().abs());
    }

    #[test]
    fn lll_on_o0() {
        // The standard maximal order O₀ is already short; LLL should still verify.
        let alg = QuatAlg::from_ui(7);
        let mut lat = QuatLattice::default();
        super::super::normeq::quat_lattice_o0_set(&mut lat);
        let mut red = ibz_mat_4x4_init();
        quat_lattice_lll(&mut red, &lat, &alg);
        let mut delta = ibq_init();
        let mut eta = ibq_init();
        quat_lll_set_ibq_parameters(&mut delta, &mut eta);
        assert_eq!(quat_lll_verify(&red, &delta, &eta, &alg), 1);
    }
}
