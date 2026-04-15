//! Modular Hermite Normal Form for 4×n matrices.
//! Ported from `hnf/{ibz_division.c, hnf_internal.c, hnf.c}`.

use super::dim4::*;
use super::intbig::*;
use super::types::*;

/// `res = x mod m`, but mapped to `(0, m]` instead of `[0, m)`.
pub fn ibz_mod_not_zero(res: &mut Ibz, x: &Ibz, modn: &Ibz) {
    let mut m = Ibz::default();
    let mut t = Ibz::default();
    ibz_mod(&mut m, x, modn);
    ibz_set(&mut t, ibz_is_zero(&m) as i32);
    let s = t.clone();
    ibz_mul(&mut t, &s, modn);
    ibz_add(res, &m, &t);
}

/// Centered remainder in `(-m/2, m/2]` (positive when on the boundary).
pub fn ibz_centered_mod(remainder: &mut Ibz, a: &Ibz, modn: &Ibz) {
    debug_assert!(ibz_cmp(modn, ibz_const_zero()) > 0);
    let mut tmp = Ibz::default();
    let mut d = Ibz::default();
    let mut t = Ibz::default();
    ibz_div_floor(&mut d, &mut tmp, modn, ibz_const_two());
    ibz_mod_not_zero(&mut tmp, a, modn);
    ibz_set(&mut t, (ibz_cmp(&tmp, &d) > 0) as i32);
    let s = t.clone();
    ibz_mul(&mut t, &s, modn);
    ibz_sub(remainder, &tmp, &t);
}

/// Extended GCD with `u != 0`, `d > 0`, and `u·x > 0` minimal.
/// See `hnf_internal.c` for the precise post-conditions.
pub fn ibz_xgcd_with_u_not_0(d: &mut Ibz, u: &mut Ibz, v: &mut Ibz, x: &Ibz, y: &Ibz) {
    if ibz_is_zero(x) & ibz_is_zero(y) {
        ibz_set(d, 1);
        ibz_set(u, 1);
        ibz_set(v, 0);
        return;
    }
    let mut q = Ibz::default();
    let mut r = Ibz::default();
    let x1 = x.clone();
    let mut y1 = y.clone();

    ibz_xgcd(d, u, v, &x1, &y1);

    // Ensure u != 0 (per GMP spec, u==0 implies y | x).
    if ibz_is_zero(u) {
        if !ibz_is_zero(&x1) {
            if ibz_is_zero(&y1) {
                ibz_set(&mut y1, 1);
            }
            ibz_div(&mut q, &mut r, &x1, &y1);
            debug_assert!(ibz_is_zero(&r));
            let t = v.clone();
            ibz_sub(v, &t, &q);
        }
        ibz_set(u, 1);
    }
    if !ibz_is_zero(&x1) {
        debug_assert!(ibz_cmp(d, ibz_const_zero()) > 0);
        ibz_mul(&mut r, &x1, &y1);
        let neg = ibz_cmp(&r, ibz_const_zero()) < 0;
        ibz_mul(&mut q, &x1, u);
        while ibz_cmp(&q, ibz_const_zero()) <= 0 {
            ibz_div(&mut q, &mut r, &y1, d);
            debug_assert!(ibz_is_zero(&r));
            if neg {
                let t = q.clone();
                ibz_neg(&mut q, &t);
            }
            let t = u.clone();
            ibz_add(u, &t, &q);
            ibz_div(&mut q, &mut r, &x1, d);
            debug_assert!(ibz_is_zero(&r));
            if neg {
                let t = q.clone();
                ibz_neg(&mut q, &t);
            }
            let t = v.clone();
            ibz_sub(v, &t, &q);

            ibz_mul(&mut q, &x1, u);
        }
    }

    #[cfg(debug_assertions)]
    {
        let mut sum = Ibz::default();
        let mut prod = Ibz::default();
        let mut res = false;
        res |= ibz_cmp(d, ibz_const_zero()) < 0;
        if ibz_is_zero(&x1) && ibz_is_zero(&y1) {
            res |= !(ibz_is_zero(v) && ibz_is_one(u) && ibz_is_one(d));
        } else {
            ibz_div(&mut sum, &mut prod, &x1, d);
            res |= !ibz_is_zero(&prod);
            ibz_div(&mut sum, &mut prod, &y1, d);
            res |= !ibz_is_zero(&prod);
            ibz_mul(&mut sum, &x1, u);
            ibz_mul(&mut prod, &y1, v);
            let t = sum.clone();
            ibz_add(&mut sum, &t, &prod);
            res |= ibz_cmp(&sum, d) != 0;
        }
        debug_assert!(!res);
    }
}

/// Check whether `mat` is in (lower-triangular column-style) HNF.
pub fn ibz_mat_4x4_is_hnf(mat: &IbzMat4x4) -> i32 {
    let mut res = true;
    let mut ind = 0;
    let zero = Ibz::default();
    for i in 0..4 {
        for j in 0..i {
            res = res && ibz_is_zero(&mat[i][j]);
        }
        let mut found = false;
        for j in i..4 {
            if found {
                res = res && ibz_cmp(&mat[i][j], &zero) >= 0;
                res = res && ibz_cmp(&mat[i][ind], &mat[i][j]) > 0;
            } else if !ibz_is_zero(&mat[i][j]) {
                found = true;
                ind = j;
                res = res && ibz_cmp(&mat[i][j], &zero) > 0;
            }
        }
    }
    // C also has a column-monotonicity loop here but it never updates
    // `linestart` (initialized to -1), so it is always vacuously true.
    // Ported faithfully as a no-op; this function is only a debug helper.
    let linestart: i32 = -1;
    for j in 0..4 {
        let mut i = 0;
        while i < 4 && ibz_is_zero(&mat[i][j]) {
            i += 1;
        }
        if i != 4 {
            res = res && (linestart < i as i32);
        }
    }
    res as i32
}

fn ibz_vec_4_linear_combination_mod(
    lc: &mut IbzVec4,
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
    modn: &Ibz,
) {
    let mut prod = Ibz::default();
    let m = modn.clone();
    let mut sums = ibz_vec_4_init();
    for i in 0..4 {
        ibz_mul(&mut sums[i], coeff_a, &vec_a[i]);
        ibz_mul(&mut prod, coeff_b, &vec_b[i]);
        let t = sums[i].clone();
        ibz_add(&mut sums[i], &t, &prod);
        let t = sums[i].clone();
        ibz_centered_mod(&mut sums[i], &t, &m);
    }
    ibz_vec_4_copy(lc, &sums);
}

fn ibz_vec_4_copy_mod(res: &mut IbzVec4, vec: &IbzVec4, modn: &Ibz) {
    let m = modn.clone();
    for i in 0..4 {
        ibz_centered_mod(&mut res[i], &vec[i], &m);
    }
}

fn ibz_vec_4_scalar_mul_mod(prod: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4, modn: &Ibz) {
    let m = modn.clone();
    let s = scalar.clone();
    for i in 0..4 {
        ibz_mul(&mut prod[i], &vec[i], &s);
        let t = prod[i].clone();
        ibz_mod(&mut prod[i], &t, &m);
    }
}

/// Modular HNF (Cohen, Algorithm 2.4.8). Input: `generator_number > 3`
/// generator vectors (rows of `generators` are `IbzVec4`); `mod` must be a
/// nonzero multiple of the lattice determinant. Output written column-major
/// into `hnf`.
pub fn ibz_mat_4xn_hnf_mod_core(
    hnf: &mut IbzMat4x4,
    generator_number: usize,
    generators: &[IbzVec4],
    modn: &Ibz,
) {
    debug_assert!(generator_number > 3);
    debug_assert!(generators.len() >= generator_number);
    debug_assert!(ibz_cmp(modn, ibz_const_zero()) > 0);
    let n = generator_number;
    let mut i: i32 = 3;
    let mut j: i32 = (n - 1) as i32;
    let mut k: i32 = (n - 1) as i32;

    let mut d = Ibz::default();
    let mut u = Ibz::default();
    let mut v = Ibz::default();
    let mut r = Ibz::default();
    let mut m = Ibz::default();
    let mut q = Ibz::default();
    let mut coeff_1 = Ibz::default();
    let mut coeff_2 = Ibz::default();
    let mut c = ibz_vec_4_init();
    let mut a: Vec<IbzVec4> = (0..n).map(|h| generators[h].clone()).collect();
    let mut w: [IbzVec4; 4] = core::array::from_fn(|_| ibz_vec_4_init());

    ibz_copy(&mut m, modn);

    while i != -1 {
        while j != 0 {
            j -= 1;
            if !ibz_is_zero(&a[j as usize][i as usize]) {
                let aki = a[k as usize][i as usize].clone();
                let aji = a[j as usize][i as usize].clone();
                ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &aki, &aji);
                let ak = a[k as usize].clone();
                let aj = a[j as usize].clone();
                ibz_vec_4_linear_combination(&mut c, &u, &ak, &v, &aj);
                ibz_div(&mut coeff_1, &mut r, &aki, &d);
                ibz_div(&mut coeff_2, &mut r, &aji, &d);
                let t = coeff_2.clone();
                ibz_neg(&mut coeff_2, &t);
                let aj = a[j as usize].clone();
                let ak = a[k as usize].clone();
                ibz_vec_4_linear_combination_mod(
                    &mut a[j as usize],
                    &coeff_1,
                    &aj,
                    &coeff_2,
                    &ak,
                    &m,
                );
                ibz_vec_4_copy_mod(&mut a[k as usize], &c, &m);
            }
        }
        let aki = a[k as usize][i as usize].clone();
        ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &aki, &m);
        let ak = a[k as usize].clone();
        ibz_vec_4_scalar_mul_mod(&mut w[i as usize], &u, &ak, &m);
        if ibz_is_zero(&w[i as usize][i as usize]) {
            ibz_copy(&mut w[i as usize][i as usize], &m);
        }
        for h in (i + 1)..4 {
            let whi = w[h as usize][i as usize].clone();
            let wii = w[i as usize][i as usize].clone();
            ibz_div_floor(&mut q, &mut r, &whi, &wii);
            let t = q.clone();
            ibz_neg(&mut q, &t);
            let wh = w[h as usize].clone();
            let wi = w[i as usize].clone();
            ibz_vec_4_linear_combination(&mut w[h as usize], ibz_const_one(), &wh, &q, &wi);
        }
        let mc = m.clone();
        ibz_div(&mut m, &mut r, &mc, &d);
        debug_assert!(ibz_is_zero(&r));
        if i != 0 {
            k -= 1;
            i -= 1;
            j = k;
            if ibz_is_zero(&a[k as usize][i as usize]) {
                ibz_copy(&mut a[k as usize][i as usize], &m);
            }
        } else {
            k -= 1;
            i -= 1;
            j = k;
        }
    }
    for jc in 0..4 {
        for ir in 0..4 {
            ibz_copy(&mut hnf[ir][jc], &w[jc][ir]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::intbig::ibz_from_i64;

    fn z(v: i64) -> Ibz {
        ibz_from_i64(v)
    }

    #[test]
    fn xgcd_basic() {
        let (mut g, mut u, mut v) = (Ibz::default(), Ibz::default(), Ibz::default());
        ibz_xgcd(&mut g, &mut u, &mut v, &z(12), &z(18));
        assert_eq!(g, z(6));
        assert_eq!(u.clone() * z(12) + v.clone() * z(18), g);
    }

    #[test]
    fn xgcd_with_u_not_0_properties() {
        for &(x, y) in &[(12, 18), (0, 5), (5, 0), (7, 14), (-6, 9), (4, -10)] {
            let (mut d, mut u, mut v) = (Ibz::default(), Ibz::default(), Ibz::default());
            ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &z(x), &z(y));
            assert!(d > z(0));
            assert!(u != z(0));
            assert_eq!(u.clone() * z(x) + v.clone() * z(y), d);
            if x != 0 {
                assert!(z(x) * &u > z(0));
            }
        }
        let (mut d, mut u, mut v) = (Ibz::default(), Ibz::default(), Ibz::default());
        ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &z(0), &z(0));
        assert_eq!(d, z(1));
        assert_eq!(u, z(1));
        assert_eq!(v, z(0));
    }

    #[test]
    fn centered_mod() {
        let m = z(7);
        for x in -20..20 {
            let mut r = Ibz::default();
            ibz_centered_mod(&mut r, &z(x), &m);
            assert!(r > z(-4) && r <= z(3), "x={x} r={r:?}");
            assert_eq!(ibz_mod_ui(&(z(x) - r), 7), 0);
        }
    }

    #[test]
    fn mod_not_zero() {
        let m = z(5);
        for x in -10..10 {
            let mut r = Ibz::default();
            ibz_mod_not_zero(&mut r, &z(x), &m);
            assert!(r > z(0) && r <= z(5));
            assert_eq!(ibz_mod_ui(&(z(x) - r), 5), 0);
        }
    }

    #[test]
    fn hnf_of_identity_scaled() {
        // HNF of 4 generators that are 5·I should be 5·I (mod 625).
        let gens: Vec<IbzVec4> = (0..4)
            .map(|i| core::array::from_fn(|j| z(if i == j { 5 } else { 0 })))
            .collect();
        let mut hnf = ibz_mat_4x4_init();
        ibz_mat_4xn_hnf_mod_core(&mut hnf, 4, &gens, &z(625));
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(hnf[i][j], if i == j { z(5) } else { z(0) }, "[{i}][{j}]");
            }
        }
        assert_eq!(ibz_mat_4x4_is_hnf(&hnf), 1);
    }

    #[test]
    fn hnf_recovers_lattice() {
        // 8 generators spanning Z^4 (the standard basis plus shifted copies).
        let gens: Vec<IbzVec4> = vec![
            [z(2), z(0), z(0), z(0)],
            [z(0), z(2), z(0), z(0)],
            [z(0), z(0), z(2), z(0)],
            [z(0), z(0), z(0), z(2)],
            [z(1), z(1), z(0), z(0)],
            [z(0), z(1), z(1), z(0)],
            [z(0), z(0), z(1), z(1)],
            [z(1), z(0), z(0), z(1)],
        ];
        let mut hnf = ibz_mat_4x4_init();
        // These span the index-2 sublattice {v : Σvᵢ even}.
        ibz_mat_4xn_hnf_mod_core(&mut hnf, 8, &gens, &z(16));
        assert_eq!(ibz_mat_4x4_is_hnf(&hnf), 1);
        let mut det = Ibz::default();
        ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &hnf);
        assert_eq!(det, z(2));
        // Adding e₀ generator yields all of Z^4.
        let mut gens2 = gens;
        gens2.push([z(1), z(0), z(0), z(0)]);
        ibz_mat_4xn_hnf_mod_core(&mut hnf, 9, &gens2, &z(16));
        ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &hnf);
        assert_eq!(det, z(1));
    }
}
