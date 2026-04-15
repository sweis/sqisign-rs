//! Rejection sampling of short lattice vectors, ported from `lat_ball.c`.

use super::algebra::*;
use super::dim4::*;
use super::intbig::*;
use super::lattice::*;
use super::lll::*;
use super::types::*;

pub fn quat_lattice_bound_parallelogram(
    boxv: &mut IbzVec4,
    u: &mut IbzMat4x4,
    g: &IbzMat4x4,
    radius: &Ibz,
) -> i32 {
    let mut denom = Ibz::default();
    let mut rem = Ibz::default();
    let mut dual_g = ibz_mat_4x4_init();

    let inv_check = ibz_mat_4x4_inv_with_det_as_denom(Some(&mut dual_g), &mut denom, g);
    debug_assert!(inv_check != 0);
    ibz_mat_4x4_identity(u);
    quat_lll_core(&mut dual_g, u);

    let mut trivial = 1i32;
    for i in 0..4 {
        ibz_mul(&mut boxv[i], &dual_g[i][i], radius);
        let t = boxv[i].clone();
        ibz_div(&mut boxv[i], &mut rem, &t, &denom);
        let t = boxv[i].clone();
        ibz_sqrt_floor(&mut boxv[i], &t);
        trivial &= ibz_is_zero(&boxv[i]);
    }

    let uc0 = u.clone();
    let inv = ibz_mat_4x4_inv_with_det_as_denom(Some(u), &mut denom, &uc0);
    let uc = u.clone();
    ibz_mat_4x4_scalar_mul(u, &denom, &uc);
    #[cfg(debug_assertions)]
    {
        debug_assert!(inv != 0);
        let mut a = Ibz::default();
        ibz_abs(&mut a, &denom);
        debug_assert!(ibz_is_one(&a) != 0);
    }
    let _ = inv;
    (trivial == 0) as i32
}

pub fn quat_lattice_sample_from_ball(
    res: &mut QuatAlgElem,
    lattice: &QuatLattice,
    alg: &QuatAlg,
    radius: &Ibz,
) -> i32 {
    debug_assert!(ibz_cmp(radius, ibz_const_zero()) > 0);
    let mut boxv = ibz_vec_4_init();
    let mut u = ibz_mat_4x4_init();
    let mut g = ibz_mat_4x4_init();
    let mut x = ibz_vec_4_init();
    let mut rad = Ibz::default();
    let mut tmp = Ibz::default();

    quat_lattice_gram(&mut g, lattice, alg);
    ibz_mul(&mut rad, radius, &lattice.denom);
    let r = rad.clone();
    ibz_mul(&mut rad, &r, &lattice.denom);
    let r = rad.clone();
    ibz_mul(&mut rad, &r, ibz_const_two());

    let mut ok = quat_lattice_bound_parallelogram(&mut boxv, &mut u, &g, &rad);
    if ok == 0 {
        return 0;
    }

    #[cfg(debug_assertions)]
    let mut cnt = 0u32;
    loop {
        for i in 0..4 {
            if ibz_is_zero(&boxv[i]) != 0 {
                ibz_copy(&mut x[i], ibz_const_zero());
            } else {
                ibz_add(&mut tmp, &boxv[i], &boxv[i]);
                ok &= ibz_rand_interval(&mut x[i], ibz_const_zero(), &tmp);
                let t = x[i].clone();
                ibz_sub(&mut x[i], &t, &boxv[i]);
                if ok == 0 {
                    return 0;
                }
            }
        }
        let xc = x.clone();
        ibz_mat_4x4_eval_t(&mut x, &xc, &u);
        quat_qf_eval(&mut tmp, &g, &x);
        #[cfg(debug_assertions)]
        {
            cnt += 1;
            if cnt.is_multiple_of(100) {
                eprintln!("Lattice sampling rejected {} times", cnt - 1);
            }
        }
        if ibz_is_zero(&tmp) == 0 && ibz_cmp(&tmp, &rad) <= 0 {
            break;
        }
    }

    ibz_mat_4x4_eval(&mut res.coord, &lattice.basis, &x);
    ibz_copy(&mut res.denom, &lattice.denom);
    quat_alg_normalize(res);

    #[cfg(debug_assertions)]
    {
        let mut nn = Ibz::default();
        let mut nd = Ibz::default();
        quat_alg_norm(&mut nn, &mut nd, res, alg);
        let mut bound = Ibz::default();
        ibz_mul(&mut bound, &nd, radius);
        debug_assert!(ibz_cmp(&nn, &bound) <= 0);
    }
    ok
}

#[cfg(test)]
mod tests {
    use super::super::normeq::quat_lattice_o0_set;
    use super::*;
    use crate::common::ctrdrbg;

    #[test]
    fn sample_in_ball() {
        ctrdrbg::randombytes_init(&[0u8; 48]);
        let alg = QuatAlg::from_ui(103);
        let mut lat = QuatLattice::default();
        quat_lattice_o0_set(&mut lat);
        let radius = Ibz::from(50);
        let mut res = QuatAlgElem::default();
        assert_eq!(
            quat_lattice_sample_from_ball(&mut res, &lat, &alg, &radius),
            1
        );
        let mut nn = Ibz::default();
        let mut nd = Ibz::default();
        quat_alg_norm(&mut nn, &mut nd, &res, &alg);
        assert!(nn <= radius * &nd);
        assert_eq!(quat_lattice_contains(None, &lat, &res), 1);
    }
}
