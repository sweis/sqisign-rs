// SPDX-License-Identifier: Apache-2.0
//! Norm-equation solving and special extremal orders, ported from `normeq.c`.

use super::algebra::*;
use super::dim4::*;
use super::ideal::*;
use super::intbig::*;
use super::integers::ibz_cornacchia_prime;
use super::lattice::*;
use super::types::*;

/// Standard maximal order O₀ ⊂ B_{p,∞}: basis (1, i, (i+j)/2, (1+ij)/2).
pub fn quat_lattice_o0_set(o0: &mut QuatLattice) {
    for i in 0..4 {
        for j in 0..4 {
            ibz_set(&mut o0.basis[i][j], 0);
        }
    }
    ibz_set(&mut o0.denom, 2);
    ibz_set(&mut o0.basis[0][0], 2);
    ibz_set(&mut o0.basis[1][1], 2);
    ibz_set(&mut o0.basis[2][2], 1);
    ibz_set(&mut o0.basis[1][2], 1);
    ibz_set(&mut o0.basis[3][3], 1);
    ibz_set(&mut o0.basis[0][3], 1);
}

pub fn quat_lattice_o0_set_extremal(o0: &mut QuatPExtremalMaximalOrder) {
    ibz_set(&mut o0.z.coord[1], 1);
    ibz_set(&mut o0.t.coord[2], 1);
    ibz_set(&mut o0.z.denom, 1);
    ibz_set(&mut o0.t.denom, 1);
    o0.q = 1;
    quat_lattice_o0_set(&mut o0.order);
}

pub fn quat_order_elem_create(
    elem: &mut QuatAlgElem,
    order: &QuatPExtremalMaximalOrder,
    coeffs: &IbzVec4,
    bpoo: &QuatAlg,
) {
    let mut quat_temp = QuatAlgElem::default();
    quat_alg_scalar(elem, &coeffs[0], ibz_const_one());
    quat_alg_scalar(&mut quat_temp, &coeffs[1], ibz_const_one());
    let qt = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.z, &qt, bpoo);
    let e = elem.clone();
    quat_alg_add(elem, &e, &quat_temp);
    quat_alg_scalar(&mut quat_temp, &coeffs[2], ibz_const_one());
    let qt = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &qt, bpoo);
    let e = elem.clone();
    quat_alg_add(elem, &e, &quat_temp);
    quat_alg_scalar(&mut quat_temp, &coeffs[3], ibz_const_one());
    let qt = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &order.t, &qt, bpoo);
    let qt = quat_temp.clone();
    quat_alg_mul(&mut quat_temp, &qt, &order.z, bpoo);
    let e = elem.clone();
    quat_alg_add(elem, &e, &quat_temp);
}

pub fn quat_represent_integer(
    gamma: &mut QuatAlgElem,
    n_gamma: &Ibz,
    non_diag: i32,
    params: &QuatRepresentIntegerParams,
) -> i32 {
    if ibz_is_even(n_gamma) != 0 {
        return 0;
    }
    let mut found = 0i32;
    let mut cornacchia_target = Ibz::new();
    let mut adjusted_n_gamma = Ibz::new();
    let mut q = Ibz::new();
    let mut bound = Ibz::new();
    let mut sq_bound = Ibz::new();
    let mut temp = Ibz::new();
    let mut coeffs = ibz_vec_4_init();

    if non_diag != 0 {
        debug_assert!(params.order.q % 4 == 1);
    }
    ibz_set(&mut q, params.order.q as i32);
    let standard_order = params.order.q == 1;

    if non_diag != 0 || standard_order {
        ibz_mul(&mut adjusted_n_gamma, n_gamma, ibz_const_two());
        let t = adjusted_n_gamma.clone();
        ibz_mul(&mut adjusted_n_gamma, &t, ibz_const_two());
    } else {
        ibz_copy(&mut adjusted_n_gamma, n_gamma);
    }

    ibz_div(
        &mut sq_bound,
        &mut bound,
        &adjusted_n_gamma,
        &params.algebra.p,
    );
    ibz_set(&mut temp, params.order.q as i32);
    let t = sq_bound.clone();
    ibz_sub(&mut sq_bound, &t, &temp);
    ibz_sqrt_floor(&mut bound, &sq_bound);

    let mut counter = Ibz::new();
    let t = temp.clone();
    ibz_mul(&mut temp, &t, &params.algebra.p);
    let t = temp.clone();
    ibz_mul(&mut temp, &t, &params.algebra.p);
    let t = temp.clone();
    ibz_sqrt_floor(&mut temp, &t);
    let tm = temp.clone();
    ibz_div(&mut counter, &mut temp, &adjusted_n_gamma, &tm);

    while found == 0 && ibz_cmp(&counter, ibz_const_zero()) != 0 {
        let c = counter.clone();
        ibz_sub(&mut counter, &c, ibz_const_one());

        ibz_rand_interval(&mut coeffs[2], ibz_const_one(), &bound);

        ibz_mul(&mut cornacchia_target, &coeffs[2], &coeffs[2]);
        ibz_mul(&mut temp, &cornacchia_target, &params.algebra.p);
        let t = adjusted_n_gamma.clone();
        let tt = temp.clone();
        ibz_sub(&mut temp, &t, &tt);
        ibz_mul(&mut sq_bound, &q, &params.algebra.p);
        let t = temp.clone();
        let sb = sq_bound.clone();
        ibz_div(&mut temp, &mut sq_bound, &t, &sb);
        let t = temp.clone();
        ibz_sqrt_floor(&mut temp, &t);

        if ibz_cmp(&temp, ibz_const_zero()) == 0 {
            continue;
        }
        ibz_rand_interval(&mut coeffs[3], ibz_const_one(), &temp);

        ibz_mul(&mut temp, &coeffs[3], &coeffs[3]);
        let t = temp.clone();
        ibz_mul(&mut temp, &q, &t);
        let ct = cornacchia_target.clone();
        ibz_add(&mut cornacchia_target, &ct, &temp);
        let ct = cornacchia_target.clone();
        ibz_mul(&mut cornacchia_target, &ct, &params.algebra.p);
        let ct = cornacchia_target.clone();
        ibz_sub(&mut cornacchia_target, &adjusted_n_gamma, &ct);
        debug_assert!(ibz_cmp(&cornacchia_target, ibz_const_zero()) > 0);

        if ibz_probab_prime(&cornacchia_target, params.primality_test_iterations) != 0 {
            let (mut c0, mut c1) = (Ibz::new(), Ibz::new());
            found = ibz_cornacchia_prime(&mut c0, &mut c1, &q, &cornacchia_target);
            coeffs[0] = c0;
            coeffs[1] = c1;
        } else {
            found = 0;
        }

        if found != 0 && non_diag != 0 && standard_order {
            if ibz_is_odd(&coeffs[0]) != ibz_is_odd(&coeffs[3]) {
                coeffs.swap(0, 1);
            }
            // C uses truncating `%`; for negative differences `(-2) % 4 == -2`, so this
            // rejects valid samples. Spec Alg 3.12 step 15 wants ≡2 (mod 4). We match C
            // for KAT determinism — see PORTING.md.
            let cond = ((ibz_get(&coeffs[0]) - ibz_get(&coeffs[3])) % 4 == 2)
                && ((ibz_get(&coeffs[1]) - ibz_get(&coeffs[2])) % 4 == 2);
            found = (found != 0 && cond) as i32;
        }
        if found != 0 {
            quat_order_elem_create(gamma, params.order, &coeffs, params.algebra);
            #[cfg(debug_assertions)]
            {
                let mut nn = Ibz::new();
                let mut nd = Ibz::new();
                quat_alg_norm(&mut nn, &mut nd, gamma, params.algebra);
                debug_assert!(ibz_is_one(&nd) != 0);
                debug_assert!(ibz_cmp(&nn, &adjusted_n_gamma) == 0);
                debug_assert!(quat_lattice_contains(None, &params.order.order, gamma) != 0);
            }
            quat_alg_make_primitive(&mut coeffs, &mut temp, gamma, &params.order.order);
            if non_diag != 0 || standard_order {
                found = (ibz_cmp(&temp, ibz_const_two()) == 0) as i32;
            } else {
                found = (ibz_cmp(&temp, ibz_const_one()) == 0) as i32;
            }
        }
    }

    if found != 0 {
        let c = coeffs.clone();
        ibz_mat_4x4_eval(&mut coeffs, &params.order.order.basis, &c);
        for i in 0..4 {
            ibz_copy(&mut gamma.coord[i], &coeffs[i]);
        }
        ibz_copy(&mut gamma.denom, &params.order.order.denom);
    }
    found
}

pub fn quat_sampling_random_ideal_o0_given_norm(
    lideal: &mut QuatLeftIdeal,
    norm: &Ibz,
    is_prime: bool,
    params: &QuatRepresentIntegerParams,
    prime_cofactor: Option<&Ibz>,
) -> i32 {
    let mut n_temp = Ibz::new();
    let mut norm_d = Ibz::new();
    let mut disc = Ibz::new();
    let mut gen = QuatAlgElem::default();
    let mut gen_rerand = QuatAlgElem::default();
    let mut found = 0i32;

    if is_prime {
        while found == 0 {
            ibz_set(&mut gen.coord[0], 0);
            ibz_sub(&mut n_temp, norm, ibz_const_one());
            for i in 1..4 {
                ibz_rand_interval(&mut gen.coord[i], ibz_const_zero(), &n_temp);
            }
            quat_alg_norm(&mut n_temp, &mut norm_d, &gen, params.algebra);
            debug_assert!(ibz_is_one(&norm_d) != 0);
            ibz_neg(&mut disc, &n_temp);
            let d = disc.clone();
            ibz_mod(&mut disc, &d, norm);
            found = ibz_sqrt_mod_p(&mut gen.coord[0], &disc, norm);
            found = (found != 0 && quat_alg_elem_is_zero(&gen) == 0) as i32;
        }
    } else {
        let pc = prime_cofactor.expect("prime_cofactor required when !is_prime");
        debug_assert!(ibz_is_zero(norm) == 0);
        ibz_mul(&mut n_temp, pc, norm);
        found = quat_represent_integer(&mut gen, &n_temp, 0, params);
        let _ = (found != 0 && quat_alg_elem_is_zero(&gen) == 0) as i32;
    }
    // C also discards `found` here before re-entering the rerandomization loop.
    found = 0;
    while found == 0 {
        for i in 0..4 {
            ibz_rand_interval(&mut gen_rerand.coord[i], ibz_const_one(), norm);
        }
        quat_alg_norm(&mut n_temp, &mut norm_d, &gen_rerand, params.algebra);
        debug_assert!(ibz_is_one(&norm_d) != 0);
        ibz_gcd(&mut disc, &n_temp, norm);
        found = (ibz_is_one(&disc) != 0 && quat_alg_elem_is_zero(&gen_rerand) == 0) as i32;
    }

    let g = gen.clone();
    quat_alg_mul(&mut gen, &g, &gen_rerand, params.algebra);
    quat_lideal_create(lideal, &gen, norm, &params.order.order, params.algebra);
    debug_assert!(ibz_cmp(norm, &lideal.norm) == 0);
    found
}

pub fn quat_change_to_o0_basis(vec: &mut IbzVec4, el: &QuatAlgElem) {
    let mut tmp = Ibz::new();
    ibz_copy(&mut vec[2], &el.coord[2]);
    let t = vec[2].clone();
    ibz_add(&mut vec[2], &t, &t);
    ibz_copy(&mut vec[3], &el.coord[3]);
    let t = vec[3].clone();
    ibz_add(&mut vec[3], &t, &t);
    ibz_sub(&mut vec[0], &el.coord[0], &el.coord[3]);
    ibz_sub(&mut vec[1], &el.coord[1], &el.coord[2]);

    debug_assert!(ibz_divides(&vec[0], &el.denom) != 0);
    debug_assert!(ibz_divides(&vec[1], &el.denom) != 0);
    debug_assert!(ibz_divides(&vec[2], &el.denom) != 0);
    debug_assert!(ibz_divides(&vec[3], &el.denom) != 0);

    for i in 0..4 {
        let t = vec[i].clone();
        ibz_div(&mut vec[i], &mut tmp, &t, &el.denom);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ctrdrbg;
    use rug::Integer;
    use std::sync::OnceLock;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }

    fn alg7() -> &'static QuatAlg {
        static A: OnceLock<QuatAlg> = OnceLock::new();
        A.get_or_init(|| QuatAlg::from_ui(7))
    }
    fn o0ext() -> &'static QuatPExtremalMaximalOrder {
        static O: OnceLock<QuatPExtremalMaximalOrder> = OnceLock::new();
        O.get_or_init(|| {
            let mut o = QuatPExtremalMaximalOrder::default();
            quat_lattice_o0_set_extremal(&mut o);
            o
        })
    }

    #[test]
    fn change_to_o0_basis_roundtrip() {
        // (1, i, (i+j)/2, (1+ij)/2) basis: take element (3+5i+7j+11ij)/2 — must lie in O₀.
        // Easier: use integer-coord element with denom 1: x = 2+4i+6j+8ij.
        let mut x = QuatAlgElem::default();
        quat_alg_elem_set(&mut x, 1, 2, 4, 6, 8);
        let mut v = ibz_vec_4_init();
        quat_change_to_o0_basis(&mut v, &x);
        // Reconstruct: O₀.basis · v / O₀.denom should equal x.coord.
        let mut o0 = QuatLattice::default();
        quat_lattice_o0_set(&mut o0);
        let mut back = ibz_vec_4_init();
        ibz_mat_4x4_eval(&mut back, &o0.basis, &v);
        for i in 0..4 {
            assert_eq!(back[i], Ibz::from(2) * &x.coord[i]); // denom=2
        }
    }

    #[test]
    fn represent_integer_finds_solution() {
        ctrdrbg::randombytes_init(&[7u8; 48]);
        let params = QuatRepresentIntegerParams {
            primality_test_iterations: 30,
            order: o0ext(),
            algebra: alg7(),
        };
        // Target: a moderately large odd integer >> p so the search space is non-trivial.
        let n = z(1_000_003);
        let mut gamma = QuatAlgElem::default();
        let ok = quat_represent_integer(&mut gamma, &n, 0, &params);
        assert_eq!(ok, 1);
        let mut nn = Ibz::new();
        let mut nd = Ibz::new();
        quat_alg_norm(&mut nn, &mut nd, &gamma, alg7());
        assert_eq!(nd, 1);
        assert_eq!(nn, n);
        assert_eq!(quat_lattice_contains(None, &o0ext().order, &gamma), 1);
    }
}
