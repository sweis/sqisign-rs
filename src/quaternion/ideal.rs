// SPDX-License-Identifier: Apache-2.0
//! Left ideals of maximal orders, ported from `ideal.c` and the
//! ideal-related functions in `lll/lll_applications.c`.

use super::algebra::*;
use super::dim4::*;
use super::intbig::*;
use super::lattice::*;
use super::lll::*;
use super::types::*;

/// Recompute and store `lideal.norm = √[order : lideal]`.
pub fn quat_lideal_norm(lideal: &mut QuatLeftIdeal) {
    quat_lattice_index(
        &mut lideal.norm,
        &lideal.lattice,
        lideal.parent_order.expect("parent_order"),
    );
    let n = lideal.norm.clone();
    let ok = ibz_sqrt(&mut lideal.norm, &n);
    debug_assert!(ok != 0);
}

pub fn quat_lideal_copy(copy: &mut QuatLeftIdeal, copied: &QuatLeftIdeal) {
    copy.parent_order = copied.parent_order;
    ibz_copy(&mut copy.norm, &copied.norm);
    ibz_copy(&mut copy.lattice.denom, &copied.lattice.denom);
    ibz_mat_4x4_copy(&mut copy.lattice.basis, &copied.lattice.basis);
}

pub fn quat_lideal_create_principal(
    lideal: &mut QuatLeftIdeal,
    x: &QuatAlgElem,
    order: &'static QuatLattice,
    alg: &QuatAlg,
) {
    debug_assert!(quat_order_is_maximal(order, alg) != 0);
    debug_assert!(quat_lattice_contains(None, order, x) != 0);
    let mut norm_n = Ibz::default();
    let mut norm_d = Ibz::default();
    quat_lattice_alg_elem_mul(&mut lideal.lattice, order, x, alg);
    let l = lideal.lattice.clone();
    quat_lattice_reduce_denom(&mut lideal.lattice, &l);
    quat_alg_norm(&mut norm_n, &mut norm_d, x, alg);
    debug_assert!(ibz_is_one(&norm_d) != 0);
    ibz_copy(&mut lideal.norm, &norm_n);
    lideal.parent_order = Some(order);
}

pub fn quat_lideal_create(
    lideal: &mut QuatLeftIdeal,
    x: &QuatAlgElem,
    n: &Ibz,
    order: &'static QuatLattice,
    alg: &QuatAlg,
) {
    debug_assert!(quat_order_is_maximal(order, alg) != 0);
    debug_assert!(quat_alg_elem_is_zero(x) == 0);
    let mut on = QuatLattice::default();
    quat_lideal_create_principal(lideal, x, order, alg);
    ibz_mat_4x4_scalar_mul(&mut on.basis, n, &order.basis);
    ibz_copy(&mut on.denom, &order.denom);
    let l = lideal.lattice.clone();
    quat_lattice_add(&mut lideal.lattice, &l, &on);
    lideal.parent_order = Some(order);
    quat_lideal_norm(lideal);
}

pub fn quat_lideal_generator(gen: &mut QuatAlgElem, lideal: &QuatLeftIdeal, alg: &QuatAlg) -> i32 {
    let mut norm_int = Ibz::default();
    let mut norm_denom = Ibz::default();
    let mut gcd = Ibz::default();
    let mut r = Ibz::default();
    let mut q = Ibz::default();
    let mut vec = ibz_vec_4_init();
    let mut int_norm = 0i32;
    loop {
        int_norm += 1;
        for a in -int_norm..=int_norm {
            let ba = int_norm - a.abs();
            for b in -ba..=ba {
                let ca = ba - b.abs();
                for c in -ca..=ca {
                    let d = ca - c.abs();
                    ibz_vec_4_set(&mut vec, a, b, c, d);
                    ibz_vec_4_content(&mut gcd, &vec);
                    if ibz_is_one(&gcd) != 0 {
                        ibz_mat_4x4_eval(&mut gen.coord, &lideal.lattice.basis, &vec);
                        ibz_copy(&mut gen.denom, &lideal.lattice.denom);
                        quat_alg_norm(&mut norm_int, &mut norm_denom, gen, alg);
                        debug_assert!(ibz_is_one(&norm_denom) != 0);
                        ibz_div(&mut q, &mut r, &norm_int, &lideal.norm);
                        debug_assert!(ibz_is_zero(&r) != 0);
                        ibz_gcd(&mut gcd, &lideal.norm, &q);
                        if ibz_cmp(&gcd, ibz_const_one()) == 0 {
                            return 1;
                        }
                    }
                }
            }
        }
    }
}

pub fn quat_lideal_mul(
    product: &mut QuatLeftIdeal,
    lideal: &QuatLeftIdeal,
    alpha: &QuatAlgElem,
    alg: &QuatAlg,
) {
    debug_assert!(quat_order_is_maximal(lideal.parent_order.unwrap(), alg) != 0);
    let mut norm = Ibz::default();
    let mut norm_d = Ibz::default();
    quat_lattice_alg_elem_mul(&mut product.lattice, &lideal.lattice, alpha, alg);
    product.parent_order = lideal.parent_order;
    quat_alg_norm(&mut norm, &mut norm_d, alpha, alg);
    ibz_mul(&mut product.norm, &lideal.norm, &norm);
    debug_assert!(ibz_divides(&product.norm, &norm_d) != 0);
    let pn = product.norm.clone();
    ibz_div(&mut product.norm, &mut norm, &pn, &norm_d);
}

pub fn quat_lideal_add(
    sum: &mut QuatLeftIdeal,
    i1: &QuatLeftIdeal,
    i2: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    debug_assert!(core::ptr::eq(
        i1.parent_order.unwrap(),
        i2.parent_order.unwrap()
    ));
    debug_assert!(quat_order_is_maximal(i2.parent_order.unwrap(), alg) != 0);
    quat_lattice_add(&mut sum.lattice, &i1.lattice, &i2.lattice);
    sum.parent_order = i1.parent_order;
    quat_lideal_norm(sum);
}

pub fn quat_lideal_inter(
    inter: &mut QuatLeftIdeal,
    i1: &QuatLeftIdeal,
    i2: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    debug_assert!(core::ptr::eq(
        i1.parent_order.unwrap(),
        i2.parent_order.unwrap()
    ));
    debug_assert!(quat_order_is_maximal(i2.parent_order.unwrap(), alg) != 0);
    quat_lattice_intersect(&mut inter.lattice, &i1.lattice, &i2.lattice);
    inter.parent_order = i1.parent_order;
    quat_lideal_norm(inter);
}

pub fn quat_lideal_equals(i1: &QuatLeftIdeal, i2: &QuatLeftIdeal, alg: &QuatAlg) -> i32 {
    debug_assert!(quat_order_is_maximal(i2.parent_order.unwrap(), alg) != 0);
    debug_assert!(quat_order_is_maximal(i1.parent_order.unwrap(), alg) != 0);
    (core::ptr::eq(i1.parent_order.unwrap(), i2.parent_order.unwrap()) as i32)
        & ((ibz_cmp(&i1.norm, &i2.norm) == 0) as i32)
        & quat_lattice_equal(&i1.lattice, &i2.lattice)
}

pub fn quat_lideal_inverse_lattice_without_hnf(
    inv: &mut QuatLattice,
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    debug_assert!(quat_order_is_maximal(lideal.parent_order.unwrap(), alg) != 0);
    quat_lattice_conjugate_without_hnf(inv, &lideal.lattice);
    let d = inv.denom.clone();
    ibz_mul(&mut inv.denom, &d, &lideal.norm);
}

pub fn quat_lideal_right_transporter(
    trans: &mut QuatLattice,
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    debug_assert!(core::ptr::eq(
        lideal1.parent_order.unwrap(),
        lideal2.parent_order.unwrap()
    ));
    let mut inv = QuatLattice::default();
    quat_lideal_inverse_lattice_without_hnf(&mut inv, lideal1, alg);
    quat_lattice_mul(trans, &inv, &lideal2.lattice, alg);
}

pub fn quat_lideal_right_order(order: &mut QuatLattice, lideal: &QuatLeftIdeal, alg: &QuatAlg) {
    quat_lideal_right_transporter(order, lideal, lideal, alg);
}

pub fn quat_lideal_class_gram(g: &mut IbzMat4x4, lideal: &QuatLeftIdeal, alg: &QuatAlg) {
    quat_lattice_gram(g, &lideal.lattice, alg);
    let mut divisor = Ibz::default();
    let mut rmd = Ibz::default();
    ibz_mul(&mut divisor, &lideal.lattice.denom, &lideal.lattice.denom);
    let d = divisor.clone();
    ibz_mul(&mut divisor, &d, &lideal.norm);
    for i in 0..4 {
        for j in 0..=i {
            let t = g[i][j].clone();
            ibz_div(&mut g[i][j], &mut rmd, &t, &divisor);
            debug_assert!(ibz_is_zero(&rmd) != 0);
        }
    }
    for i in 0..4 {
        for j in 0..i {
            let t = g[i][j].clone();
            ibz_copy(&mut g[j][i], &t);
        }
    }
}

pub fn quat_lideal_conjugate_without_hnf(
    conj: &mut QuatLeftIdeal,
    new_parent_order: &'static QuatLattice,
    lideal: &QuatLeftIdeal,
    _alg: &QuatAlg,
) {
    // C computes the right order into `new_parent_order`; in Rust we accept it
    // precomputed since it must outlive `conj` (`&'static`).
    quat_lattice_conjugate_without_hnf(&mut conj.lattice, &lideal.lattice);
    conj.parent_order = Some(new_parent_order);
    ibz_copy(&mut conj.norm, &lideal.norm);
}

/// Variant that mirrors C's behaviour exactly: computes the right order
/// into `right_order_out` and sets `conj.parent_order = None` (since
/// `parent_order` requires `'static` in Rust). Callers that need the
/// resulting parent order pass `right_order_out` explicitly to the
/// `*_dyn` family below.
pub fn quat_lideal_conjugate_without_hnf_dyn(
    conj: &mut QuatLeftIdeal,
    right_order_out: &mut QuatLattice,
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    quat_lideal_right_order(right_order_out, lideal, alg);
    quat_lattice_conjugate_without_hnf(&mut conj.lattice, &lideal.lattice);
    conj.parent_order = None;
    ibz_copy(&mut conj.norm, &lideal.norm);
}

pub fn quat_order_discriminant(disc: &mut Ibz, order: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut det = Ibz::default();
    let mut sqr = Ibz::default();
    let mut div = Ibz::default();
    let mut transposed = ibz_mat_4x4_init();
    let mut norm = ibz_mat_4x4_init();
    let mut prod = ibz_mat_4x4_init();
    ibz_mat_4x4_transpose(&mut transposed, &order.basis);
    ibz_mat_4x4_identity(&mut norm);
    ibz_copy(&mut norm[2][2], &alg.p);
    ibz_copy(&mut norm[3][3], &alg.p);
    let n = norm.clone();
    ibz_mat_4x4_scalar_mul(&mut norm, ibz_const_two(), &n);
    ibz_mat_4x4_mul(&mut prod, &transposed, &norm);
    let p = prod.clone();
    ibz_mat_4x4_mul(&mut prod, &p, &order.basis);
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &prod);
    ibz_mul(&mut div, &order.denom, &order.denom);
    let d = div.clone();
    ibz_mul(&mut div, &d, &d);
    let d = div.clone();
    ibz_mul(&mut div, &d, &d);
    let dv = div.clone();
    ibz_div(&mut sqr, &mut div, &det, &dv);
    let mut ok = ibz_is_zero(&div);
    ok &= ibz_sqrt(disc, &sqr);
    ok
}

pub fn quat_order_is_maximal(order: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut disc = Ibz::default();
    quat_order_discriminant(&mut disc, order, alg);
    (ibz_cmp(&disc, &alg.p) == 0) as i32
}

// ---------------------------------------------------------------------------
// LLL-based reductions (lll_applications.c).

pub fn quat_lideal_reduce_basis(
    reduced: &mut IbzMat4x4,
    gram: &mut IbzMat4x4,
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    #[cfg(debug_assertions)]
    if let Some(po) = lideal.parent_order {
        debug_assert!(quat_order_is_maximal(po, alg) != 0);
    }
    let mut gram_corrector = Ibz::default();
    ibz_mul(
        &mut gram_corrector,
        &lideal.lattice.denom,
        &lideal.lattice.denom,
    );
    quat_lideal_class_gram(gram, lideal, alg);
    ibz_mat_4x4_copy(reduced, &lideal.lattice.basis);
    quat_lll_core(gram, reduced);
    let g = gram.clone();
    ibz_mat_4x4_scalar_mul(gram, &gram_corrector, &g);
    for i in 0..4 {
        let t = gram[i][i].clone();
        ibz_div_2exp(&mut gram[i][i], &t, 1);
        for j in (i + 1)..4 {
            ibz_set(&mut gram[i][j], 0);
        }
    }
}

pub fn quat_lideal_lideal_mul_reduced(
    prod: &mut QuatLeftIdeal,
    gram: &mut IbzMat4x4,
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    let mut red = ibz_mat_4x4_init();
    quat_lattice_mul(&mut prod.lattice, &lideal1.lattice, &lideal2.lattice, alg);
    prod.parent_order = lideal1.parent_order;
    quat_lideal_norm(prod);
    quat_lideal_reduce_basis(&mut red, gram, prod, alg);
    ibz_mat_4x4_copy(&mut prod.lattice.basis, &red);
}

/// Variant where `lideal1.parent_order` is not `'static` and is supplied
/// explicitly. Computes the norm via `parent_order` and leaves
/// `prod.parent_order = None`.
pub fn quat_lideal_lideal_mul_reduced_dyn(
    prod: &mut QuatLeftIdeal,
    gram: &mut IbzMat4x4,
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    parent_order: &QuatLattice,
    alg: &QuatAlg,
) {
    let mut red = ibz_mat_4x4_init();
    quat_lattice_mul(&mut prod.lattice, &lideal1.lattice, &lideal2.lattice, alg);
    prod.parent_order = None;
    quat_lattice_index(&mut prod.norm, &prod.lattice, parent_order);
    let n = prod.norm.clone();
    let ok = ibz_sqrt(&mut prod.norm, &n);
    debug_assert!(ok != 0);
    quat_lideal_reduce_basis(&mut red, gram, prod, alg);
    ibz_mat_4x4_copy(&mut prod.lattice.basis, &red);
}

pub fn quat_lideal_prime_norm_reduced_equivalent(
    lideal: &mut QuatLeftIdeal,
    alg: &QuatAlg,
    primality_num_iter: i32,
    equiv_bound_coeff: i32,
) -> i32 {
    let mut gram = ibz_mat_4x4_init();
    let mut red = ibz_mat_4x4_init();
    let li = lideal.clone();
    quat_lideal_reduce_basis(&mut red, &mut gram, &li, alg);

    let mut new_alpha = QuatAlgElem::default();
    let mut tmp = Ibz::default();
    let mut remainder = Ibz::default();
    let mut adjusted_norm = Ibz::default();
    ibz_mul(
        &mut adjusted_norm,
        &lideal.lattice.denom,
        &lideal.lattice.denom,
    );

    debug_assert!(equiv_bound_coeff < (1 << 20));
    // (2B+1)^4 with B<2^20 can overflow i64; the actual SQIsign parameter is small
    // (QUAT_equiv_bound_coeff = 4 at lvl1) but use saturating arithmetic for safety.
    let base = (2 * equiv_bound_coeff as i64 + 1).saturating_pow(2);
    let equiv_num_iter = base.saturating_mul(base);

    let mut found = 0;
    let mut ctr: i64 = 0;
    while found == 0 && ctr < equiv_num_iter {
        ctr += 1;
        for i in 0..4 {
            ibz_rand_interval_minm_m(&mut new_alpha.coord[i], equiv_bound_coeff);
        }
        quat_qf_eval(&mut tmp, &gram, &new_alpha.coord);
        let t = tmp.clone();
        ibz_div(&mut tmp, &mut remainder, &t, &adjusted_norm);
        debug_assert!(ibz_is_zero(&remainder) != 0);
        if ibz_probab_prime(&tmp, primality_num_iter) != 0 {
            let coord = new_alpha.coord.clone();
            ibz_mat_4x4_eval(&mut new_alpha.coord, &red, &coord);
            ibz_copy(&mut new_alpha.denom, &lideal.lattice.denom);
            debug_assert!(quat_lattice_contains(None, &lideal.lattice, &new_alpha) != 0);
            let na = new_alpha.clone();
            quat_alg_conj(&mut new_alpha, &na);
            let d = new_alpha.denom.clone();
            ibz_mul(&mut new_alpha.denom, &d, &lideal.norm);
            let li = lideal.clone();
            quat_lideal_mul(lideal, &li, &new_alpha, alg);
            debug_assert!(ibz_probab_prime(&lideal.norm, primality_num_iter) != 0);
            found = 1;
        }
    }
    debug_assert!(found != 0);
    found
}

#[cfg(test)]
mod tests {
    use super::super::normeq::quat_lattice_o0_set;
    use super::*;
    use crate::quaternion::intbig::ibz_from_i64;
    use std::sync::OnceLock;

    fn z(v: i64) -> Ibz {
        ibz_from_i64(v)
    }

    fn o0() -> &'static QuatLattice {
        static O0: OnceLock<QuatLattice> = OnceLock::new();
        O0.get_or_init(|| {
            let mut o = QuatLattice::default();
            quat_lattice_o0_set(&mut o);
            o
        })
    }

    #[test]
    fn o0_is_maximal() {
        let alg = QuatAlg::from_ui(7);
        assert_eq!(quat_order_is_maximal(o0(), &alg), 1);
        let alg2 = QuatAlg::from_ui(103);
        assert_eq!(quat_order_is_maximal(o0(), &alg2), 1);
    }

    #[test]
    fn create_and_generator_roundtrip() {
        let alg = QuatAlg::from_ui(7);
        // x = 1 + i + j (in O₀; note (i+j)/2 ∈ O₀ but use integer coords).
        let mut x = QuatAlgElem::default();
        quat_alg_elem_set(&mut x, 1, 1, 1, 1, 0);
        assert_eq!(quat_lattice_contains(None, o0(), &x), 1);
        let n = z(3);
        let mut lideal = QuatLeftIdeal::default();
        quat_lideal_create(&mut lideal, &x, &n, o0(), &alg);
        // Nrd(x) = 1+1+7 = 9 = 3², gcd with N=3 → norm should divide 3.
        assert!(lideal.norm == 3 || lideal.norm == 1);
        // Find a generator and recreate.
        let mut gen = QuatAlgElem::default();
        assert_eq!(quat_lideal_generator(&mut gen, &lideal, &alg), 1);
        let mut lideal2 = QuatLeftIdeal::default();
        quat_lideal_create(&mut lideal2, &gen, &lideal.norm, o0(), &alg);
        assert_eq!(quat_lideal_equals(&lideal, &lideal2, &alg), 1);
    }

    #[test]
    fn right_order_is_maximal() {
        let alg = QuatAlg::from_ui(103);
        let mut x = QuatAlgElem::default();
        quat_alg_elem_set(&mut x, 1, 2, 0, 1, 0);
        let n = z(5);
        let mut lideal = QuatLeftIdeal::default();
        quat_lideal_create(&mut lideal, &x, &n, o0(), &alg);
        let mut ro = QuatLattice::default();
        quat_lideal_right_order(&mut ro, &lideal, &alg);
        assert_eq!(quat_order_is_maximal(&ro, &alg), 1);
    }
}
