//! Dimension-2 ideal-to-isogeny conversion (CLAPOTIS).
//! Ported from `src/id2iso/ref/lvlx/dim2id2iso.c`.

use crate::ec::biextension::weil;
use crate::ec::*;
use crate::gf::Fp2;
use crate::hd::*;
use crate::mp::Digit;
use crate::precomp::{
    connecting_ideals, curves_with_endomorphisms, extremal_orders, quatalg_pinfty, FINDUV_BOX_SIZE,
    FINDUV_CUBE_SIZE, NUM_ALTERNATE_EXTREMAL_ORDERS, NWORDS_ORDER, QUAT_REPRES_BOUND_INPUT,
    TORSION_EVEN_POWER,
};
use crate::quaternion::*;

use super::{endomorphism_application_even_basis, torsion_plus_2power};

/// `ALTERNATE_CONNECTING_IDEALS[i]` is `CONNECTING_IDEALS[i+1]`.
#[inline]
fn alternate_connecting_ideal(i: usize) -> &'static QuatLeftIdeal {
    &connecting_ideals()[i + 1]
}

// ---------------------------------------------------------------------------
// fixed_degree_isogeny_and_eval

fn fixed_degree_isogeny_impl(
    lideal: &mut QuatLeftIdeal,
    u: &Ibz,
    small: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    index_alternate_order: usize,
) -> i32 {
    let mut two_pow = Ibz::default();
    let mut tmp = Ibz::default();
    let mut theta = QuatAlgElem::default();

    let mut e0 = curves_with_endomorphisms()[index_alternate_order].curve;
    e0.normalize_a24();

    let u_bitsize = ibz_bitsize(u);

    let length: u32 = if small {
        let l = ibz_bitsize(&quatalg_pinfty().p) as u32 + QUAT_REPRES_BOUND_INPUT as u32
            - u_bitsize as u32;
        debug_assert!(u_bitsize < l as i32);
        debug_assert!(l < (TORSION_EVEN_POWER as u32) - HD_EXTRA_TORSION);
        l
    } else {
        (TORSION_EVEN_POWER as u32) - HD_EXTRA_TORSION
    };
    debug_assert!(length > 0);

    ibz_pow(&mut two_pow, ibz_const_two(), length);
    tmp.clone_from(u);
    debug_assert!(ibz_cmp(&two_pow, &tmp) > 0);
    debug_assert!(!ibz_is_even(&tmp));

    // theta has norm u·(2^length − u).
    let t = tmp.clone();
    ibz_sub(&mut tmp, &two_pow, &t);
    tmp *= u;
    debug_assert!(!ibz_is_even(&tmp));

    let ri_params = QuatRepresentIntegerParams {
        primality_test_iterations: super::quat_represent_integer_params().primality_test_iterations,
        order: &extremal_orders()[index_alternate_order],
        algebra: quatalg_pinfty(),
    };

    #[cfg(debug_assertions)]
    {
        debug_assert!(quat_lattice_contains(None, &ri_params.order.order, &ri_params.order.z) != 0);
        debug_assert!(quat_lattice_contains(None, &ri_params.order.order, &ri_params.order.t) != 0);
    }

    let ret = quat_represent_integer(&mut theta, &tmp, 1, &ri_params);
    debug_assert!(!ibz_is_even(&tmp));
    if ret == 0 {
        #[cfg(debug_assertions)]
        eprintln!(
            "represent integer failed for the alternate order number {} and for a target of size \
             {} for a u of size {} with length = {}",
            index_alternate_order,
            ibz_bitsize(&tmp),
            ibz_bitsize(u),
            length
        );
        return 0;
    }
    quat_lideal_create(
        lideal,
        &theta,
        u,
        &extremal_orders()[index_alternate_order].order,
        quatalg_pinfty(),
    );

    #[cfg(debug_assertions)]
    {
        let mut test_norm = Ibz::default();
        let mut test_denom = Ibz::default();
        quat_alg_norm(&mut test_norm, &mut test_denom, &theta, quatalg_pinfty());
        debug_assert!(ibz_is_one(&test_denom));
        debug_assert!(ibz_cmp(&test_norm, &tmp) == 0);
        debug_assert!(!ibz_is_even(&tmp));
        debug_assert!(
            quat_lattice_contains(
                None,
                &extremal_orders()[index_alternate_order].order,
                &theta
            ) != 0
        );
    }

    let mut b0_two = curves_with_endomorphisms()[index_alternate_order].basis_even;
    debug_assert!(test_basis_order_twof(
        &b0_two,
        &e0,
        TORSION_EVEN_POWER as i32
    ));
    let s = b0_two;
    ec_dbl_iter_basis(
        &mut b0_two,
        TORSION_EVEN_POWER as i32 - length as i32 - HD_EXTRA_TORSION as i32,
        &s,
        &mut e0,
    );
    debug_assert!(test_basis_order_twof(
        &b0_two,
        &e0,
        (length + HD_EXTRA_TORSION) as i32
    ));

    // theta ← theta · u⁻¹ mod 2^{length+2}
    let tp = two_pow.clone();
    ibz_mul(&mut two_pow, &tp, ibz_const_two());
    let tp = two_pow.clone();
    ibz_mul(&mut two_pow, &tp, ibz_const_two());
    tmp.clone_from(u);
    let t = tmp.clone();
    ibz_invmod(&mut tmp, &t, &two_pow);
    debug_assert!(!ibz_is_even(&tmp));
    for i in 0..4 {
        theta.coord[i] *= &tmp;
    }

    let mut b0_two_theta = b0_two;
    endomorphism_application_even_basis(
        &mut b0_two_theta,
        index_alternate_order,
        &e0,
        &theta,
        (length + HD_EXTRA_TORSION) as i32,
    );
    debug_assert!(test_basis_order_twof(
        &b0_two_theta,
        &e0,
        (length + HD_EXTRA_TORSION) as i32
    ));

    let mut e00 = ThetaCoupleCurve { e1: e0, e2: e0 };
    let dim_two_ker = ThetaKernelCouplePoints::from_bases(&b0_two, &b0_two_theta);

    if !theta_chain_compute_and_eval(length, &mut e00, &dim_two_ker, true, e34, p12) {
        return 0;
    }
    debug_assert!(length > 0);
    length as i32
}

/// Compute a 2D isogeny of degree `2^length` embedding a degree-`u`
/// 1D isogeny from the alternate starting curve `index_alternate_order`,
/// and evaluate it at `p12`. Returns the chain length or 0 on failure.
pub fn fixed_degree_isogeny_and_eval(
    lideal: &mut QuatLeftIdeal,
    u: &Ibz,
    small: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    index_alternate_order: usize,
) -> i32 {
    fixed_degree_isogeny_impl(lideal, u, small, e34, p12, index_alternate_order)
}

// ---------------------------------------------------------------------------
// LLL post-processing and hypercube enumeration.

fn post_lll_basis_treatment(
    gram: &mut IbzMat4x4,
    reduced: &mut IbzMat4x4,
    _norm: &Ibz,
    is_special_order: bool,
) {
    if !is_special_order {
        return;
    }
    if ibz_cmp(&gram[0][0], &gram[2][2]) == 0 {
        for i in 0..4 {
            let (a, b) = (reduced[i][1].clone(), reduced[i][2].clone());
            reduced[i][1] = b;
            reduced[i][2] = a;
        }
        swap_gram(gram, 0, 2, 0, 1);
        swap_gram(gram, 2, 0, 1, 0);
        swap_gram(gram, 3, 2, 3, 1);
        swap_gram(gram, 2, 3, 1, 3);
        swap_gram(gram, 2, 2, 1, 1);
    } else if ibz_cmp(&gram[0][0], &gram[3][3]) == 0 {
        for i in 0..4 {
            let (a, b) = (reduced[i][1].clone(), reduced[i][3].clone());
            reduced[i][1] = b;
            reduced[i][3] = a;
        }
        swap_gram(gram, 0, 3, 0, 1);
        swap_gram(gram, 3, 0, 1, 0);
        swap_gram(gram, 2, 3, 2, 1);
        swap_gram(gram, 3, 2, 1, 2);
        swap_gram(gram, 3, 3, 1, 1);
    } else if ibz_cmp(&gram[1][1], &gram[3][3]) == 0 {
        for i in 0..4 {
            let (a, b) = (reduced[i][1].clone(), reduced[i][2].clone());
            reduced[i][1] = b;
            reduced[i][2] = a;
        }
        swap_gram(gram, 0, 2, 0, 1);
        swap_gram(gram, 2, 0, 1, 0);
        swap_gram(gram, 3, 2, 3, 1);
        swap_gram(gram, 2, 3, 1, 3);
        swap_gram(gram, 2, 2, 1, 1);
    }

    if ibz_cmp(&reduced[0][0], &reduced[1][1]) != 0 {
        for i in 0..4 {
            let r = reduced[i][1].clone();
            ibz_neg(&mut reduced[i][1], &r);
            let g = gram[i][1].clone();
            ibz_neg(&mut gram[i][1], &g);
            let g = gram[1][i].clone();
            ibz_neg(&mut gram[1][i], &g);
        }
    }
    if ibz_cmp(&reduced[0][2], &reduced[1][3]) != 0 {
        for i in 0..4 {
            let r = reduced[i][3].clone();
            ibz_neg(&mut reduced[i][3], &r);
            let g = gram[i][3].clone();
            ibz_neg(&mut gram[i][3], &g);
            let g = gram[3][i].clone();
            ibz_neg(&mut gram[3][i], &g);
        }
    }
}

#[inline]
fn swap_gram(m: &mut IbzMat4x4, r1: usize, c1: usize, r2: usize, c2: usize) {
    let a = m[r1][c1].clone();
    let b = m[r2][c2].clone();
    m[r1][c1] = b;
    m[r2][c2] = a;
}

/// Enumerate vectors of ∞-norm ≤ m in the lattice with given Gram matrix,
/// recording those with odd reduced norm. Returns `count - 1` (matches C).
fn enumerate_hypercube(
    vecs: &mut [IbzVec4],
    norms: &mut [Ibz],
    m: i32,
    gram: &IbzMat4x4,
    adjusted_norm: &Ibz,
) -> i32 {
    let mut remain = Ibz::default();
    let mut norm = Ibz::default();
    let mut point = ibz_vec_4_init();

    debug_assert!(m > 0);
    let mut count: usize = 0;
    let dim = 2 * m + 1;
    let dim2 = dim * dim;
    let dim3 = dim2 * dim;

    let need_remove_symmetry =
        ibz_cmp(&gram[0][0], &gram[1][1]) == 0 && ibz_cmp(&gram[3][3], &gram[2][2]) == 0;

    for x in -m..=0 {
        for y in -m..=m {
            if x == 0 && y > 0 {
                break;
            }
            for z in -m..=m {
                if x == 0 && y == 0 && z > 0 {
                    break;
                }
                for w in -m..=m {
                    if x == 0 && y == 0 && z == 0 && w >= 0 {
                        break;
                    }
                    if (x | y | z | w) & 1 == 0 {
                        continue;
                    }
                    if x % 3 == 0 && y % 3 == 0 && z % 3 == 0 && w % 3 == 0 {
                        continue;
                    }

                    let check1 = (m + w) + dim * (m + z) + dim2 * (m + y) + dim3 * (m + x);
                    let check2 = (m - z) + dim * (m + w) + dim2 * (m - x) + dim3 * (m + y);
                    let check3 = (m + z) + dim * (m - w) + dim2 * (m + x) + dim3 * (m - y);

                    if !need_remove_symmetry || (check1 <= check2 && check1 <= check3) {
                        ibz_set(&mut point[0], x);
                        ibz_set(&mut point[1], y);
                        ibz_set(&mut point[2], z);
                        ibz_set(&mut point[3], w);
                        quat_qf_eval(&mut norm, gram, &point);
                        let n = norm.clone();
                        ibz_div(&mut norm, &mut remain, &n, adjusted_norm);
                        debug_assert!(ibz_is_zero(&remain));

                        if ibz_mod_ui(&norm, 2) == 1 {
                            ibz_set(&mut vecs[count][0], x);
                            ibz_set(&mut vecs[count][1], y);
                            ibz_set(&mut vecs[count][2], z);
                            ibz_set(&mut vecs[count][3], w);
                            norms[count].clone_from(&norm);
                            count += 1;
                        }
                    }
                }
            }
        }
    }
    count as i32 - 1
}

#[allow(clippy::too_many_arguments)]
fn find_uv_from_lists(
    au: &mut Ibz,
    bu: &mut Ibz,
    av: &mut Ibz,
    bv: &mut Ibz,
    u: &mut Ibz,
    v: &mut Ibz,
    index_sol1: &mut usize,
    index_sol2: &mut usize,
    target: &Ibz,
    small_norms1: &[Ibz],
    small_norms2: &[Ibz],
    quotients: &[Ibz],
    index1: usize,
    index2: usize,
    is_diagonal: bool,
    number_sum_square: i32,
) -> i32 {
    let mut n = Ibz::default();
    let mut remain = Ibz::default();
    let mut adjusted_norm = Ibz::default();
    let mut found = 0;
    n.clone_from(target);

    'outer: for i1 in 0..index1 {
        ibz_mod(&mut adjusted_norm, &n, &small_norms1[i1]);
        let starting_index2 = if is_diagonal { i1 } else { 0 };
        for i2 in starting_index2..index2 {
            if ibz_invmod(&mut remain, &small_norms2[i2], &small_norms1[i1]) == 0 {
                continue;
            }
            ibz_mul(v, &remain, &adjusted_norm);
            let vv = v.clone();
            ibz_mod(v, &vv, &small_norms1[i1]);
            let mut cmp = ibz_cmp(v, &quotients[i2]);
            while found == 0 && cmp < 0 {
                if number_sum_square > 0 {
                    found = ibz_cornacchia_prime(av, bv, ibz_const_one(), v) as i32;
                } else if number_sum_square == 0 {
                    found = 1;
                }
                if found != 0 {
                    ibz_mul(&mut remain, v, &small_norms2[i2]);
                    (au).clone_from(&n);
                    let au_v = au.clone();
                    ibz_sub(u, &au_v, &remain);
                    debug_assert!(ibz_cmp(u, ibz_const_zero()) > 0);
                    let uv = u.clone();
                    ibz_div(u, &mut remain, &uv, &small_norms1[i1]);
                    debug_assert!(ibz_is_zero(&remain));
                    found = i32::from(found != 0 && ibz_get(u) != 0 && ibz_get(v) != 0);
                    if number_sum_square == 2 {
                        found = ibz_cornacchia_prime(au, bu, ibz_const_one(), u) as i32;
                    }
                }
                if found == 0 {
                    let vv = v.clone();
                    ibz_add(v, &vv, &small_norms1[i1]);
                    cmp = ibz_cmp(v, &quotients[i2]);
                }
            }
            if found != 0 {
                *index_sol1 = i1;
                *index_sol2 = i2;
                break 'outer;
            }
        }
    }
    found
}

/// SuitableIdeals from the spec: find u, v, β₁, β₂, d₁, d₂ such that
/// u·d₁ + v·d₂ = target and βᵢ ∈ Jᵢ with reduced norm dᵢ·n(lideal).
#[allow(clippy::too_many_arguments)]
pub fn find_uv(
    u: &mut Ibz,
    v: &mut Ibz,
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    d1: &mut Ibz,
    d2: &mut Ibz,
    index_alternate_order_1: &mut usize,
    index_alternate_order_2: &mut usize,
    target: &Ibz,
    lideal: &QuatLeftIdeal,
    bpoo: &QuatAlg,
    num_alternate_order: usize,
) -> i32 {
    let mut au = Ibz::default();
    let mut bu = Ibz::default();
    let mut av = Ibz::default();
    let mut bv = Ibz::default();
    let mut n = Ibz::default();
    let mut remain = Ibz::default();
    n.clone_from(target);

    let no = num_alternate_order + 1;
    let mut adjusted_norm: Vec<Ibz> = (0..no).map(|_| Ibz::default()).collect();
    let mut gram: Vec<IbzMat4x4> = (0..no).map(|_| ibz_mat_4x4_init()).collect();
    let mut reduced: Vec<IbzMat4x4> = (0..no).map(|_| ibz_mat_4x4_init()).collect();
    let mut ideal: Vec<QuatLeftIdeal> = (0..no).map(|_| QuatLeftIdeal::default()).collect();

    ideal[0].clone_from(lideal);
    {
        let id0 = ideal[0].clone();
        quat_lideal_reduce_basis(&mut reduced[0], &mut gram[0], &id0, bpoo);
    }
    ideal[0].lattice.basis.clone_from(&reduced[0].clone());
    ibz_set(&mut adjusted_norm[0], 1);
    adjusted_norm[0] *= &ideal[0].lattice.denom;
    adjusted_norm[0] *= &ideal[0].lattice.denom;
    {
        let norm0 = ideal[0].norm.clone();
        post_lll_basis_treatment(&mut gram[0], &mut reduced[0], &norm0, true);
    }

    let mut reduced_id = ideal[0].clone();
    let mut delta = QuatAlgElem::default();
    ibz_set(&mut delta.coord[0], 1);
    ibz_set(&mut delta.coord[1], 0);
    ibz_set(&mut delta.coord[2], 0);
    ibz_set(&mut delta.coord[3], 0);
    delta.denom.clone_from(&reduced_id.lattice.denom);
    let dc = delta.coord.clone();
    ibz_mat_4x4_eval(&mut delta.coord, &reduced[0], &dc);
    debug_assert!(quat_lattice_contains(None, &reduced_id.lattice, &delta) != 0);

    let dd = delta.clone();
    quat_alg_conj(&mut delta, &dd);
    delta.denom *= &ideal[0].norm;
    let lat = reduced_id.lattice.clone();
    quat_lattice_alg_elem_mul(&mut reduced_id.lattice, &lat, &delta, bpoo);
    reduced_id.norm.clone_from(&gram[0][0][0]);
    let rn = reduced_id.norm.clone();
    ibz_div(&mut reduced_id.norm, &mut remain, &rn, &adjusted_norm[0]);
    debug_assert!(ibz_cmp(&remain, ibz_const_zero()) == 0);

    let mut right_order = QuatLattice::default();
    let mut conj_ideal = QuatLeftIdeal::default();
    quat_lideal_conjugate_without_hnf_dyn(&mut conj_ideal, &mut right_order, &reduced_id, bpoo);

    for i in 1..no {
        let aci = alternate_connecting_ideal(i - 1).clone();
        quat_lideal_lideal_mul_reduced_dyn(
            &mut ideal[i],
            &mut gram[i],
            &conj_ideal,
            &aci,
            &right_order,
            bpoo,
        );
        let basis_i = ideal[i].lattice.basis.clone();
        reduced[i].clone_from(&basis_i);
        ibz_set(&mut adjusted_norm[i], 1);
        adjusted_norm[i] *= &ideal[i].lattice.denom;
        adjusted_norm[i] *= &ideal[i].lattice.denom;
        let norm_i = ideal[i].norm.clone();
        post_lll_basis_treatment(&mut gram[i], &mut reduced[i], &norm_i, false);
    }

    let m = FINDUV_BOX_SIZE as i32;
    let m4 = FINDUV_CUBE_SIZE;

    let mut small_vecs: Vec<Vec<IbzVec4>> = (0..no)
        .map(|_| (0..m4).map(|_| ibz_vec_4_init()).collect())
        .collect();
    let mut small_norms: Vec<Vec<Ibz>> = (0..no)
        .map(|_| (0..m4).map(|_| Ibz::default()).collect())
        .collect();
    let mut quotients: Vec<Vec<Ibz>> = (0..no)
        .map(|_| (0..m4).map(|_| Ibz::default()).collect())
        .collect();
    let mut indices: Vec<usize> = vec![0; no];

    for j in 0..no {
        let cnt = enumerate_hypercube(
            &mut small_vecs[j],
            &mut small_norms[j],
            m,
            &gram[j],
            &adjusted_norm[j],
        );
        // C uses `count - 1`; clamp to ≥ 0 for safety.
        indices[j] = cnt.max(0) as usize;

        // Stable sort by (norm, original_index), matching qsort with idx tiebreak.
        let mut order: Vec<usize> = (0..indices[j]).collect();
        order.sort_by(|&a, &b| {
            let c = ibz_cmp(&small_norms[j][a], &small_norms[j][b]);
            if c != 0 {
                if c < 0 {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        });
        let svs: Vec<IbzVec4> = order
            .iter()
            .map(|&i| core::mem::take(&mut small_vecs[j][i]))
            .collect();
        let sns: Vec<Ibz> = order
            .iter()
            .map(|&i| core::mem::take(&mut small_norms[j][i]))
            .collect();
        for (i, (sv, sn)) in svs.into_iter().zip(sns).enumerate() {
            small_vecs[j][i] = sv;
            small_norms[j][i] = sn;
        }
        #[cfg(debug_assertions)]
        for i in 1..indices[j] {
            debug_assert!(ibz_cmp(&small_norms[j][i - 1], &small_norms[j][i]) <= 0);
        }

        for i in 0..indices[j] {
            ibz_div(&mut quotients[j][i], &mut remain, &n, &small_norms[j][i]);
        }
    }

    let mut found = 0;
    let mut i1 = 0usize;
    let mut i2 = 0usize;
    'search: for j1 in 0..no {
        for j2 in j1..no {
            let is_diago = j1 == j2;
            found = find_uv_from_lists(
                &mut au,
                &mut bu,
                &mut av,
                &mut bv,
                u,
                v,
                &mut i1,
                &mut i2,
                target,
                &small_norms[j1],
                &small_norms[j2],
                &quotients[j2],
                indices[j1],
                indices[j2],
                is_diago,
                0,
            );
            if found != 0 {
                beta1.denom.clone_from(&ideal[j1].lattice.denom);
                beta2.denom.clone_from(&ideal[j2].lattice.denom);
                (d1).clone_from(&small_norms[j1][i1]);
                (d2).clone_from(&small_norms[j2][i2]);
                ibz_mat_4x4_eval(&mut beta1.coord, &reduced[j1], &small_vecs[j1][i1]);
                ibz_mat_4x4_eval(&mut beta2.coord, &reduced[j2], &small_vecs[j2][i2]);
                debug_assert!(quat_lattice_contains(None, &ideal[j1].lattice, beta1) != 0);
                debug_assert!(quat_lattice_contains(None, &ideal[j2].lattice, beta2) != 0);

                if j1 != 0 || j2 != 0 {
                    let dd = delta.denom.clone();
                    ibz_div(&mut delta.denom, &mut remain, &dd, &lideal.norm);
                    debug_assert!(ibz_cmp(&remain, ibz_const_zero()) == 0);
                    delta.denom *= &conj_ideal.norm;
                }
                if j1 != 0 {
                    let b = beta1.clone();
                    quat_alg_mul(beta1, &delta, &b, bpoo);
                    quat_alg_normalize(beta1);
                }
                if j2 != 0 {
                    let b = beta2.clone();
                    quat_alg_mul(beta2, &delta, &b, bpoo);
                    quat_alg_normalize(beta2);
                }
                if j1 != 0 {
                    let b = beta1.clone();
                    quat_alg_conj(beta1, &b);
                }
                if j2 != 0 {
                    let b = beta2.clone();
                    quat_alg_conj(beta2, &b);
                }

                #[cfg(debug_assertions)]
                {
                    let mut nrm = Ibz::default();
                    let mut nd = Ibz::default();
                    let mut chk = Ibz::default();
                    quat_alg_norm(&mut nrm, &mut nd, beta1, quatalg_pinfty());
                    debug_assert!(ibz_is_one(&nd));
                    ibz_mul(&mut chk, d1, &ideal[0].norm);
                    if j1 > 0 {
                        let c = chk.clone();
                        ibz_mul(&mut chk, &c, &alternate_connecting_ideal(j1 - 1).norm);
                    }
                    debug_assert!(ibz_cmp(&chk, &nrm) == 0);
                    quat_alg_norm(&mut nrm, &mut nd, beta2, quatalg_pinfty());
                    debug_assert!(ibz_is_one(&nd));
                    ibz_mul(&mut chk, d2, &ideal[0].norm);
                    if j2 > 0 {
                        let c = chk.clone();
                        ibz_mul(&mut chk, &c, &alternate_connecting_ideal(j2 - 1).norm);
                    }
                    debug_assert!(ibz_cmp(&chk, &nrm) == 0);
                    debug_assert!(quat_lattice_contains(None, &ideal[0].lattice, beta1) != 0);
                    debug_assert!(quat_lattice_contains(None, &ideal[0].lattice, beta2) != 0);
                }

                *index_alternate_order_1 = j1;
                *index_alternate_order_2 = j2;
                break 'search;
            }
        }
    }
    found
}

// ---------------------------------------------------------------------------
// CLAPOTIS

#[allow(clippy::too_many_arguments)]
pub fn dim2id2iso_ideal_to_isogeny_clapotis(
    beta1: &mut QuatAlgElem,
    beta2: &mut QuatAlgElem,
    u: &mut Ibz,
    v: &mut Ibz,
    d1: &mut Ibz,
    d2: &mut Ibz,
    codomain: &mut EcCurve,
    basis: &mut EcBasis,
    lideal: &QuatLeftIdeal,
    bpoo: &QuatAlg,
) -> i32 {
    let mut tmp = Ibz::default();
    let mut two_pow = Ibz::default();
    let mut test1 = Ibz::default();
    let mut theta = QuatAlgElem::default();

    let mut index_order1 = 0usize;
    let mut index_order2 = 0usize;

    let ret = find_uv(
        u,
        v,
        beta1,
        beta2,
        d1,
        d2,
        &mut index_order1,
        &mut index_order2,
        torsion_plus_2power(),
        lideal,
        bpoo,
        NUM_ALTERNATE_EXTREMAL_ORDERS,
    );
    if ret == 0 {
        return 0;
    }
    debug_assert!(ibz_is_odd(d1) && ibz_is_odd(d2));

    ibz_gcd(&mut tmp, u, v);
    debug_assert!(ibz_cmp(&tmp, ibz_const_zero()) != 0);
    let exp_gcd = ibz_two_adic(&tmp);
    let exp = TORSION_EVEN_POWER as i32 - exp_gcd;
    let uv = u.clone();
    ibz_div(u, &mut test1, &uv, &tmp);
    debug_assert!(ibz_cmp(&test1, ibz_const_zero()) == 0);
    let vv = v.clone();
    ibz_div(v, &mut test1, &vv, &tmp);
    debug_assert!(ibz_cmp(&test1, ibz_const_zero()) == 0);

    #[cfg(debug_assertions)]
    {
        let mut pow_check = Ibz::default();
        let mut tmp_check = Ibz::default();
        ibz_pow(&mut pow_check, ibz_const_two(), exp as u32);
        ibz_mul(&mut tmp_check, d1, u);
        pow_check -= &tmp_check;
        ibz_mul(&mut tmp_check, v, d2);
        pow_check -= &tmp_check;
        debug_assert!(ibz_cmp(&pow_check, ibz_const_zero()) == 0);
    }

    let mut e1 = curves_with_endomorphisms()[index_order1].curve;
    let bas1 = curves_with_endomorphisms()[index_order1].basis_even;
    let mut bas2 = curves_with_endomorphisms()[index_order2].basis_even;

    let mut e01 = ThetaCoupleCurve::default();
    let mut ker = ThetaKernelCouplePoints::default();
    let mut bas_u = EcBasis::default();

    // theta = β₂ · conj(β₁) / n(lideal).
    ibz_set(&mut theta.denom, 1);
    let b1 = beta1.clone();
    quat_alg_conj(&mut theta, &b1);
    let th = theta.clone();
    quat_alg_mul(&mut theta, beta2, &th, quatalg_pinfty());
    theta.denom *= &lideal.norm;

    let mut idealu = QuatLeftIdeal::default();
    let mut idealv = QuatLeftIdeal::default();
    let mut fu_codomain = ThetaCoupleCurve::default();
    let mut fv_codomain = ThetaCoupleCurve::default();
    let mut pushed = [ThetaCouplePoint::default(); 3];

    pushed[0] = ThetaCouplePoint {
        p1: bas1.p,
        p2: EcPoint::default(),
    };
    pushed[1] = ThetaCouplePoint {
        p1: bas1.q,
        p2: EcPoint::default(),
    };
    pushed[2] = ThetaCouplePoint {
        p1: bas1.pmq,
        p2: EcPoint::default(),
    };
    pushed[0].p2 = EcPoint::IDENTITY;
    pushed[1].p2 = EcPoint::IDENTITY;
    pushed[2].p2 = EcPoint::IDENTITY;

    let ret = fixed_degree_isogeny_and_eval(
        &mut idealu,
        u,
        true,
        &mut fu_codomain,
        &mut pushed,
        index_order1,
    );
    if ret == 0 {
        return 0;
    }
    debug_assert!(test_point_order_twof(
        &pushed[0].p1,
        &fu_codomain.e1,
        TORSION_EVEN_POWER as i32
    ));
    debug_assert!(test_point_order_twof(
        &pushed[0].p2,
        &fu_codomain.e2,
        TORSION_EVEN_POWER as i32
    ));

    bas_u.p = pushed[0].p1;
    bas_u.q = pushed[1].p1;
    bas_u.pmq = pushed[2].p1;

    ker.t1.p1 = bas_u.p;
    ker.t2.p1 = bas_u.q;
    ker.t1m2.p1 = bas_u.pmq;
    e01.e1 = fu_codomain.e1;

    pushed[0] = ThetaCouplePoint {
        p1: bas2.p,
        p2: EcPoint::default(),
    };
    pushed[1] = ThetaCouplePoint {
        p1: bas2.q,
        p2: EcPoint::default(),
    };
    pushed[2] = ThetaCouplePoint {
        p1: bas2.pmq,
        p2: EcPoint::default(),
    };
    pushed[0].p2 = EcPoint::IDENTITY;
    pushed[1].p2 = EcPoint::IDENTITY;
    pushed[2].p2 = EcPoint::IDENTITY;

    let ret = fixed_degree_isogeny_and_eval(
        &mut idealv,
        v,
        true,
        &mut fv_codomain,
        &mut pushed,
        index_order2,
    );
    if ret == 0 {
        return 0;
    }
    debug_assert!(test_point_order_twof(
        &pushed[0].p1,
        &fv_codomain.e1,
        TORSION_EVEN_POWER as i32
    ));
    debug_assert!(test_point_order_twof(
        &pushed[0].p2,
        &fv_codomain.e2,
        TORSION_EVEN_POWER as i32
    ));

    bas2.p = pushed[0].p1;
    bas2.q = pushed[1].p1;
    bas2.pmq = pushed[2].p1;

    // theta ← theta / (d₁ · n(connecting_ideal[index_order2])) mod 2^TORSION_EVEN_POWER
    ibz_pow(&mut two_pow, ibz_const_two(), TORSION_EVEN_POWER as u32);
    tmp.clone_from(d1);
    if index_order2 > 0 {
        let t = tmp.clone();
        ibz_mul(
            &mut tmp,
            &t,
            &alternate_connecting_ideal(index_order2 - 1).norm,
        );
    }
    let t = tmp.clone();
    ibz_invmod(&mut tmp, &t, &two_pow);
    for i in 0..4 {
        theta.coord[i] *= &tmp;
    }

    endomorphism_application_even_basis(
        &mut bas2,
        0,
        &fv_codomain.e1,
        &theta,
        TORSION_EVEN_POWER as i32,
    );
    debug_assert!(test_basis_order_twof(
        &bas2,
        &fv_codomain.e1,
        TORSION_EVEN_POWER as i32
    ));

    ker.t1.p2 = bas2.p;
    ker.t2.p2 = bas2.q;
    ker.t1m2.p2 = bas2.pmq;
    e01.e2 = fv_codomain.e1;

    let dbl_n = (TORSION_EVEN_POWER as i32 - exp) as u32;
    let s = ker.t1;
    double_couple_point_iter(&mut ker.t1, dbl_n, &s, &e01);
    let s = ker.t2;
    double_couple_point_iter(&mut ker.t2, dbl_n, &s, &e01);
    let s = ker.t1m2;
    double_couple_point_iter(&mut ker.t1m2, dbl_n, &s, &e01);

    debug_assert!(test_point_order_twof(&ker.t1.p1, &e01.e1, exp));
    debug_assert!(test_point_order_twof(&ker.t1m2.p2, &e01.e2, exp));
    debug_assert!(ibz_is_odd(u));
    debug_assert!(test_basis_order_twof(
        &bas_u,
        &e01.e1,
        TORSION_EVEN_POWER as i32
    ));

    pushed[0].p1 = bas_u.p;
    pushed[1].p1 = bas_u.q;
    pushed[2].p1 = bas_u.pmq;
    pushed[0].p2 = EcPoint::IDENTITY;
    pushed[1].p2 = EcPoint::IDENTITY;
    pushed[2].p2 = EcPoint::IDENTITY;

    let mut theta_codomain = ThetaCoupleCurve::default();
    if !theta_chain_compute_and_eval_randomized(
        exp as u32,
        &mut e01,
        &ker,
        false,
        &mut theta_codomain,
        &mut pushed,
    ) {
        return 0;
    }

    let t1 = pushed[0];
    let t2 = pushed[1];
    let t1m2 = pushed[2];
    debug_assert!(test_point_order_twof(
        &t1.p2,
        &theta_codomain.e2,
        TORSION_EVEN_POWER as i32
    ));
    debug_assert!(test_point_order_twof(
        &t1.p1,
        &theta_codomain.e1,
        TORSION_EVEN_POWER as i32
    ));
    debug_assert!(test_point_order_twof(
        &t1m2.p2,
        &theta_codomain.e2,
        TORSION_EVEN_POWER as i32
    ));

    basis.p = t1.p1;
    basis.q = t2.p1;
    basis.pmq = t1m2.p1;
    *codomain = theta_codomain.e1;

    // Weil-pairing check to pick the correct factor.
    let mut w0 = Fp2::default();
    let mut w1 = Fp2::default();
    let mut codomain_tmp = *codomain;
    weil(
        &mut w0,
        TORSION_EVEN_POWER as u32,
        &bas1.p,
        &bas1.q,
        &bas1.pmq,
        &mut e1,
    );
    weil(
        &mut w1,
        TORSION_EVEN_POWER as u32,
        &basis.p,
        &basis.q,
        &basis.pmq,
        &mut codomain_tmp,
    );

    let mut digit_d = [0 as Digit; NWORDS_ORDER];
    ibz_mul(&mut tmp, d1, u);
    tmp *= &*u;
    let t = tmp.clone();
    ibz_mod(&mut tmp, &t, torsion_plus_2power());
    ibz_to_digits(&mut digit_d, &tmp);
    let test_pow = w0.pow_vartime(&digit_d);

    if w1 != test_pow {
        basis.p = t1.p2;
        basis.q = t2.p2;
        basis.pmq = t1m2.p2;
        *codomain = theta_codomain.e2;
        #[cfg(debug_assertions)]
        {
            let mut codomain_tmp2 = *codomain;
            weil(
                &mut w1,
                TORSION_EVEN_POWER as u32,
                &basis.p,
                &basis.q,
                &basis.pmq,
                &mut codomain_tmp2,
            );
            debug_assert!(w0.pow_vartime(&digit_d) == w1);
        }
    }

    // Apply M / (u·d₁) where M ~ β₁.
    ibz_mul(&mut tmp, u, d1);
    if index_order1 != 0 {
        let t = tmp.clone();
        ibz_mul(&mut tmp, &t, &connecting_ideals()[index_order1].norm);
    }
    let t = tmp.clone();
    ibz_invmod(&mut tmp, &t, torsion_plus_2power());
    for i in 0..4 {
        beta1.coord[i] *= &tmp;
    }
    endomorphism_application_even_basis(basis, 0, codomain, beta1, TORSION_EVEN_POWER as i32);

    1
}

/// Convenience wrapper: compute `basis` = image of canonical E₀ basis
/// under the isogeny corresponding to `lideal`, and `codomain` = its
/// codomain curve.
pub fn dim2id2iso_arbitrary_isogeny_evaluation(
    basis: &mut EcBasis,
    codomain: &mut EcCurve,
    lideal: &QuatLeftIdeal,
) -> i32 {
    let mut beta1 = QuatAlgElem::default();
    let mut beta2 = QuatAlgElem::default();
    let (mut u, mut v, mut d1, mut d2) = (
        Ibz::default(),
        Ibz::default(),
        Ibz::default(),
        Ibz::default(),
    );
    dim2id2iso_ideal_to_isogeny_clapotis(
        &mut beta1,
        &mut beta2,
        &mut u,
        &mut v,
        &mut d1,
        &mut d2,
        codomain,
        basis,
        lideal,
        quatalg_pinfty(),
    )
}
