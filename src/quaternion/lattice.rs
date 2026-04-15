//! Full-rank dimension-4 lattices, ported from `lattice.c`.

use super::algebra::*;
use super::dim4::*;
use super::hnf::*;
use super::intbig::*;
use super::types::*;

pub fn quat_lattice_reduce_denom(reduced: &mut QuatLattice, lat: &QuatLattice) {
    let mut gcd = Ibz::default();
    ibz_mat_4x4_gcd(&mut gcd, &lat.basis);
    let g = gcd.clone();
    ibz_gcd(&mut gcd, &g, &lat.denom);
    ibz_mat_4x4_scalar_div(&mut reduced.basis, &gcd, &lat.basis);
    let mut r = Ibz::default();
    ibz_div(&mut reduced.denom, &mut r, &lat.denom, &gcd);
    let d = reduced.denom.clone();
    ibz_abs(&mut reduced.denom, &d);
}

impl PartialEq for QuatLattice {
    fn eq(&self, other: &Self) -> bool {
        quat_lattice_equal(self, other) != 0
    }
}
impl Eq for QuatLattice {}

pub fn quat_lattice_equal(lat1: &QuatLattice, lat2: &QuatLattice) -> i32 {
    let mut a = QuatLattice::default();
    let mut b = QuatLattice::default();
    quat_lattice_reduce_denom(&mut a, lat1);
    quat_lattice_reduce_denom(&mut b, lat2);
    let d = a.denom.clone();
    ibz_abs(&mut a.denom, &d);
    let d = b.denom.clone();
    ibz_abs(&mut b.denom, &d);
    quat_lattice_hnf(&mut a);
    quat_lattice_hnf(&mut b);
    ((ibz_cmp(&a.denom, &b.denom) == 0) && ibz_mat_4x4_equal(&a.basis, &b.basis) != 0) as i32
}

pub fn quat_lattice_inclusion(sublat: &QuatLattice, overlat: &QuatLattice) -> i32 {
    let mut sum = QuatLattice::default();
    quat_lattice_add(&mut sum, overlat, sublat);
    quat_lattice_equal(&sum, overlat)
}

pub fn quat_lattice_conjugate_without_hnf(conj: &mut QuatLattice, lat: &QuatLattice) {
    conj.basis.clone_from(&lat.basis);
    conj.denom.clone_from(&lat.denom);
    for row in 1..4 {
        for col in 0..4 {
            let t = conj.basis[row][col].clone();
            ibz_neg(&mut conj.basis[row][col], &t);
        }
    }
}

pub fn quat_lattice_dual_without_hnf(dual: &mut QuatLattice, lat: &QuatLattice) {
    let mut inv = ibz_mat_4x4_init();
    let mut det = Ibz::default();
    ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), &mut det, &lat.basis);
    let invc = inv.clone();
    ibz_mat_4x4_transpose(&mut inv, &invc);
    ibz_mat_4x4_scalar_mul(&mut dual.basis, &lat.denom, &inv);
    dual.denom.clone_from(&det);
}

pub fn quat_lattice_add(res: &mut QuatLattice, lat1: &QuatLattice, lat2: &QuatLattice) {
    let mut generators: Vec<IbzVec4> = (0..8).map(|_| ibz_vec_4_init()).collect();
    let mut tmp = ibz_mat_4x4_init();
    let mut det1 = Ibz::default();
    let mut det2 = Ibz::default();
    let mut detprod = Ibz::default();

    ibz_mat_4x4_scalar_mul(&mut tmp, &lat1.denom, &lat2.basis);
    for i in 0..4 {
        for j in 0..4 {
            generators[j][i].clone_from(&tmp[i][j]);
        }
    }
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det1, &tmp);
    ibz_mat_4x4_scalar_mul(&mut tmp, &lat2.denom, &lat1.basis);
    for i in 0..4 {
        for j in 0..4 {
            generators[4 + j][i].clone_from(&tmp[i][j]);
        }
    }
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det2, &tmp);
    debug_assert!(!ibz_is_zero(&det1));
    debug_assert!(!ibz_is_zero(&det2));
    ibz_gcd(&mut detprod, &det1, &det2);
    ibz_mat_4xn_hnf_mod_core(&mut res.basis, 8, &generators, &detprod);
    ibz_mul(&mut res.denom, &lat1.denom, &lat2.denom);
    let r = res.clone();
    quat_lattice_reduce_denom(res, &r);
}

pub fn quat_lattice_intersect(res: &mut QuatLattice, lat1: &QuatLattice, lat2: &QuatLattice) {
    let mut dual1 = QuatLattice::default();
    let mut dual2 = QuatLattice::default();
    let mut dual_res = QuatLattice::default();
    quat_lattice_dual_without_hnf(&mut dual1, lat1);
    quat_lattice_dual_without_hnf(&mut dual2, lat2);
    quat_lattice_add(&mut dual_res, &dual1, &dual2);
    quat_lattice_dual_without_hnf(res, &dual_res);
    quat_lattice_hnf(res);
}

pub fn quat_lattice_mat_alg_coord_mul_without_hnf(
    prod: &mut IbzMat4x4,
    lat: &IbzMat4x4,
    coord: &IbzVec4,
    alg: &QuatAlg,
) {
    let mut p = ibz_vec_4_init();
    let mut a = ibz_vec_4_init();
    for i in 0..4 {
        ibz_vec_4_copy_ibz(&mut a, &lat[0][i], &lat[1][i], &lat[2][i], &lat[3][i]);
        quat_alg_coord_mul(&mut p, &a, coord, alg);
        for r in 0..4 {
            prod[r][i].clone_from(&p[r]);
        }
    }
}

pub fn quat_lattice_alg_elem_mul(
    prod: &mut QuatLattice,
    lat: &QuatLattice,
    elem: &QuatAlgElem,
    alg: &QuatAlg,
) {
    quat_lattice_mat_alg_coord_mul_without_hnf(&mut prod.basis, &lat.basis, &elem.coord, alg);
    ibz_mul(&mut prod.denom, &lat.denom, &elem.denom);
    quat_lattice_hnf(prod);
}

pub fn quat_lattice_mul(
    res: &mut QuatLattice,
    lat1: &QuatLattice,
    lat2: &QuatLattice,
    alg: &QuatAlg,
) {
    let mut elem1 = ibz_vec_4_init();
    let mut elem2 = ibz_vec_4_init();
    let mut elem_res = ibz_vec_4_init();
    let mut generators: Vec<IbzVec4> = (0..16).map(|_| ibz_vec_4_init()).collect();
    let mut detmat = ibz_mat_4x4_init();
    let mut det = Ibz::default();
    for k in 0..4 {
        ibz_vec_4_copy_ibz(
            &mut elem1,
            &lat1.basis[0][k],
            &lat1.basis[1][k],
            &lat1.basis[2][k],
            &lat1.basis[3][k],
        );
        for i in 0..4 {
            ibz_vec_4_copy_ibz(
                &mut elem2,
                &lat2.basis[0][i],
                &lat2.basis[1][i],
                &lat2.basis[2][i],
                &lat2.basis[3][i],
            );
            quat_alg_coord_mul(&mut elem_res, &elem1, &elem2, alg);
            for j in 0..4 {
                if k == 0 {
                    detmat[i][j].clone_from(&elem_res[j]);
                }
                generators[4 * k + i][j].clone_from(&elem_res[j]);
            }
        }
    }
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &detmat);
    let d = det.clone();
    ibz_abs(&mut det, &d);
    ibz_mat_4xn_hnf_mod_core(&mut res.basis, 16, &generators, &det);
    ibz_mul(&mut res.denom, &lat1.denom, &lat2.denom);
    let r = res.clone();
    quat_lattice_reduce_denom(res, &r);
}

/// Returns 1 if `x ∈ lat`, writing its coordinates in `lat`'s basis to `coord`.
pub fn quat_lattice_contains(
    coord: Option<&mut IbzVec4>,
    lat: &QuatLattice,
    x: &QuatAlgElem,
) -> i32 {
    let mut work_coord = ibz_vec_4_init();
    let mut inv = ibz_mat_4x4_init();
    let mut det = Ibz::default();
    let mut prod = Ibz::default();
    ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), &mut det, &lat.basis);
    debug_assert!(!ibz_is_zero(&det));
    ibz_mat_4x4_eval(&mut work_coord, &inv, &x.coord);
    let wc = work_coord.clone();
    ibz_vec_4_scalar_mul(&mut work_coord, &lat.denom, &wc);
    ibz_mul(&mut prod, &x.denom, &det);
    let wc = work_coord.clone();
    let divisible = ibz_vec_4_scalar_div(&mut work_coord, &prod, &wc);
    if divisible {
        if let Some(coord) = coord {
            (coord).clone_from(&work_coord);
        }
    }
    divisible as i32
}

pub fn quat_lattice_index(index: &mut Ibz, sublat: &QuatLattice, overlat: &QuatLattice) {
    let mut tmp = Ibz::default();
    let mut det = Ibz::default();
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &sublat.basis);
    ibz_mul(&mut tmp, &overlat.denom, &overlat.denom);
    tmp *= tmp.clone();
    ibz_mul(index, &det, &tmp);
    ibz_mul(&mut tmp, &sublat.denom, &sublat.denom);
    tmp *= tmp.clone();
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &overlat.basis);
    tmp *= &det;
    let idx = index.clone();
    let tm = tmp.clone();
    ibz_div(index, &mut tmp, &idx, &tm);
    debug_assert!(ibz_is_zero(&tmp));
    let t = index.clone();
    ibz_abs(index, &t);
}

pub fn quat_lattice_hnf(lat: &mut QuatLattice) {
    let mut modn = Ibz::default();
    ibz_mat_4x4_inv_with_det_as_denom(None, &mut modn, &lat.basis);
    let m = modn.clone();
    ibz_abs(&mut modn, &m);
    let mut generators: [IbzVec4; 4] = Default::default();
    for i in 0..4 {
        for j in 0..4 {
            generators[j][i].clone_from(&lat.basis[i][j]);
        }
    }
    ibz_mat_4xn_hnf_mod_core(&mut lat.basis, 4, &generators, &modn);
    let l = lat.clone();
    quat_lattice_reduce_denom(lat, &l);
}

/// Gram matrix of `lattice` for the bilinear form `2·(x₀y₀ + x₁y₁ + p·x₂y₂ + p·x₃y₃)`.
pub fn quat_lattice_gram(g: &mut IbzMat4x4, lattice: &QuatLattice, alg: &QuatAlg) {
    let mut tmp = Ibz::default();
    for i in 0..4 {
        for j in 0..=i {
            ibz_set(&mut g[i][j], 0);
            for k in 0..4 {
                ibz_mul(&mut tmp, &lattice.basis[k][i], &lattice.basis[k][j]);
                if k >= 2 {
                    tmp *= &alg.p;
                }
                g[i][j] += &tmp;
            }
            let t = g[i][j].clone();
            ibz_mul(&mut g[i][j], &t, ibz_const_two());
        }
    }
    for i in 0..4 {
        for j in (i + 1)..4 {
            let (lo, hi) = g.split_at_mut(j);
            lo[i][j].clone_from(&hi[0][i]);
        }
    }
}

/// Finishes the deferred port from `algebra.c`: returns coords of `x` in
/// `order` divided by their content.
pub fn quat_alg_make_primitive(
    primitive_x: &mut IbzVec4,
    content: &mut Ibz,
    x: &QuatAlgElem,
    order: &QuatLattice,
) {
    let ok = quat_lattice_contains(Some(primitive_x), order, x);
    debug_assert!(ok != 0);
    ibz_vec_4_content(content, primitive_x);
    let mut r = Ibz::default();
    for i in 0..4 {
        let t = primitive_x[i].clone();
        ibz_div(&mut primitive_x[i], &mut r, &t, content);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::intbig::ibz_from_i64;

    fn z(v: i64) -> Ibz {
        ibz_from_i64(v)
    }

    fn lat_identity() -> QuatLattice {
        let mut l = QuatLattice::default();
        ibz_mat_4x4_identity(&mut l.basis);
        l
    }

    #[test]
    fn add_intersect_index() {
        let id = lat_identity();
        let mut l2 = QuatLattice::default();
        for i in 0..4 {
            ibz_set(&mut l2.basis[i][i], 2);
        }
        // 2Z^4 + Z^4 = Z^4
        let mut sum = QuatLattice::default();
        quat_lattice_add(&mut sum, &l2, &id);
        assert_eq!(quat_lattice_equal(&sum, &id), 1);
        // Z^4 ∩ 2Z^4 = 2Z^4
        let mut inter = QuatLattice::default();
        quat_lattice_intersect(&mut inter, &id, &l2);
        assert_eq!(quat_lattice_equal(&inter, &l2), 1);
        // [Z^4 : 2Z^4] = 16
        let mut idx = Ibz::default();
        quat_lattice_index(&mut idx, &l2, &id);
        assert_eq!(idx, z(16));
        assert_eq!(quat_lattice_inclusion(&l2, &id), 1);
        assert_eq!(quat_lattice_inclusion(&id, &l2), 0);
    }

    #[test]
    fn contains_and_make_primitive() {
        let id = lat_identity();
        let mut x = QuatAlgElem::default();
        quat_alg_elem_set(&mut x, 1, 4, 6, 10, 14);
        let mut coord = ibz_vec_4_init();
        assert_eq!(quat_lattice_contains(Some(&mut coord), &id, &x), 1);
        assert_eq!(coord, [z(4), z(6), z(10), z(14)]);
        let mut prim = ibz_vec_4_init();
        let mut content = Ibz::default();
        quat_alg_make_primitive(&mut prim, &mut content, &x, &id);
        assert_eq!(content, z(2));
        assert_eq!(prim, [z(2), z(3), z(5), z(7)]);
        // Not contained when denom doesn't divide.
        let mut y = QuatAlgElem::default();
        quat_alg_elem_set(&mut y, 3, 1, 0, 0, 0);
        assert_eq!(quat_lattice_contains(None, &id, &y), 0);
    }

    #[test]
    fn gram_matrix() {
        let alg = QuatAlg::from_ui(7);
        let id = lat_identity();
        let mut g = ibz_mat_4x4_init();
        quat_lattice_gram(&mut g, &id, &alg);
        assert_eq!(g[0][0], z(2));
        assert_eq!(g[1][1], z(2));
        assert_eq!(g[2][2], z(14));
        assert_eq!(g[3][3], z(14));
        assert_eq!(g[0][1], z(0));
    }

    #[test]
    fn hnf_idempotent() {
        let mut l = QuatLattice::default();
        // Non-HNF basis of Z^4.
        l.basis = core::array::from_fn(|i| core::array::from_fn(|j| z(if j >= i { 1 } else { 0 })));
        quat_lattice_hnf(&mut l);
        let l1 = l.clone();
        quat_lattice_hnf(&mut l);
        assert_eq!(ibz_mat_4x4_equal(&l.basis, &l1.basis), 1);
        let mut det = Ibz::default();
        ibz_mat_4x4_inv_with_det_as_denom(None, &mut det, &l.basis);
        assert_eq!(det, z(1));
    }
}
