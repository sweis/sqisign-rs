// SPDX-License-Identifier: Apache-2.0
//! SQIsign key generation and signing.
//! Ported from `src/signature/ref/lvlx/{keygen,sign,encode_signature}.c` and
//! the top-level wrappers in `src/sqisign.c`.

use std::sync::OnceLock;

use crate::ec::*;
use crate::gf::{fp2_copy, fp2_is_one, Fp2};
use crate::hd::*;
use crate::id2iso::*;
use crate::mp::Digit;
use crate::params::*;
use crate::precomp::*;
use crate::quaternion::*;
use crate::verification::{
    hash_to_challenge, public_key_from_bytes, public_key_to_bytes, signature_to_bytes, PublicKey,
    Signature,
};

#[cfg(debug_assertions)]
use crate::ec::biextension::weil;

// ---------------------------------------------------------------------------
// Torsion-constant accessors as Ibz (the precomp module stores raw limb arrays
// to keep the verify path GMP-free; convert lazily here).
// ---------------------------------------------------------------------------

fn ibz_from_limbs(limbs: &[u64]) -> Ibz {
    let mut x = Ibz::default();
    ibz_copy_digits(&mut x, limbs);
    x
}

fn sec_degree() -> &'static Ibz {
    static V: OnceLock<Ibz> = OnceLock::new();
    V.get_or_init(|| ibz_from_limbs(SEC_DEGREE))
}
fn com_degree() -> &'static Ibz {
    static V: OnceLock<Ibz> = OnceLock::new();
    V.get_or_init(|| ibz_from_limbs(COM_DEGREE))
}
fn torsion_plus_2power_ibz() -> &'static Ibz {
    static V: OnceLock<Ibz> = OnceLock::new();
    V.get_or_init(|| ibz_from_limbs(TORSION_PLUS_2POWER))
}

// ---------------------------------------------------------------------------
// Secret key type
// ---------------------------------------------------------------------------

#[derive(Clone, Default)]
pub struct SecretKey {
    pub curve: EcCurve,
    pub secret_ideal: QuatLeftIdeal,
    pub mat_ba_can_to_ba0_two: IbzMat2x2,
    pub canonical_basis: EcBasis,
}

impl core::fmt::Debug for SecretKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SecretKey").finish_non_exhaustive()
    }
}

// Best-effort zeroization. `rug::Integer` (inside QuatLeftIdeal/IbzMat2x2) does
// not zero its heap allocation on drop, so the ideal/matrix limbs may persist.
// See PORTING.md "Security hardening" for details.
impl Drop for SecretKey {
    fn drop(&mut self) {
        use zeroize::Zeroize;
        self.curve.a.re.0.zeroize();
        self.curve.a.im.0.zeroize();
        self.curve.c.re.0.zeroize();
        self.curve.c.im.0.zeroize();
        for p in [
            &mut self.canonical_basis.p,
            &mut self.canonical_basis.q,
            &mut self.canonical_basis.pmq,
        ] {
            p.x.re.0.zeroize();
            p.x.im.0.zeroize();
            p.z.re.0.zeroize();
            p.z.im.0.zeroize();
        }
        // Scrub bigint limb buffers. The rug and crypto-bigint backends scrub the
        // current allocation; the malachite backend is best-effort (no public
        // mutable limb access). Buffers from prior reallocations are never
        // covered — see PORTING.md.
        ibz_secure_clear(&mut self.secret_ideal.norm);
        for row in &mut self.secret_ideal.lattice.basis {
            for e in row.iter_mut() {
                ibz_secure_clear(e);
            }
        }
        ibz_secure_clear(&mut self.secret_ideal.lattice.denom);
        for row in &mut self.mat_ba_can_to_ba0_two {
            for e in row.iter_mut() {
                ibz_secure_clear(e);
            }
        }
    }
}

pub fn secret_key_init(sk: &mut SecretKey) {
    *sk = SecretKey::default();
    ec_curve_init(&mut sk.curve);
}

// ---------------------------------------------------------------------------
// keygen.c
// ---------------------------------------------------------------------------

pub fn protocols_keygen(pk: &mut PublicKey, sk: &mut SecretKey) -> i32 {
    let mut found = 0;
    let mut b0_two = EcBasis::default();
    let params = quat_represent_integer_params();

    while found == 0 {
        found = quat_sampling_random_ideal_o0_given_norm(
            &mut sk.secret_ideal,
            sec_degree(),
            true,
            &params,
            None,
        );

        found = (found != 0
            && quat_lideal_prime_norm_reduced_equivalent(
                &mut sk.secret_ideal,
                quatalg_pinfty(),
                QUAT_PRIMALITY_NUM_ITER as i32,
                QUAT_EQUIV_BOUND_COEFF as i32,
            ) != 0) as i32;

        found = (found != 0
            && dim2id2iso_arbitrary_isogeny_evaluation(
                &mut b0_two,
                &mut sk.curve,
                &sk.secret_ideal,
            ) != 0) as i32;
    }

    debug_assert!(test_basis_order_twof(&b0_two, &sk.curve, TORSION_EVEN_POWER as i32) != 0);

    pk.hint_pk = ec_curve_to_basis_2f_to_hint(
        &mut sk.canonical_basis,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
    );

    debug_assert!(
        test_basis_order_twof(&sk.canonical_basis, &sk.curve, TORSION_EVEN_POWER as i32) != 0
    );

    change_of_basis_matrix_tate(
        &mut sk.mat_ba_can_to_ba0_two,
        &sk.canonical_basis,
        &b0_two,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
    );

    copy_curve(&mut pk.curve, &sk.curve);
    pk.curve.is_a24_computed_and_normalized = false;
    debug_assert!(fp2_is_one(&pk.curve.c) == 0xFFFF_FFFF);

    found
}

// ---------------------------------------------------------------------------
// encode_signature.c
// ---------------------------------------------------------------------------

use crate::mp::{decode_digits, encode_digits};

fn ibz_to_bytes(enc: &mut [u8], x: &Ibz, nbytes: usize, sgn: bool) -> usize {
    #[cfg(debug_assertions)]
    {
        let mut bnd = Ibz::default();
        let mut absv = Ibz::default();
        ibz_pow(
            &mut bnd,
            ibz_const_two(),
            (8 * nbytes - sgn as usize) as u32,
        );
        ibz_abs(&mut absv, x);
        debug_assert!(ibz_cmp(&absv, &bnd) < 0);
    }
    let ndigits = nbytes.div_ceil(8);
    let mut d = zeroize::Zeroizing::new(vec![0u64; ndigits]);
    if ibz_cmp(x, ibz_const_zero()) >= 0 {
        ibz_to_digits(&mut d, x);
    } else {
        debug_assert!(sgn);
        let mut tmp = Ibz::default();
        ibz_neg(&mut tmp, x);
        let s = tmp.clone();
        ibz_sub(&mut tmp, &s, ibz_const_one());
        ibz_to_digits(&mut d, &tmp);
        ibz_secure_clear(&mut tmp);
        for di in d.iter_mut() {
            *di = !*di;
        }
    }
    encode_digits(enc, &d, nbytes);
    nbytes
}

fn ibz_from_bytes(x: &mut Ibz, enc: &[u8], nbytes: usize, sgn: bool) -> usize {
    debug_assert!(nbytes > 0);
    let ndigits = nbytes.div_ceil(8);
    let mut d = zeroize::Zeroizing::new(vec![0u64; ndigits]);
    decode_digits(&mut d, enc, nbytes);
    if sgn && (enc[nbytes - 1] >> 7) != 0 {
        let s = 7 - (ndigits * 8 - nbytes);
        debug_assert!(s < 8);
        d[ndigits - 1] |= ((!0u64) >> (8 * s)) << (8 * s);
        for di in d.iter_mut() {
            *di = !*di;
        }
        ibz_copy_digits(x, &d);
        let s = x.clone();
        ibz_add(x, &s, ibz_const_one());
        let s = x.clone();
        ibz_neg(x, &s);
    } else {
        ibz_copy_digits(x, &d);
    }
    nbytes
}

pub fn secret_key_to_bytes(enc: &mut [u8], sk: &SecretKey, pk: &PublicKey) {
    debug_assert_eq!(enc.len(), SECRETKEY_BYTES);
    let mut p = 0;
    public_key_to_bytes(&mut enc[..PUBLICKEY_BYTES], pk);
    p += PUBLICKEY_BYTES;

    p += ibz_to_bytes(
        &mut enc[p..],
        &sk.secret_ideal.norm,
        FP_ENCODED_BYTES,
        false,
    );

    {
        let mut gen = QuatAlgElem::default();
        let ok = quat_lideal_generator(&mut gen, &sk.secret_ideal, quatalg_pinfty());
        debug_assert!(ok != 0);
        let _ = ok;
        for k in 0..4 {
            p += ibz_to_bytes(&mut enc[p..], &gen.coord[k], FP_ENCODED_BYTES, true);
        }
    }

    for i in 0..2 {
        for j in 0..2 {
            p += ibz_to_bytes(
                &mut enc[p..],
                &sk.mat_ba_can_to_ba0_two[i][j],
                TORSION_2POWER_BYTES,
                false,
            );
        }
    }
    debug_assert_eq!(p, SECRETKEY_BYTES);
}

pub fn secret_key_from_bytes(sk: &mut SecretKey, pk: &mut PublicKey, enc: &[u8]) {
    debug_assert_eq!(enc.len(), SECRETKEY_BYTES);
    let mut p = 0;
    let ok = public_key_from_bytes(pk, &enc[..PUBLICKEY_BYTES]);
    debug_assert!(ok, "non-canonical pk encoding inside sk");
    p += PUBLICKEY_BYTES;

    {
        let mut norm = Ibz::default();
        let mut gen = QuatAlgElem::default();
        p += ibz_from_bytes(&mut norm, &enc[p..], FP_ENCODED_BYTES, false);
        for k in 0..4 {
            p += ibz_from_bytes(&mut gen.coord[k], &enc[p..], FP_ENCODED_BYTES, true);
        }
        quat_lideal_create(
            &mut sk.secret_ideal,
            &gen,
            &norm,
            maxord_o0(),
            quatalg_pinfty(),
        );
    }

    for i in 0..2 {
        for j in 0..2 {
            p += ibz_from_bytes(
                &mut sk.mat_ba_can_to_ba0_two[i][j],
                &enc[p..],
                TORSION_2POWER_BYTES,
                false,
            );
        }
    }
    debug_assert_eq!(p, SECRETKEY_BYTES);

    sk.curve = pk.curve;
    ec_curve_to_basis_2f_from_hint(
        &mut sk.canonical_basis,
        &mut sk.curve,
        TORSION_EVEN_POWER as i32,
        pk.hint_pk,
    );
}

// ---------------------------------------------------------------------------
// sign.c
// ---------------------------------------------------------------------------

fn commit(
    e_com: &mut EcCurve,
    basis_even_com: &mut EcBasis,
    lideal_com: &mut QuatLeftIdeal,
) -> bool {
    let params = quat_represent_integer_params();
    let mut found =
        quat_sampling_random_ideal_o0_given_norm(lideal_com, com_degree(), true, &params, None)
            != 0;
    found = found
        && quat_lideal_prime_norm_reduced_equivalent(
            lideal_com,
            quatalg_pinfty(),
            QUAT_PRIMALITY_NUM_ITER as i32,
            QUAT_EQUIV_BOUND_COEFF as i32,
        ) != 0;
    found =
        found && dim2id2iso_arbitrary_isogeny_evaluation(basis_even_com, e_com, lideal_com) != 0;
    found
}

fn compute_challenge_ideal_signature(
    lideal_chall_two: &mut QuatLeftIdeal,
    sig: &Signature,
    sk: &SecretKey,
) {
    let mut vec = ibz_vec_2_init();
    ibz_set(&mut vec[0], 1);
    ibz_copy_digits(&mut vec[1], &sig.chall_coeff);

    let v_in = vec.clone();
    ibz_mat_2x2_eval(&mut vec, &sk.mat_ba_can_to_ba0_two, &v_in);

    id2iso_kernel_dlogs_to_ideal_even(lideal_chall_two, &vec, TORSION_EVEN_POWER as i32);
    debug_assert_eq!(
        ibz_cmp(&lideal_chall_two.norm, torsion_plus_2power_ibz()),
        0
    );
}

fn sample_response(x: &mut QuatAlgElem, lattice: &QuatLattice, lattice_content: &Ibz) {
    let mut bound = Ibz::default();
    ibz_pow(&mut bound, ibz_const_two(), SQISIGN_RESPONSE_LENGTH as u32);
    let s = bound.clone();
    ibz_sub(&mut bound, &s, ibz_const_one());
    let s = bound.clone();
    ibz_mul(&mut bound, &s, lattice_content);

    let ok = quat_lattice_sample_from_ball(x, lattice, quatalg_pinfty(), &bound);
    debug_assert!(ok != 0);
    let _ = ok;
}

fn compute_response_quat_element(
    resp_quat: &mut QuatAlgElem,
    lattice_content: &mut Ibz,
    sk: &SecretKey,
    lideal_chall_two: &QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
) {
    let mut lideal_chall_secret = QuatLeftIdeal::default();
    let mut lat_commit = QuatLattice::default();
    let mut lattice_hom_chall_to_com = QuatLattice::default();

    quat_lideal_inter(
        &mut lideal_chall_secret,
        lideal_chall_two,
        &sk.secret_ideal,
        quatalg_pinfty(),
    );

    quat_lattice_conjugate_without_hnf(&mut lat_commit, &lideal_commit.lattice);
    quat_lattice_intersect(
        &mut lattice_hom_chall_to_com,
        &lideal_chall_secret.lattice,
        &lat_commit,
    );

    ibz_mul(
        lattice_content,
        &lideal_chall_secret.norm,
        &lideal_commit.norm,
    );
    sample_response(resp_quat, &lattice_hom_chall_to_com, lattice_content);
}

fn compute_backtracking_signature(
    sig: &mut Signature,
    resp_quat: &mut QuatAlgElem,
    lattice_content: &mut Ibz,
    remain: &mut Ibz,
) {
    let mut tmp = Ibz::default();
    let mut dummy_coord = ibz_vec_4_init();

    quat_alg_make_primitive(&mut dummy_coord, &mut tmp, resp_quat, maxord_o0());
    let s = resp_quat.denom.clone();
    ibz_mul(&mut resp_quat.denom, &s, &tmp);
    debug_assert!(quat_lattice_contains(None, maxord_o0(), resp_quat) != 0);

    let backtracking = ibz_two_adic(&tmp);
    sig.backtracking = backtracking as u8;

    ibz_pow(&mut tmp, ibz_const_two(), backtracking as u32);
    let lc = lattice_content.clone();
    ibz_div(lattice_content, remain, &lc, &tmp);
}

fn compute_random_aux_norm_and_helpers(
    sig: &mut Signature,
    random_aux_norm: &mut Ibz,
    degree_resp_inv: &mut Ibz,
    remain: &mut Ibz,
    lattice_content: &Ibz,
    resp_quat: &mut QuatAlgElem,
    lideal_com_resp: &mut QuatLeftIdeal,
    lideal_commit: &QuatLeftIdeal,
) -> u8 {
    let mut degree_full_resp = Ibz::default();
    let mut degree_odd_resp = Ibz::default();
    let mut norm_d = Ibz::default();
    let mut tmp = Ibz::default();

    quat_alg_norm(
        &mut degree_full_resp,
        &mut norm_d,
        resp_quat,
        quatalg_pinfty(),
    );
    debug_assert!(ibz_is_one(&norm_d) != 0);
    let dfr = degree_full_resp.clone();
    ibz_div(&mut degree_full_resp, remain, &dfr, lattice_content);
    debug_assert_eq!(ibz_cmp(remain, ibz_const_zero()), 0);

    let exp_diadic_val_full_resp = ibz_two_adic(&degree_full_resp);
    sig.two_resp_length = exp_diadic_val_full_resp as u8;

    ibz_pow(&mut tmp, ibz_const_two(), exp_diadic_val_full_resp as u32);
    ibz_div(&mut degree_odd_resp, remain, &degree_full_resp, &tmp);
    debug_assert_eq!(ibz_cmp(remain, ibz_const_zero()), 0);

    let rq = resp_quat.clone();
    quat_alg_conj(resp_quat, &rq);

    ibz_mul(&mut tmp, &lideal_commit.norm, &degree_odd_resp);
    quat_lideal_create(
        lideal_com_resp,
        resp_quat,
        &tmp,
        maxord_o0(),
        quatalg_pinfty(),
    );

    let pow_dim2_deg_resp =
        SQISIGN_RESPONSE_LENGTH as i32 - exp_diadic_val_full_resp - sig.backtracking as i32;
    ibz_pow(remain, ibz_const_two(), pow_dim2_deg_resp as u32);
    ibz_sub(random_aux_norm, remain, &degree_odd_resp);

    ibz_pow(
        remain,
        ibz_const_two(),
        (pow_dim2_deg_resp + HD_EXTRA_TORSION as i32) as u32,
    );
    ibz_invmod(degree_resp_inv, &degree_odd_resp, remain);

    pow_dim2_deg_resp as u8
}

fn evaluate_random_aux_isogeny_signature(
    e_aux: &mut EcCurve,
    b_aux: &mut EcBasis,
    norm: &Ibz,
    lideal_com_resp: &QuatLeftIdeal,
) -> i32 {
    let mut lideal_aux = QuatLeftIdeal::default();
    let mut lideal_aux_resp_com = QuatLeftIdeal::default();
    let params = quat_represent_integer_params();

    let mut found = quat_sampling_random_ideal_o0_given_norm(
        &mut lideal_aux,
        norm,
        false,
        &params,
        Some(quat_prime_cofactor()),
    );
    if found != 0 {
        quat_lideal_inter(
            &mut lideal_aux_resp_com,
            lideal_com_resp,
            &lideal_aux,
            quatalg_pinfty(),
        );
        found = dim2id2iso_arbitrary_isogeny_evaluation(b_aux, e_aux, &lideal_aux_resp_com);
    }
    found
}

fn compute_dim2_isogeny_challenge(
    codomain: &mut ThetaCoupleCurveWithBasis,
    domain: &ThetaCoupleCurveWithBasis,
    degree_resp_inv: &Ibz,
    pow_dim2_deg_resp: i32,
    exp_diadic_val_full_resp: i32,
    reduced_order: i32,
) -> i32 {
    let mut e_com_x_aux = ThetaCoupleCurve::default();
    copy_curve(&mut e_com_x_aux.e1, &domain.e1);
    copy_curve(&mut e_com_x_aux.e2, &domain.e2);

    let mut dim_two_ker = ThetaKernelCouplePoints::default();
    copy_bases_to_kernel(&mut dim_two_ker, &domain.b1, &domain.b2);

    let mut scalar = [0 as Digit; NWORDS_ORDER];
    ibz_to_digits(&mut scalar, degree_resp_inv);

    let p = dim_two_ker.t1.p2;
    ec_mul(
        &mut dim_two_ker.t1.p2,
        &scalar,
        reduced_order,
        &p,
        &mut e_com_x_aux.e2,
    );
    let p = dim_two_ker.t2.p2;
    ec_mul(
        &mut dim_two_ker.t2.p2,
        &scalar,
        reduced_order,
        &p,
        &mut e_com_x_aux.e2,
    );
    let p = dim_two_ker.t1m2.p2;
    ec_mul(
        &mut dim_two_ker.t1m2.p2,
        &scalar,
        reduced_order,
        &p,
        &mut e_com_x_aux.e2,
    );

    let n = exp_diadic_val_full_resp as u32;
    let s = dim_two_ker.t1;
    double_couple_point_iter(&mut dim_two_ker.t1, n, &s, &e_com_x_aux);
    let s = dim_two_ker.t2;
    double_couple_point_iter(&mut dim_two_ker.t2, n, &s, &e_com_x_aux);
    let s = dim_two_ker.t1m2;
    double_couple_point_iter(&mut dim_two_ker.t1m2, n, &s, &e_com_x_aux);

    let mut pushed_points = [ThetaCouplePoint::default(); 3];
    copy_point(&mut pushed_points[0].p1, &domain.b1.p);
    copy_point(&mut pushed_points[1].p1, &domain.b1.q);
    copy_point(&mut pushed_points[2].p1, &domain.b1.pmq);
    ec_point_init(&mut pushed_points[0].p2);
    ec_point_init(&mut pushed_points[1].p2);
    ec_point_init(&mut pushed_points[2].p2);

    let mut codomain_product = ThetaCoupleCurve::default();

    if theta_chain_compute_and_eval_randomized(
        pow_dim2_deg_resp as u32,
        &mut e_com_x_aux,
        &dim_two_ker,
        true,
        &mut codomain_product,
        &mut pushed_points,
    ) == 0
    {
        return 0;
    }

    debug_assert!(
        test_couple_point_order_twof(&pushed_points[0], &codomain_product, reduced_order) != 0
    );

    copy_curve(&mut codomain.e1, &codomain_product.e2);
    copy_curve(&mut codomain.e2, &codomain_product.e1);

    copy_point(&mut codomain.b1.p, &pushed_points[0].p2);
    copy_point(&mut codomain.b1.q, &pushed_points[1].p2);
    copy_point(&mut codomain.b1.pmq, &pushed_points[2].p2);

    copy_point(&mut codomain.b2.p, &pushed_points[0].p1);
    copy_point(&mut codomain.b2.q, &pushed_points[1].p1);
    copy_point(&mut codomain.b2.pmq, &pushed_points[2].p1);
    1
}

fn compute_small_chain_isogeny_signature(
    e_chall_2: &mut EcCurve,
    b_chall_2: &mut EcBasis,
    resp_quat: &QuatAlgElem,
    pow_dim2_deg_resp: i32,
    length: i32,
) -> i32 {
    let mut ret = 1;
    let mut two_pow = Ibz::default();
    let mut vec_resp_two = ibz_vec_2_init();
    let mut lideal_resp_two = QuatLeftIdeal::default();

    ibz_pow(&mut two_pow, ibz_const_two(), length as u32);
    quat_lideal_create(
        &mut lideal_resp_two,
        resp_quat,
        &two_pow,
        maxord_o0(),
        quatalg_pinfty(),
    );
    id2iso_ideal_to_kernel_dlogs_even(&mut vec_resp_two, &lideal_resp_two);

    let mut points = [b_chall_2.p, b_chall_2.q, b_chall_2.pmq];

    let bcc = *b_chall_2;
    ec_dbl_iter_basis(
        b_chall_2,
        pow_dim2_deg_resp + HD_EXTRA_TORSION as i32,
        &bcc,
        e_chall_2,
    );
    debug_assert!(test_basis_order_twof(b_chall_2, e_chall_2, length) != 0);

    let mut ker = EcPoint::default();
    ec_biscalar_mul_ibz_vec(&mut ker, &vec_resp_two, length, b_chall_2, e_chall_2);
    debug_assert!(test_point_order_twof(&ker, e_chall_2, length) != 0);

    if ec_eval_small_chain(e_chall_2, &ker, length, &mut points, true) != 0 {
        ret = 0;
    }
    copy_point(&mut b_chall_2.p, &points[0]);
    copy_point(&mut b_chall_2.q, &points[1]);
    copy_point(&mut b_chall_2.pmq, &points[2]);

    ret
}

fn compute_challenge_codomain_signature(
    sig: &Signature,
    sk: &mut SecretKey,
    e_chall: &mut EcCurve,
    _e_chall_2: &EcCurve,
    b_chall_2: &mut EcBasis,
) -> i32 {
    let mut phi_chall = EcIsogEven::default();
    let bas_sk = sk.canonical_basis;

    phi_chall.curve = sk.curve;
    phi_chall.length = (TORSION_EVEN_POWER - sig.backtracking as usize) as u32;
    debug_assert!(test_basis_order_twof(&bas_sk, &sk.curve, TORSION_EVEN_POWER as i32) != 0);

    ec_ladder3pt(
        &mut phi_chall.kernel,
        &sig.chall_coeff,
        &bas_sk.p,
        &bas_sk.q,
        &bas_sk.pmq,
        &sk.curve,
    );
    debug_assert!(
        test_point_order_twof(&phi_chall.kernel, &sk.curve, TORSION_EVEN_POWER as i32) != 0
    );

    let pk = phi_chall.kernel;
    ec_dbl_iter(
        &mut phi_chall.kernel,
        sig.backtracking as i32,
        &pk,
        &mut sk.curve,
    );

    if ec_eval_even(e_chall, &phi_chall, &mut []) != 0 {
        return 0;
    }

    #[cfg(debug_assertions)]
    {
        let mut j_chall = Fp2::default();
        let mut j_codomain = Fp2::default();
        ec_j_inv(&mut j_codomain, _e_chall_2);
        ec_j_inv(&mut j_chall, e_chall);
        debug_assert!(crate::gf::fp2_is_equal(&j_chall, &j_codomain) != 0);
    }

    let mut isom = EcIsom::default();
    if ec_isomorphism(&mut isom, _e_chall_2, e_chall) != 0 {
        return 0;
    }
    ec_iso_eval(&mut b_chall_2.p, &isom);
    ec_iso_eval(&mut b_chall_2.q, &isom);
    ec_iso_eval(&mut b_chall_2.pmq, &isom);
    1
}

fn set_aux_curve_signature(sig: &mut Signature, e_aux: &mut EcCurve) {
    ec_normalize_curve(e_aux);
    fp2_copy(&mut sig.e_aux_a, &e_aux.a);
}

fn compute_and_set_basis_change_matrix(
    sig: &mut Signature,
    b_aux_2: &EcBasis,
    b_chall_2: &mut EcBasis,
    e_aux_2: &mut EcCurve,
    e_chall: &mut EcCurve,
    f: i32,
) {
    let mut mat_baux2_to_baux2_can = ibz_mat_2x2_init();
    let mut mat_bchall_can_to_bchall = ibz_mat_2x2_init();

    let mut b_can_chall = EcBasis::default();
    let mut b_aux_2_can = EcBasis::default();
    sig.hint_chall =
        ec_curve_to_basis_2f_to_hint(&mut b_can_chall, e_chall, TORSION_EVEN_POWER as i32);
    sig.hint_aux =
        ec_curve_to_basis_2f_to_hint(&mut b_aux_2_can, e_aux_2, TORSION_EVEN_POWER as i32);

    #[cfg(debug_assertions)]
    {
        debug_assert!(test_basis_order_twof(&b_aux_2_can, e_aux_2, TORSION_EVEN_POWER as i32) != 0);
        debug_assert!(test_basis_order_twof(b_aux_2, e_aux_2, f) != 0);
        let mut w0 = Fp2::default();
        weil(
            &mut w0,
            f as u32,
            &b_aux_2.p,
            &b_aux_2.q,
            &b_aux_2.pmq,
            e_aux_2,
        );
    }

    change_of_basis_matrix_tate_invert(
        &mut mat_baux2_to_baux2_can,
        &b_aux_2_can,
        b_aux_2,
        e_aux_2,
        f,
    );

    matrix_application_even_basis(b_chall_2, e_chall, &mut mat_baux2_to_baux2_can, f);

    debug_assert!(test_basis_order_twof(&b_can_chall, e_chall, TORSION_EVEN_POWER as i32) != 0);

    change_of_basis_matrix_tate(
        &mut mat_bchall_can_to_bchall,
        b_chall_2,
        &b_can_chall,
        e_chall,
        f,
    );

    for i in 0..2 {
        for j in 0..2 {
            debug_assert!(
                ibz_bitsize(&mat_bchall_can_to_bchall[i][j])
                    <= (SQISIGN_RESPONSE_LENGTH + HD_EXTRA_TORSION as usize) as i32
            );
            ibz_to_digits(
                &mut sig.mat_bchall_can_to_b_chall[i][j],
                &mat_bchall_can_to_bchall[i][j],
            );
        }
    }
}

pub fn protocols_sign(sig: &mut Signature, pk: &PublicKey, sk: &mut SecretKey, m: &[u8]) -> i32 {
    let mut ret = 0;
    let mut reduced_order: i32 = 0;
    let mut pow_dim2_deg_resp: u8;

    let mut remain = Ibz::default();
    let mut lattice_content = Ibz::default();
    let mut random_aux_norm = Ibz::default();
    let mut degree_resp_inv = Ibz::default();

    let mut resp_quat = QuatAlgElem::default();
    let mut lideal_commit = QuatLeftIdeal::default();
    let mut lideal_com_resp = QuatLeftIdeal::default();

    let mut ecom_eaux = ThetaCoupleCurveWithBasis::default();
    let mut eaux2_echall2 = ThetaCoupleCurveWithBasis::default();

    let mut e_chall = sk.curve;

    ec_curve_init(&mut ecom_eaux.e1);
    ec_curve_init(&mut ecom_eaux.e2);

    while ret == 0 {
        ret = commit(&mut ecom_eaux.e1, &mut ecom_eaux.b1, &mut lideal_commit) as i32;
        if ret == 0 {
            continue;
        }

        hash_to_challenge(&mut sig.chall_coeff, pk, &ecom_eaux.e1, m);

        {
            let mut lideal_chall_two = QuatLeftIdeal::default();
            compute_challenge_ideal_signature(&mut lideal_chall_two, sig, sk);
            compute_response_quat_element(
                &mut resp_quat,
                &mut lattice_content,
                sk,
                &lideal_chall_two,
                &lideal_commit,
            );
        }

        compute_backtracking_signature(sig, &mut resp_quat, &mut lattice_content, &mut remain);

        pow_dim2_deg_resp = compute_random_aux_norm_and_helpers(
            sig,
            &mut random_aux_norm,
            &mut degree_resp_inv,
            &mut remain,
            &lattice_content,
            &mut resp_quat,
            &mut lideal_com_resp,
            &lideal_commit,
        );

        if pow_dim2_deg_resp > 0 {
            ret = evaluate_random_aux_isogeny_signature(
                &mut ecom_eaux.e2,
                &mut ecom_eaux.b2,
                &random_aux_norm,
                &lideal_com_resp,
            );
            if ret == 0 {
                continue;
            }

            debug_assert!(
                test_basis_order_twof(&ecom_eaux.b1, &ecom_eaux.e1, TORSION_EVEN_POWER as i32) != 0
            );
            debug_assert!(
                test_basis_order_twof(&ecom_eaux.b2, &ecom_eaux.e2, TORSION_EVEN_POWER as i32) != 0
            );

            reduced_order =
                pow_dim2_deg_resp as i32 + HD_EXTRA_TORSION as i32 + sig.two_resp_length as i32;
            let bb = ecom_eaux.b1;
            ec_dbl_iter_basis(
                &mut ecom_eaux.b1,
                TORSION_EVEN_POWER as i32 - reduced_order,
                &bb,
                &mut ecom_eaux.e1,
            );
            let bb = ecom_eaux.b2;
            ec_dbl_iter_basis(
                &mut ecom_eaux.b2,
                TORSION_EVEN_POWER as i32 - reduced_order,
                &bb,
                &mut ecom_eaux.e2,
            );

            ret = compute_dim2_isogeny_challenge(
                &mut eaux2_echall2,
                &ecom_eaux,
                &degree_resp_inv,
                pow_dim2_deg_resp as i32,
                sig.two_resp_length as i32,
                reduced_order,
            );
            if ret == 0 {
                continue;
            }
        } else {
            copy_curve(&mut eaux2_echall2.e1, &ecom_eaux.e1);
            copy_curve(&mut eaux2_echall2.e2, &ecom_eaux.e1);
            reduced_order = sig.two_resp_length as i32;
            let bb = ecom_eaux.b1;
            ec_dbl_iter_basis(
                &mut eaux2_echall2.b1,
                TORSION_EVEN_POWER as i32 - reduced_order,
                &bb,
                &mut ecom_eaux.e1,
            );
            // C duplicates this line; faithful port.
            let bb = ecom_eaux.b1;
            ec_dbl_iter_basis(
                &mut eaux2_echall2.b1,
                TORSION_EVEN_POWER as i32 - reduced_order,
                &bb,
                &mut ecom_eaux.e1,
            );
            copy_basis(&mut eaux2_echall2.b2, &eaux2_echall2.b1);
        }

        if sig.two_resp_length > 0 {
            let ok = compute_small_chain_isogeny_signature(
                &mut eaux2_echall2.e2,
                &mut eaux2_echall2.b2,
                &resp_quat,
                pow_dim2_deg_resp as i32,
                sig.two_resp_length as i32,
            );
            debug_assert!(ok != 0);
        }

        let ok = compute_challenge_codomain_signature(
            sig,
            sk,
            &mut e_chall,
            &eaux2_echall2.e2,
            &mut eaux2_echall2.b2,
        );
        debug_assert!(ok != 0);
        let _ = ok;
    }

    set_aux_curve_signature(sig, &mut eaux2_echall2.e1);

    let b_aux_2 = eaux2_echall2.b1;
    compute_and_set_basis_change_matrix(
        sig,
        &b_aux_2,
        &mut eaux2_echall2.b2,
        &mut eaux2_echall2.e1,
        &mut e_chall,
        reduced_order,
    );

    ret
}

// ---------------------------------------------------------------------------
// Top-level wrappers (sqisign.c)
// ---------------------------------------------------------------------------

pub fn sqisign_keypair() -> Result<([u8; CRYPTO_PUBLICKEYBYTES], [u8; CRYPTO_SECRETKEYBYTES]), ()> {
    let mut sk = SecretKey::default();
    let mut pk = PublicKey::default();
    secret_key_init(&mut sk);
    if protocols_keygen(&mut pk, &mut sk) == 0 {
        return Err(());
    }
    let mut sk_bytes = [0u8; CRYPTO_SECRETKEYBYTES];
    let mut pk_bytes = [0u8; CRYPTO_PUBLICKEYBYTES];
    secret_key_to_bytes(&mut sk_bytes, &sk, &pk);
    public_key_to_bytes(&mut pk_bytes, &pk);
    Ok((pk_bytes, sk_bytes))
}

pub fn sqisign_sign(m: &[u8], sk_bytes: &[u8]) -> Result<[u8; CRYPTO_BYTES], ()> {
    if sk_bytes.len() != CRYPTO_SECRETKEYBYTES {
        return Err(());
    }
    let mut sk = SecretKey::default();
    let mut pk = PublicKey::default();
    secret_key_init(&mut sk);
    secret_key_from_bytes(&mut sk, &mut pk, sk_bytes);

    let mut sig = Signature::default();
    if protocols_sign(&mut sig, &pk, &mut sk, m) == 0 {
        return Err(());
    }
    let mut sig_bytes = [0u8; CRYPTO_BYTES];
    signature_to_bytes(&mut sig_bytes, &sig);
    Ok(sig_bytes)
}
