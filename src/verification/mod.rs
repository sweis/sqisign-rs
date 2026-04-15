// SPDX-License-Identifier: Apache-2.0
//! SQIsign signature verification.
//! Ported from `src/verification/ref/lvlx/{common,encode_verification,verify}.c`.

use crate::common::fips202::Shake256Inc;
use crate::ec::*;
use crate::gf::{
    fp2_decode, fp2_encode, fp2_inv, fp2_is_one, fp2_is_zero, fp2_mul, fp2_set_one, Fp2,
};
use crate::hd::{
    copy_bases_to_kernel, theta_chain_compute_and_eval_verify, ThetaCoupleCurve,
    ThetaKernelCouplePoints, HD_EXTRA_TORSION,
};
use crate::mp::{mp_compare, mp_is_even, mp_mod_2exp, mp_sub, multiple_mp_shiftl, Digit, RADIX};
use crate::precomp::{
    FP2_ENCODED_BYTES, HASH_ITERATIONS, NWORDS_ORDER, PUBLICKEY_BYTES, SECURITY_BITS,
    SIGNATURE_BYTES, SQISIGN_RESPONSE_LENGTH, TORSION_EVEN_POWER,
};

pub type Scalar = [Digit; NWORDS_ORDER];
pub type ScalarMtx2x2 = [[Scalar; 2]; 2];

#[derive(Clone, Copy, Default, Debug)]
pub struct Signature {
    pub e_aux_a: Fp2,
    pub backtracking: u8,
    pub two_resp_length: u8,
    pub mat_bchall_can_to_b_chall: ScalarMtx2x2,
    pub chall_coeff: Scalar,
    pub hint_aux: u8,
    pub hint_chall: u8,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct PublicKey {
    pub curve: EcCurve,
    pub hint_pk: u8,
}

pub fn public_key_init(pk: &mut PublicKey) {
    ec_curve_init(&mut pk.curve);
}

// ---------------------------------------------------------------------------
// hash_to_challenge (common.c)
// ---------------------------------------------------------------------------

pub fn hash_to_challenge(scalar: &mut Scalar, pk: &PublicKey, com_curve: &EcCurve, message: &[u8]) {
    let mut buf = [0u8; 2 * FP2_ENCODED_BYTES];
    {
        let mut j1 = Fp2::default();
        let mut j2 = Fp2::default();
        ec_j_inv(&mut j1, &pk.curve);
        ec_j_inv(&mut j2, com_curve);
        fp2_encode(&mut buf[..FP2_ENCODED_BYTES], &j1);
        fp2_encode(&mut buf[FP2_ENCODED_BYTES..], &j2);
    }

    let mut hash_bytes = (2 * SECURITY_BITS).div_ceil(8);
    let mut limbs = hash_bytes.div_ceil(8);
    let mut bits = (2 * SECURITY_BITS) % RADIX as usize;
    let mut mask = (!0u64) >> ((RADIX as usize - bits) % RADIX as usize);

    // Reuse the scalar buffer as bytes (LE) for the iteration loop.
    // C does this via punning; we go via a stack buffer of u64s mirrored to bytes.
    let mut acc = [0u8; NWORDS_ORDER * 8];

    let mut ctx = Shake256Inc::new();
    ctx.absorb(&buf);
    ctx.absorb(message);
    ctx.finalize();
    ctx.squeeze(&mut acc[..hash_bytes]);
    mask_top_limb(&mut acc, limbs, mask);

    for _ in 2..HASH_ITERATIONS {
        let mut ctx = Shake256Inc::new();
        ctx.absorb(&acc[..hash_bytes]);
        ctx.finalize();
        ctx.squeeze(&mut acc[..hash_bytes]);
        mask_top_limb(&mut acc, limbs, mask);
    }

    let mut ctx = Shake256Inc::new();
    ctx.absorb(&acc[..hash_bytes]);
    ctx.finalize();

    hash_bytes = (TORSION_EVEN_POWER - SQISIGN_RESPONSE_LENGTH).div_ceil(8);
    limbs = hash_bytes.div_ceil(8);
    bits = (TORSION_EVEN_POWER - SQISIGN_RESPONSE_LENGTH) % RADIX as usize;
    mask = (!0u64) >> ((RADIX as usize - bits) % RADIX as usize);

    acc.fill(0);
    ctx.squeeze(&mut acc[..hash_bytes]);
    mask_top_limb(&mut acc, limbs, mask);

    for (i, limb) in scalar.iter_mut().enumerate() {
        *limb = u64::from_le_bytes(acc[i * 8..i * 8 + 8].try_into().unwrap());
    }
    mp_mod_2exp(scalar, SECURITY_BITS as u32, NWORDS_ORDER);
}

#[inline]
fn mask_top_limb(acc: &mut [u8], limbs: usize, mask: u64) {
    let off = (limbs - 1) * 8;
    let mut top = u64::from_le_bytes(acc[off..off + 8].try_into().unwrap());
    top &= mask;
    acc[off..off + 8].copy_from_slice(&top.to_le_bytes());
}

// ---------------------------------------------------------------------------
// encode_verification.c
// ---------------------------------------------------------------------------

use crate::mp::{decode_digits, encode_digits};

pub fn public_key_to_bytes(enc: &mut [u8], pk: &PublicKey) {
    debug_assert_eq!(enc.len(), PUBLICKEY_BYTES);
    let mut tmp = pk.curve.c;
    debug_assert!(fp2_is_zero(&tmp) == 0);
    fp2_inv(&mut tmp);
    let s = tmp;
    fp2_mul(&mut tmp, &pk.curve.a, &s);
    fp2_encode(&mut enc[..FP2_ENCODED_BYTES], &tmp);
    enc[FP2_ENCODED_BYTES] = pk.hint_pk;
}

/// Decode a public key. Returns `false` if the field-element encoding is
/// non-canonical (limb value ≥ p); the C reference discards this status.
pub fn public_key_from_bytes(pk: &mut PublicKey, enc: &[u8]) -> bool {
    debug_assert_eq!(enc.len(), PUBLICKEY_BYTES);
    *pk = PublicKey::default();
    let ok = fp2_decode(&mut pk.curve.a, &enc[..FP2_ENCODED_BYTES]);
    fp2_set_one(&mut pk.curve.c);
    pk.hint_pk = enc[FP2_ENCODED_BYTES];
    ok == 0xFFFF_FFFF
}

pub fn signature_to_bytes(enc: &mut [u8], sig: &Signature) {
    debug_assert_eq!(enc.len(), SIGNATURE_BYTES);
    let mut p = 0;
    fp2_encode(&mut enc[p..p + FP2_ENCODED_BYTES], &sig.e_aux_a);
    p += FP2_ENCODED_BYTES;
    enc[p] = sig.backtracking;
    p += 1;
    enc[p] = sig.two_resp_length;
    p += 1;
    let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
    for i in 0..2 {
        for j in 0..2 {
            encode_digits(&mut enc[p..], &sig.mat_bchall_can_to_b_chall[i][j], nbytes);
            p += nbytes;
        }
    }
    let nbytes = SECURITY_BITS / 8;
    encode_digits(&mut enc[p..], &sig.chall_coeff, nbytes);
    p += nbytes;
    enc[p] = sig.hint_aux;
    p += 1;
    enc[p] = sig.hint_chall;
    p += 1;
    debug_assert_eq!(p, SIGNATURE_BYTES);
}

/// Decode a signature. Returns `false` if the field-element encoding is
/// non-canonical (limb value ≥ p); the C reference discards this status.
pub fn signature_from_bytes(sig: &mut Signature, enc: &[u8]) -> bool {
    debug_assert_eq!(enc.len(), SIGNATURE_BYTES);
    let mut p = 0;
    let ok = fp2_decode(&mut sig.e_aux_a, &enc[..FP2_ENCODED_BYTES]);
    p += FP2_ENCODED_BYTES;
    sig.backtracking = enc[p];
    p += 1;
    sig.two_resp_length = enc[p];
    p += 1;
    let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
    for i in 0..2 {
        for j in 0..2 {
            decode_digits(&mut sig.mat_bchall_can_to_b_chall[i][j], &enc[p..], nbytes);
            p += nbytes;
        }
    }
    let nbytes = SECURITY_BITS / 8;
    decode_digits(&mut sig.chall_coeff, &enc[p..], nbytes);
    p += nbytes;
    sig.hint_aux = enc[p];
    p += 1;
    sig.hint_chall = enc[p];
    p += 1;
    debug_assert_eq!(p, SIGNATURE_BYTES);
    ok == 0xFFFF_FFFF
}

// ---------------------------------------------------------------------------
// verify.c
// ---------------------------------------------------------------------------

fn check_canonical_basis_change_matrix(sig: &Signature) -> bool {
    let mut aux: Scalar = [0; NWORDS_ORDER];
    aux[0] = 1;
    let shift = SQISIGN_RESPONSE_LENGTH as i32 + HD_EXTRA_TORSION as i32 - sig.backtracking as i32;
    // Hardening vs C: a malicious `backtracking ≥ SQISIGN_RESPONSE_LENGTH+2` makes
    // `shift ≤ 0`. C (and a faithful port) then either hits UB (`>> 64` at shift=0)
    // or wraps to ~4e9 and burns ~18ms in `multiple_mp_shiftl`. Either way the
    // signature is rejected later by the `pow_dim2_deg_resp < 0` check, so bailing
    // here is equivalent and avoids the panic/amplification.
    if shift <= 0 {
        return false;
    }
    multiple_mp_shiftl(&mut aux, shift as u32, NWORDS_ORDER);
    for i in 0..2 {
        for j in 0..2 {
            if mp_compare(&aux, &sig.mat_bchall_can_to_b_chall[i][j], NWORDS_ORDER) <= 0 {
                return false;
            }
        }
    }
    true
}

fn compute_challenge_verify(
    e_chall: &mut EcCurve,
    sig: &Signature,
    epk: &EcCurve,
    hint_pk: u8,
) -> bool {
    let mut bas_ea = EcBasis::default();
    let mut phi_chall = EcIsogEven::default();
    copy_curve(&mut phi_chall.curve, epk);
    phi_chall.length = (TORSION_EVEN_POWER - sig.backtracking as usize) as u32;

    if ec_curve_to_basis_2f_from_hint(
        &mut bas_ea,
        &mut phi_chall.curve,
        TORSION_EVEN_POWER as i32,
        hint_pk,
    ) == 0
    {
        return false;
    }

    if ec_ladder3pt(
        &mut phi_chall.kernel,
        &sig.chall_coeff,
        &bas_ea.p,
        &bas_ea.q,
        &bas_ea.pmq,
        &phi_chall.curve,
    ) == 0
    {
        return false;
    }

    let ker = phi_chall.kernel;
    ec_dbl_iter(
        &mut phi_chall.kernel,
        sig.backtracking as i32,
        &ker,
        &mut phi_chall.curve,
    );

    copy_curve(e_chall, &phi_chall.curve);
    if ec_eval_even(e_chall, &phi_chall, &mut []) != 0 {
        return false;
    }
    true
}

fn matrix_scalar_application_even_basis(
    bas: &mut EcBasis,
    e: &EcCurve,
    mat: &ScalarMtx2x2,
    f: i32,
) -> bool {
    let mut scalar0: Scalar = [0; NWORDS_ORDER];
    let mut scalar1: Scalar = [0; NWORDS_ORDER];
    let tmp_bas = *bas;

    if ec_biscalar_mul(&mut bas.p, &mat[0][0], &mat[1][0], f, &tmp_bas, e) == 0 {
        return false;
    }
    if ec_biscalar_mul(&mut bas.q, &mat[0][1], &mat[1][1], f, &tmp_bas, e) == 0 {
        return false;
    }
    mp_sub(&mut scalar0, &mat[0][0], &mat[0][1], NWORDS_ORDER);
    mp_mod_2exp(&mut scalar0, f as u32, NWORDS_ORDER);
    mp_sub(&mut scalar1, &mat[1][0], &mat[1][1], NWORDS_ORDER);
    mp_mod_2exp(&mut scalar1, f as u32, NWORDS_ORDER);
    ec_biscalar_mul(&mut bas.pmq, &scalar0, &scalar1, f, &tmp_bas, e) != 0
}

fn challenge_and_aux_basis_verify(
    b_chall_can: &mut EcBasis,
    b_aux_can: &mut EcBasis,
    e_chall: &mut EcCurve,
    e_aux: &mut EcCurve,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> bool {
    if ec_curve_to_basis_2f_from_hint(
        b_chall_can,
        e_chall,
        TORSION_EVEN_POWER as i32,
        sig.hint_chall,
    ) == 0
    {
        return false;
    }
    let bcc = *b_chall_can;
    ec_dbl_iter_basis(
        b_chall_can,
        TORSION_EVEN_POWER as i32
            - pow_dim2_deg_resp
            - HD_EXTRA_TORSION as i32
            - sig.two_resp_length as i32,
        &bcc,
        e_chall,
    );

    if ec_curve_to_basis_2f_from_hint(b_aux_can, e_aux, TORSION_EVEN_POWER as i32, sig.hint_aux)
        == 0
    {
        return false;
    }
    let bac = *b_aux_can;
    ec_dbl_iter_basis(
        b_aux_can,
        TORSION_EVEN_POWER as i32 - pow_dim2_deg_resp - HD_EXTRA_TORSION as i32,
        &bac,
        e_aux,
    );

    #[cfg(debug_assertions)]
    if test_basis_order_twof(
        b_chall_can,
        e_chall,
        HD_EXTRA_TORSION as i32 + pow_dim2_deg_resp + sig.two_resp_length as i32,
    ) == 0
    {
        eprintln!("canonical basis has wrong order, expect something to fail");
    }

    matrix_scalar_application_even_basis(
        b_chall_can,
        e_chall,
        &sig.mat_bchall_can_to_b_chall,
        pow_dim2_deg_resp + HD_EXTRA_TORSION as i32 + sig.two_resp_length as i32,
    )
}

fn two_response_isogeny_verify(
    e_chall: &mut EcCurve,
    b_chall_can: &mut EcBasis,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> bool {
    let mut ker = EcPoint::default();
    let mut points = [EcPoint::default(); 3];

    if mp_is_even(&sig.mat_bchall_can_to_b_chall[0][0])
        && mp_is_even(&sig.mat_bchall_can_to_b_chall[1][0])
    {
        copy_point(&mut ker, &b_chall_can.q);
    } else {
        copy_point(&mut ker, &b_chall_can.p);
    }
    copy_point(&mut points[0], &b_chall_can.p);
    copy_point(&mut points[1], &b_chall_can.q);
    copy_point(&mut points[2], &b_chall_can.pmq);

    let k = ker;
    ec_dbl_iter(
        &mut ker,
        pow_dim2_deg_resp + HD_EXTRA_TORSION as i32,
        &k,
        e_chall,
    );

    if ec_eval_small_chain(
        e_chall,
        &ker,
        sig.two_resp_length as i32,
        &mut points,
        false,
    ) != 0
    {
        return false;
    }

    copy_point(&mut b_chall_can.p, &points[0]);
    copy_point(&mut b_chall_can.q, &points[1]);
    copy_point(&mut b_chall_can.pmq, &points[2]);
    true
}

fn compute_commitment_curve_verify(
    e_com: &mut EcCurve,
    b_chall_can: &EcBasis,
    b_aux_can: &EcBasis,
    e_chall: &EcCurve,
    e_aux: &EcCurve,
    pow_dim2_deg_resp: i32,
) -> bool {
    let mut e12 = ThetaCoupleCurve::default();
    copy_curve(&mut e12.e1, e_chall);
    copy_curve(&mut e12.e2, e_aux);

    let mut dim_two_ker = ThetaKernelCouplePoints::default();
    copy_bases_to_kernel(&mut dim_two_ker, b_chall_can, b_aux_can);

    let mut codomain = ThetaCoupleCurve::default();
    ec_curve_init(&mut codomain.e1);
    ec_curve_init(&mut codomain.e2);

    let codomain_splits = if pow_dim2_deg_resp == 0 {
        copy_curve(&mut codomain.e1, &e12.e1);
        copy_curve(&mut codomain.e2, &e12.e2);
        if ec_is_basis_four_torsion(b_chall_can, e_chall) == 0 {
            return false;
        }
        1
    } else {
        theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp as u32,
            &mut e12,
            &dim_two_ker,
            true,
            &mut codomain,
            &mut [],
        )
    };

    copy_curve(e_com, &codomain.e1);
    codomain_splits != 0
}

pub fn protocols_verify(sig: &Signature, pk: &PublicKey, m: &[u8]) -> bool {
    if !check_canonical_basis_change_matrix(sig) {
        return false;
    }

    let pow_dim2_deg_resp =
        SQISIGN_RESPONSE_LENGTH as i32 - sig.two_resp_length as i32 - sig.backtracking as i32;
    if pow_dim2_deg_resp < 0 || pow_dim2_deg_resp == 1 {
        return false;
    }

    if ec_curve_verify_a(&pk.curve.a) == 0 {
        return false;
    }

    let mut e_aux = EcCurve::default();
    if ec_curve_init_from_a(&mut e_aux, &sig.e_aux_a) == 0 {
        return false;
    }

    debug_assert!(
        fp2_is_one(&pk.curve.c) == 0xFFFF_FFFF && !pk.curve.is_a24_computed_and_normalized
    );

    let mut e_chall = EcCurve::default();
    if !compute_challenge_verify(&mut e_chall, sig, &pk.curve, pk.hint_pk) {
        return false;
    }

    let mut b_chall_can = EcBasis::default();
    let mut b_aux_can = EcBasis::default();
    if !challenge_and_aux_basis_verify(
        &mut b_chall_can,
        &mut b_aux_can,
        &mut e_chall,
        &mut e_aux,
        sig,
        pow_dim2_deg_resp,
    ) {
        return false;
    }

    if sig.two_resp_length > 0
        && !two_response_isogeny_verify(&mut e_chall, &mut b_chall_can, sig, pow_dim2_deg_resp)
    {
        return false;
    }

    let mut e_com = EcCurve::default();
    if !compute_commitment_curve_verify(
        &mut e_com,
        &b_chall_can,
        &b_aux_can,
        &e_chall,
        &e_aux,
        pow_dim2_deg_resp,
    ) {
        return false;
    }

    let mut chk_chall: Scalar = [0; NWORDS_ORDER];
    hash_to_challenge(&mut chk_chall, pk, &e_com, m);

    mp_compare(&sig.chall_coeff, &chk_chall, NWORDS_ORDER) == 0
}

/// Top-level verification matching the NIST API contract.
pub fn sqisign_verify(m: &[u8], sig_bytes: &[u8], pk_bytes: &[u8]) -> bool {
    if sig_bytes.len() != SIGNATURE_BYTES || pk_bytes.len() != PUBLICKEY_BYTES {
        return false;
    }
    let mut pkt = PublicKey::default();
    let mut sigt = Signature::default();
    if !public_key_from_bytes(&mut pkt, pk_bytes) {
        return false;
    }
    if !signature_from_bytes(&mut sigt, sig_bytes) {
        return false;
    }
    protocols_verify(&sigt, &pkt, m)
}

#[cfg(all(test, feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
mod tests {
    use super::*;
    use hex_literal::hex;

    // Golden values from tools/dump_verify_intermediates.c on KAT 0.
    const KAT0_PK: [u8; PUBLICKEY_BYTES] = hex!(
        "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202\
         9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B"
    );
    const KAT0_SM: [u8; 181] = hex!(
        "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04\
         55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01\
         0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04\
         DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A448402674\
         71D86EFF3447018ADB0A6551EE8322AB30010202D81C4D8D734FCBFBEADE3D3F\
         8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8"
    );

    #[test]
    fn decode_kat0_pk_sig() {
        let mut pk = PublicKey::default();
        public_key_from_bytes(&mut pk, &KAT0_PK);
        assert_eq!(pk.hint_pk, 11);

        let mut sig = Signature::default();
        signature_from_bytes(&mut sig, &KAT0_SM[..SIGNATURE_BYTES]);
        assert_eq!(sig.backtracking, 0);
        assert_eq!(sig.two_resp_length, 1);
        assert_eq!(sig.hint_aux, 2);
        assert_eq!(sig.hint_chall, 2);
        assert_eq!(
            sig.mat_bchall_can_to_b_chall[0][0],
            [0x57f61f67cab574d9, 0x144aa8a562b4d808, 0, 0]
        );
        assert_eq!(
            sig.chall_coeff,
            [0x0adb8a014734ff6e, 0x0130ab2283ee5165, 0, 0]
        );

        // Roundtrip.
        let mut roundtrip = [0u8; SIGNATURE_BYTES];
        signature_to_bytes(&mut roundtrip, &sig);
        assert_eq!(&roundtrip[..], &KAT0_SM[..SIGNATURE_BYTES]);
    }

    #[test]
    fn reject_wrong_order_kernel() {
        // Fuzzer-found regression: a signature whose basis-change matrix
        // produces kernel points without the expected 2³-torsion order
        // previously hit a debug_assert in `gluing_compute`. Must reject
        // gracefully (also in debug builds), not panic. The C reference
        // crashes here in debug mode (`#ifndef NDEBUG assert(...)`).
        let pk: [u8; PUBLICKEY_BYTES] = hex!(
            "1584B4356B7B4E181EDCF26987F87BE2AC278594EFE19F4A19B2B64E35D93000\
             93FB5A2203F31BA91E3D544B3D84451DC05C95006D6D5643E461CC68F749B60102"
        );
        let sig: [u8; SIGNATURE_BYTES] = hex!(
            "E0B82F57B5E76BCCE7CD8EC810AF242BE46BF0289726C9E962E1D1B21CDCEE02\
             69019AD486AE5C386201EFEB356409EDC6D83B94DB5443BD7A2BD83415623701\
             00014E35D9300093FB5A2203F31BA91E3D544B3D84451DC05C95006D6D5643E4\
             61CC68F749B60102E0B82F57B5E76BCCE7CDD26E57E315332ADA8EDDB914D0A7\
             A3470D1EC18F38E6DF5EBE2D0FE37E1D6B030B0B"
        );
        assert!(!sqisign_verify(b"", &sig, &pk));
    }

    #[test]
    fn reject_adversarial_backtracking() {
        // Regression test for the DoS/panic at backtracking ≥ 128 (negative
        // shift wraps to ~4e9 in `multiple_mp_shiftl`; see PORTING.md). All
        // these must reject quickly without panic, including in debug mode.
        let m = &KAT0_SM[SIGNATURE_BYTES..];
        for bt in [127, 128, 129, 200, 255] {
            let mut sm = KAT0_SM;
            sm[64] = bt;
            let t = std::time::Instant::now();
            assert!(!sqisign_verify(m, &sm[..SIGNATURE_BYTES], &KAT0_PK));
            assert!(
                t.elapsed().as_millis() < 5,
                "bt={bt} took {:?}",
                t.elapsed()
            );
        }
    }

    #[test]
    fn verify_kat0() {
        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let sig = &KAT0_SM[..SIGNATURE_BYTES];
        assert!(sqisign_verify(m, sig, &KAT0_PK));
    }

    /// Parse (pk, sm) for the given KAT index from the embedded .rsp file.
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    fn kat_vector(idx: usize) -> (Vec<u8>, Vec<u8>) {
        const RSP: &str = include_str!("../../KAT/PQCsignKAT_353_SQIsign_lvl1.rsp");
        let mut count = usize::MAX;
        let mut pk = vec![];
        for line in RSP.lines() {
            if let Some(v) = line.strip_prefix("count = ") {
                count = v.parse().unwrap();
            } else if count == idx {
                if let Some(v) = line.strip_prefix("pk = ") {
                    pk = (0..v.len() / 2)
                        .map(|i| u8::from_str_radix(&v[2 * i..2 * i + 2], 16).unwrap())
                        .collect();
                } else if let Some(v) = line.strip_prefix("sm = ") {
                    let sm = (0..v.len() / 2)
                        .map(|i| u8::from_str_radix(&v[2 * i..2 * i + 2], 16).unwrap())
                        .collect();
                    return (pk, sm);
                }
            }
        }
        panic!("KAT {idx} not found");
    }

    /// Verify a small set of KAT vectors covering diverse signature shapes
    /// (backtracking>0, two_resp_length=0, both-even matrix parity), so that
    /// branch arithmetic in protocols_verify and the isogeny chain are
    /// exercised by lib tests, not only by the integration suite.
    #[test]
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    fn verify_diverse_kats() {
        // 0: bt=0 trl=1 (baseline); 1: bt=1; 5: m00,m10 both even;
        // 7: trl=0; 9: bt=2; 2: trl=6 (long two-response chain).
        for idx in [0, 1, 2, 5, 7, 9] {
            let (pk, sm) = kat_vector(idx);
            assert!(
                sqisign_verify(&sm[SIGNATURE_BYTES..], &sm[..SIGNATURE_BYTES], &pk),
                "KAT{idx} should verify"
            );
            // Tampered signature byte (E_aux_A) is rejected.
            let mut bad = sm.clone();
            bad[10] ^= 1;
            assert!(
                !sqisign_verify(&bad[SIGNATURE_BYTES..], &bad[..SIGNATURE_BYTES], &pk),
                "KAT{idx} tamper should reject"
            );
        }
    }

    #[test]
    fn reject_malformed_inputs() {
        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let sig = &KAT0_SM[..SIGNATURE_BYTES];

        // Wrong sizes.
        assert!(!sqisign_verify(m, &sig[..SIGNATURE_BYTES - 1], &KAT0_PK));
        assert!(!sqisign_verify(m, sig, &KAT0_PK[..PUBLICKEY_BYTES - 1]));

        // Tampered message.
        let mut bad_m = m.to_vec();
        bad_m[0] ^= 1;
        assert!(!sqisign_verify(&bad_m, sig, &KAT0_PK));

        // Non-canonical basis change matrix: with backtracking=10 the bound is
        // 2^118; set bit 119 of mat[0][0] so it exceeds the bound and must be
        // rejected by check_canonical_basis_change_matrix.
        let mut bad = KAT0_SM;
        bad[64] = 10; // backtracking
        let mat00_byte_119 = 64 + 2 + 14;
        bad[mat00_byte_119] = 0x80;
        assert!(!sqisign_verify(m, &bad[..SIGNATURE_BYTES], &KAT0_PK));
        // ...and with the matrix in range, the same backtracking=10 still
        // rejects (computation differs), proving the matrix check itself fired.
        let mut bad2 = KAT0_SM;
        bad2[64] = 10;
        for i in 0..4 {
            let off = 64 + 2 + i * ((SQISIGN_RESPONSE_LENGTH + 9) / 8) + 14;
            bad2[off] &= 0x3F;
            bad2[off + 1] = 0;
        }
        assert!(!sqisign_verify(m, &bad2[..SIGNATURE_BYTES], &KAT0_PK));

        // E_aux_A = 2 (invalid Montgomery coefficient).
        let mut bad = KAT0_SM;
        bad[..FP2_ENCODED_BYTES].fill(0);
        bad[0] = 2;
        assert!(!sqisign_verify(m, &bad[..SIGNATURE_BYTES], &KAT0_PK));

        // pk curve A = 2.
        let mut bad_pk = KAT0_PK;
        bad_pk[..FP2_ENCODED_BYTES].fill(0);
        bad_pk[0] = 2;
        assert!(!sqisign_verify(m, sig, &bad_pk));

        // two_resp_length such that pow_dim2_deg_resp = 1 (rejected).
        let mut bad = KAT0_SM;
        bad[64] = 0; // backtracking
        bad[65] = (SQISIGN_RESPONSE_LENGTH - 1) as u8; // two_resp_length
        assert!(!sqisign_verify(m, &bad[..SIGNATURE_BYTES], &KAT0_PK));

        // two_resp_length=0 takes the no-two-response branch.
        // KAT0 has two_resp_length=1; setting it to 0 changes the computation
        // and must reject (since the matrix/chall_coeff no longer match).
        let mut bad = KAT0_SM;
        bad[65] = 0;
        assert!(!sqisign_verify(m, &bad[..SIGNATURE_BYTES], &KAT0_PK));
    }

    #[test]
    fn reject_adversarial_nqr_dos() {
        // Regression for the find_nqr_factor infinite-loop DoS: with a public-key
        // curve A satisfying Re(A²)=2 and hint_pk=1 (hint_a=true, hint_p=0), the
        // entangled-basis NQR search never terminates because Im(ib·A²-(1+ib)²)=0
        // for every b. The C reference loops unbounded; we cap the search and
        // reject. A = √2 ∈ Fp works since p ≡ 7 (mod 8).
        use crate::gf::{fp_is_square, fp_set_small, fp_sqrt, FP2_ENCODED_BYTES};
        let mut a = Fp2::default();
        fp_set_small(&mut a.re, 2);
        assert!(fp_is_square(&a.re) != 0);
        fp_sqrt(&mut a.re);
        let mut pk = [0u8; PUBLICKEY_BYTES];
        fp2_encode(&mut pk[..FP2_ENCODED_BYTES], &a);
        pk[FP2_ENCODED_BYTES] = 1; // hint_pk = 1

        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let t = std::time::Instant::now();
        assert!(!sqisign_verify(m, &KAT0_SM[..SIGNATURE_BYTES], &pk));
        assert!(t.elapsed() < std::time::Duration::from_millis(100));

        // Same attack via sig.e_aux_a: take KAT0 sig, replace e_aux_a with √2,
        // set hint_aux=1.
        let mut sig = KAT0_SM;
        fp2_encode(&mut sig[..FP2_ENCODED_BYTES], &a);
        sig[SIGNATURE_BYTES - 2] = 1; // hint_aux
        let t = std::time::Instant::now();
        assert!(!sqisign_verify(m, &sig[..SIGNATURE_BYTES], &KAT0_PK));
        assert!(t.elapsed() < std::time::Duration::from_millis(100));
    }

    #[test]
    fn reject_noncanonical_encoding() {
        // PK or sig field-element bytes ≥ p must be rejected at decode time
        // rather than silently mapping to A=0 (j=1728).
        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let mut bad_pk = KAT0_PK;
        bad_pk[..32].fill(0xFF);
        assert!(!sqisign_verify(m, &KAT0_SM[..SIGNATURE_BYTES], &bad_pk));

        let mut bad_sig = KAT0_SM;
        bad_sig[..32].fill(0xFF);
        assert!(!sqisign_verify(m, &bad_sig[..SIGNATURE_BYTES], &KAT0_PK));
    }

    #[test]
    fn reject_zero_matrix_gracefully() {
        // An all-zero basis-change matrix collapses the basis to O. This is
        // not rejected at decode time (honest signatures can have det≡0 mod 2),
        // but the lift_basis guard and theta-chain order checks must reject it
        // gracefully rather than corrupt the curve or hang.
        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let mut bad = KAT0_SM;
        for b in bad[66..66 + 4 * 17].iter_mut() {
            *b = 0;
        }
        assert!(!sqisign_verify(m, &bad[..SIGNATURE_BYTES], &KAT0_PK));
    }

    #[test]
    fn hash_to_challenge_stable() {
        // Golden test for the SHAKE256 iteration loop: hash of KAT0 (pk, com=pk)
        // with empty message must be a specific value. Any change to mask/loop
        // arithmetic breaks this.
        let mut pk = PublicKey::default();
        public_key_from_bytes(&mut pk, &KAT0_PK);
        let mut s: Scalar = [0; NWORDS_ORDER];
        hash_to_challenge(&mut s, &pk, &pk.curve, b"");
        // Value pinned from current implementation (matches C; KAT-validated).
        assert_eq!(s, [16658541885460340183, 53469704974404856, 0, 0]);
    }
}
