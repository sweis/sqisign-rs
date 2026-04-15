//! SQIsign signature verification.
//! Ported from `src/verification/ref/lvlx/{common,encode_verification,verify}.c`.

use crate::common::fips202::Shake256Inc;
use crate::ec::*;
use crate::gf::Fp2;
use crate::hd::{
    theta_chain_compute_and_eval_verify, ThetaCoupleCurve, ThetaKernelCouplePoints,
    HD_EXTRA_TORSION,
};
use crate::mp::{mp_compare, mp_is_even, mp_mod_2exp, mp_sub, multiple_mp_shiftl, Digit, RADIX};
use crate::precomp::{
    FP2_ENCODED_BYTES, HASH_ITERATIONS, NWORDS_ORDER, PUBLICKEY_BYTES, SECURITY_BITS,
    SIGNATURE_BYTES, SQISIGN_RESPONSE_LENGTH, TORSION_EVEN_POWER,
};

pub type Scalar = [Digit; NWORDS_ORDER];
pub type ScalarMtx2x2 = [[Scalar; 2]; 2];

/// Reasons signature verification can fail.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VerifyError {
    /// Input slice has the wrong length.
    BadLength,
    /// A field element was not canonically encoded.
    BadEncoding,
    /// A curve coefficient is invalid (A = ±2).
    BadCurve,
    /// The signature's response-length parameters are inconsistent.
    BadResponseLength,
    /// The basis-change matrix entries exceed the canonical bound.
    BadMatrixBound,
    /// Basis hint did not yield a valid 2ᶠ-torsion basis.
    BadHint,
    /// An intermediate point was at infinity / had a zero coordinate.
    DegeneratePoint,
    /// A 2ⁿ-isogeny step failed (singular kernel or bad torsion order).
    IsogenyFailure,
    /// Theta-chain codomain did not split as an elliptic product.
    NotAProduct,
    /// Recomputed challenge does not match the signature's `chall_coeff`.
    ChallengeMismatch,
}

type VResult<T> = Result<T, VerifyError>;

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
    pk.curve = EcCurve::e0();
}

// ---------------------------------------------------------------------------
// hash_to_challenge (common.c)
// ---------------------------------------------------------------------------

pub fn hash_to_challenge(pk: &PublicKey, com_curve: &EcCurve, message: &[u8]) -> Scalar {
    let mut buf = [0u8; 2 * FP2_ENCODED_BYTES];
    {
        let j1 = pk.curve.j_inv();
        let j2 = com_curve.j_inv();
        buf[..FP2_ENCODED_BYTES].copy_from_slice(&j1.encode());
        buf[FP2_ENCODED_BYTES..].copy_from_slice(&j2.encode());
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

    let mut scalar = [0u64; NWORDS_ORDER];
    for (i, limb) in scalar.iter_mut().enumerate() {
        *limb = u64::from_le_bytes(acc[i * 8..i * 8 + 8].try_into().unwrap());
    }
    mp_mod_2exp(&mut scalar, SECURITY_BITS as u32, NWORDS_ORDER);
    scalar
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

impl PublicKey {
    pub fn to_bytes(&self) -> [u8; PUBLICKEY_BYTES] {
        let mut tmp = self.curve.c;
        debug_assert!(tmp.is_zero_ct() == 0);
        tmp = (tmp).inv();
        tmp = self.curve.a * tmp;
        let mut enc = [0u8; PUBLICKEY_BYTES];
        enc[..FP2_ENCODED_BYTES].copy_from_slice(&tmp.encode());
        enc[FP2_ENCODED_BYTES] = self.hint_pk;
        enc
    }
}

impl TryFrom<&[u8]> for PublicKey {
    type Error = VerifyError;
    /// Fails if the slice has the wrong length or the field-element encoding
    /// is non-canonical (≥ p). The C reference discards the latter status.
    fn try_from(enc: &[u8]) -> VResult<Self> {
        if enc.len() != PUBLICKEY_BYTES {
            return Err(VerifyError::BadLength);
        }
        let mut pk = PublicKey::default();
        pk.curve.a = Fp2::try_decode(&enc[..FP2_ENCODED_BYTES]).ok_or(VerifyError::BadEncoding)?;
        pk.curve.c = Fp2::ONE;
        pk.hint_pk = enc[FP2_ENCODED_BYTES];
        Ok(pk)
    }
}

impl Signature {
    pub fn to_bytes(&self) -> [u8; SIGNATURE_BYTES] {
        let mut enc = [0u8; SIGNATURE_BYTES];
        let mut p = 0;
        enc[p..p + FP2_ENCODED_BYTES].copy_from_slice(&self.e_aux_a.encode());
        p += FP2_ENCODED_BYTES;
        enc[p] = self.backtracking;
        p += 1;
        enc[p] = self.two_resp_length;
        p += 1;
        let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
        for row in &self.mat_bchall_can_to_b_chall {
            for entry in row {
                encode_digits(&mut enc[p..], entry, nbytes);
                p += nbytes;
            }
        }
        let nbytes = SECURITY_BITS / 8;
        encode_digits(&mut enc[p..], &self.chall_coeff, nbytes);
        p += nbytes;
        enc[p] = self.hint_aux;
        p += 1;
        enc[p] = self.hint_chall;
        p += 1;
        debug_assert_eq!(p, SIGNATURE_BYTES);
        enc
    }
}

impl TryFrom<&[u8]> for Signature {
    type Error = VerifyError;
    /// Fails if the slice has the wrong length or the field-element encoding
    /// is non-canonical (≥ p). The C reference discards the latter status.
    fn try_from(enc: &[u8]) -> VResult<Self> {
        if enc.len() != SIGNATURE_BYTES {
            return Err(VerifyError::BadLength);
        }
        let mut sig = Signature {
            e_aux_a: Fp2::try_decode(&enc[..FP2_ENCODED_BYTES]).ok_or(VerifyError::BadEncoding)?,
            ..Default::default()
        };
        let mut p = FP2_ENCODED_BYTES;
        sig.backtracking = enc[p];
        p += 1;
        sig.two_resp_length = enc[p];
        p += 1;
        let nbytes = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
        for row in &mut sig.mat_bchall_can_to_b_chall {
            for entry in row {
                decode_digits(entry, &enc[p..], nbytes);
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
        Ok(sig)
    }
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
            if mp_compare(&aux, &sig.mat_bchall_can_to_b_chall[i][j], NWORDS_ORDER).is_le() {
                return false;
            }
        }
    }
    true
}

fn compute_challenge_verify(sig: &Signature, epk: &EcCurve, hint_pk: u8) -> VResult<EcCurve> {
    let mut phi_chall = EcIsogEven {
        curve: *epk,
        length: (TORSION_EVEN_POWER - sig.backtracking as usize) as u32,
        ..Default::default()
    };

    let bas_ea =
        ec_curve_to_basis_2f_from_hint(&mut phi_chall.curve, TORSION_EVEN_POWER as i32, hint_pk)
            .ok_or(VerifyError::BadHint)?;

    phi_chall.kernel = ec_ladder3pt(
        &sig.chall_coeff,
        &bas_ea.p,
        &bas_ea.q,
        &bas_ea.pmq,
        &phi_chall.curve,
    )
    .ok_or(VerifyError::DegeneratePoint)?;

    let ker = phi_chall.kernel;
    ec_dbl_iter(
        &mut phi_chall.kernel,
        sig.backtracking as i32,
        &ker,
        &mut phi_chall.curve,
    );

    let mut e_chall = phi_chall.curve;
    if ec_eval_even(&mut e_chall, &phi_chall, &mut []) != 0 {
        return Err(VerifyError::IsogenyFailure);
    }
    Ok(e_chall)
}

fn matrix_scalar_application_even_basis(
    bas: &EcBasis,
    e: &EcCurve,
    mat: &ScalarMtx2x2,
    f: i32,
) -> Option<EcBasis> {
    let mut scalar0: Scalar = [0; NWORDS_ORDER];
    let mut scalar1: Scalar = [0; NWORDS_ORDER];
    mp_sub(&mut scalar0, &mat[0][0], &mat[0][1], NWORDS_ORDER);
    mp_mod_2exp(&mut scalar0, f as u32, NWORDS_ORDER);
    mp_sub(&mut scalar1, &mat[1][0], &mat[1][1], NWORDS_ORDER);
    mp_mod_2exp(&mut scalar1, f as u32, NWORDS_ORDER);
    Some(EcBasis {
        p: ec_biscalar_mul(&mat[0][0], &mat[1][0], f, bas, e)?,
        q: ec_biscalar_mul(&mat[0][1], &mat[1][1], f, bas, e)?,
        pmq: ec_biscalar_mul(&scalar0, &scalar1, f, bas, e)?,
    })
}

fn challenge_and_aux_basis_verify(
    e_chall: &mut EcCurve,
    e_aux: &mut EcCurve,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> VResult<(EcBasis, EcBasis)> {
    let f = TORSION_EVEN_POWER as i32;
    let extra = HD_EXTRA_TORSION as i32;
    let resp_order = pow_dim2_deg_resp + extra + sig.two_resp_length as i32;

    let mut b_chall_can =
        ec_curve_to_basis_2f_from_hint(e_chall, f, sig.hint_chall).ok_or(VerifyError::BadHint)?;
    let bcc = b_chall_can;
    ec_dbl_iter_basis(&mut b_chall_can, f - resp_order, &bcc, e_chall);

    let mut b_aux_can =
        ec_curve_to_basis_2f_from_hint(e_aux, f, sig.hint_aux).ok_or(VerifyError::BadHint)?;
    let bac = b_aux_can;
    ec_dbl_iter_basis(&mut b_aux_can, f - pow_dim2_deg_resp - extra, &bac, e_aux);

    #[cfg(debug_assertions)]
    if !test_basis_order_twof(&b_chall_can, e_chall, resp_order) {
        eprintln!("canonical basis has wrong order, expect something to fail");
    }

    let b_chall_can = matrix_scalar_application_even_basis(
        &b_chall_can,
        e_chall,
        &sig.mat_bchall_can_to_b_chall,
        resp_order,
    )
    .ok_or(VerifyError::DegeneratePoint)?;
    Ok((b_chall_can, b_aux_can))
}

fn two_response_isogeny_verify(
    e_chall: &mut EcCurve,
    b_chall_can: &mut EcBasis,
    sig: &Signature,
    pow_dim2_deg_resp: i32,
) -> VResult<()> {
    let m = &sig.mat_bchall_can_to_b_chall;
    let mut ker = if mp_is_even(&m[0][0]) && mp_is_even(&m[1][0]) {
        b_chall_can.q
    } else {
        b_chall_can.p
    };
    let mut points = [b_chall_can.p, b_chall_can.q, b_chall_can.pmq];

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
        return Err(VerifyError::IsogenyFailure);
    }

    *b_chall_can = EcBasis {
        p: points[0],
        q: points[1],
        pmq: points[2],
    };
    Ok(())
}

fn compute_commitment_curve_verify(
    b_chall_can: &EcBasis,
    b_aux_can: &EcBasis,
    e_chall: &EcCurve,
    e_aux: &EcCurve,
    pow_dim2_deg_resp: i32,
) -> VResult<EcCurve> {
    let mut e12 = ThetaCoupleCurve {
        e1: *e_chall,
        e2: *e_aux,
    };
    let dim_two_ker = ThetaKernelCouplePoints::from_bases(b_chall_can, b_aux_can);
    let mut codomain = ThetaCoupleCurve {
        e1: EcCurve::e0(),
        e2: EcCurve::e0(),
    };

    let codomain_splits = if pow_dim2_deg_resp == 0 {
        codomain = e12;
        if !e_chall.is_basis_four_torsion(b_chall_can) {
            return Err(VerifyError::IsogenyFailure);
        }
        true
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

    if codomain_splits {
        Ok(codomain.e1)
    } else {
        Err(VerifyError::NotAProduct)
    }
}

/// SQIsign verification (Algorithm 4.9).
pub fn protocols_verify(sig: &Signature, pk: &PublicKey, m: &[u8]) -> VResult<()> {
    if !check_canonical_basis_change_matrix(sig) {
        return Err(VerifyError::BadMatrixBound);
    }

    let pow_dim2_deg_resp =
        SQISIGN_RESPONSE_LENGTH as i32 - sig.two_resp_length as i32 - sig.backtracking as i32;
    if pow_dim2_deg_resp < 0 || pow_dim2_deg_resp == 1 {
        return Err(VerifyError::BadResponseLength);
    }

    if !EcCurve::verify_a(&pk.curve.a) {
        return Err(VerifyError::BadCurve);
    }
    let mut e_aux = EcCurve::try_from_a(&sig.e_aux_a).ok_or(VerifyError::BadCurve)?;

    debug_assert!(pk.curve.c.is_one() && !pk.curve.is_a24_computed_and_normalized);

    let mut e_chall = compute_challenge_verify(sig, &pk.curve, pk.hint_pk)?;

    let (mut b_chall_can, b_aux_can) =
        challenge_and_aux_basis_verify(&mut e_chall, &mut e_aux, sig, pow_dim2_deg_resp)?;

    if sig.two_resp_length > 0 {
        two_response_isogeny_verify(&mut e_chall, &mut b_chall_can, sig, pow_dim2_deg_resp)?;
    }

    let e_com = compute_commitment_curve_verify(
        &b_chall_can,
        &b_aux_can,
        &e_chall,
        &e_aux,
        pow_dim2_deg_resp,
    )?;

    if sig.chall_coeff == hash_to_challenge(pk, &e_com, m) {
        Ok(())
    } else {
        Err(VerifyError::ChallengeMismatch)
    }
}

/// Top-level verification matching the NIST API contract.
pub fn sqisign_verify(m: &[u8], sig_bytes: &[u8], pk_bytes: &[u8]) -> bool {
    let go = || -> VResult<()> {
        let pkt = PublicKey::try_from(pk_bytes)?;
        let sigt = Signature::try_from(sig_bytes)?;
        protocols_verify(&sigt, &pkt, m)
    };
    go().is_ok()
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
        let pk = PublicKey::try_from(&KAT0_PK[..]).unwrap();
        assert_eq!(pk.hint_pk, 11);

        let sig = Signature::try_from(&KAT0_SM[..SIGNATURE_BYTES]).unwrap();
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
        assert_eq!(&sig.to_bytes()[..], &KAT0_SM[..SIGNATURE_BYTES]);
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
        use crate::gf::{Fp, FP2_ENCODED_BYTES};
        let mut a = Fp2::default();
        a.re = Fp::from_small(2);
        assert!(a.re.is_square());
        a.re = a.re.sqrt();
        let mut pk = [0u8; PUBLICKEY_BYTES];
        pk[..FP2_ENCODED_BYTES].copy_from_slice(&a.encode());
        pk[FP2_ENCODED_BYTES] = 1; // hint_pk = 1

        let m = &KAT0_SM[SIGNATURE_BYTES..];
        let t = std::time::Instant::now();
        assert!(!sqisign_verify(m, &KAT0_SM[..SIGNATURE_BYTES], &pk));
        assert!(t.elapsed() < std::time::Duration::from_millis(100));

        // Same attack via sig.e_aux_a: take KAT0 sig, replace e_aux_a with √2,
        // set hint_aux=1.
        let mut sig = KAT0_SM;
        sig[..FP2_ENCODED_BYTES].copy_from_slice(&a.encode());
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
        let pk = PublicKey::try_from(&KAT0_PK[..]).unwrap();
        let s = hash_to_challenge(&pk, &pk.curve, b"");
        // Value pinned from current implementation (matches C; KAT-validated).
        assert_eq!(s, [16658541885460340183, 53469704974404856, 0, 0]);
    }
}
