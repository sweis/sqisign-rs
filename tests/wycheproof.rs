//! Wycheproof-style edge-case tests for SQIsign verification.
//!
//! These tests probe the attacker-facing surface (`crypto_sign_open`) with
//! malformed, boundary, and adversarially-crafted inputs in the spirit of
//! <https://github.com/C2SP/wycheproof>. Cases are generated programmatically
//! by mutating KAT vector 0.
//!
//! Run:      cargo test --no-default-features --features lvl1 --release wycheproof
//! Export:   WYCHEPROOF_EXPORT=1 cargo test ... wycheproof_export -- --nocapture
//!
//! Any case that fails the asserted expectation indicates either a bug in this
//! implementation or an under-specified aspect of SQIsign verification.

use sqisign_rs::nistapi::crypto_sign_open;
use sqisign_rs::params::{CRYPTO_ALGNAME, CRYPTO_BYTES, CRYPTO_PUBLICKEYBYTES};
use sqisign_rs::precomp::{
    FP2_ENCODED_BYTES, FP_ENCODED_BYTES, SECURITY_BITS, SQISIGN_RESPONSE_LENGTH,
};
use sqisign_rs::verification::Signature;

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

// ===========================================================================
// Signature byte layout (matches verification::signature_{from,to}_bytes).
// All offsets are in the SIGNATURE portion of `sm` (i.e. sm[..CRYPTO_BYTES]).
// ===========================================================================

const SIG_E_AUX_A: usize = 0; // ..FP2_ENCODED_BYTES
const SIG_BACKTRACKING: usize = FP2_ENCODED_BYTES;
const SIG_TWO_RESP_LEN: usize = SIG_BACKTRACKING + 1;
const SIG_MAT: usize = SIG_TWO_RESP_LEN + 1;
const NBYTES_MAT: usize = (SQISIGN_RESPONSE_LENGTH + 9) / 8;
const SIG_CHALL_COEFF: usize = SIG_MAT + 4 * NBYTES_MAT;
const NBYTES_CHALL: usize = SECURITY_BITS / 8;
const SIG_HINT_AUX: usize = SIG_CHALL_COEFF + NBYTES_CHALL;
const SIG_HINT_CHALL: usize = SIG_HINT_AUX + 1;

const PK_CURVE_A: usize = 0; // ..FP2_ENCODED_BYTES
const PK_HINT: usize = FP2_ENCODED_BYTES;

// HD_EXTRA_TORSION = 2 (constant across levels).
const HD_EXTRA_TORSION: usize = 2;

// ===========================================================================
// Expectation for a test case.
// ===========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Expect {
    /// Verification MUST succeed.
    Valid,
    /// Verification MUST reject.
    Invalid,
    /// Verification MUST reject and complete within a small time bound
    /// (no DoS amplification on adversarial inputs).
    InvalidFast,
    /// Either accept or reject is allowed. Used for inputs that are
    /// distinct encodings of valid signatures (malleability) where the
    /// SQIsign spec/reference does not mandate strict canonical checking.
    Acceptable,
}

#[derive(Debug, Clone)]
struct WycheproofCase {
    tcid: u32,
    category: &'static str,
    comment: String,
    flags: Vec<&'static str>,
    pk: Vec<u8>,
    sm: Vec<u8>,
    expected: Expect,
}

// ===========================================================================
// KAT loader (minimal; same format as tests/kat.rs).
// ===========================================================================

#[derive(Clone)]
struct Kat {
    msg: Vec<u8>,
    pk: Vec<u8>,
    sm: Vec<u8>,
}

fn load_kats(n: usize) -> Vec<Kat> {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("KAT");
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    p.push("PQCsignKAT_353_SQIsign_lvl1.rsp");
    #[cfg(feature = "lvl3")]
    p.push("PQCsignKAT_529_SQIsign_lvl3.rsp");
    #[cfg(feature = "lvl5")]
    p.push("PQCsignKAT_701_SQIsign_lvl5.rsp");

    let content = fs::read_to_string(p).expect("KAT file");
    let mut out = Vec::new();
    let (mut msg, mut pk) = (Vec::new(), Vec::new());
    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("msg = ") {
            msg = hex::decode(rest).unwrap();
        } else if let Some(rest) = line.strip_prefix("pk = ") {
            pk = hex::decode(rest).unwrap();
        } else if let Some(rest) = line.strip_prefix("sm = ") {
            let sm = hex::decode(rest).unwrap();
            out.push(Kat {
                msg: msg.clone(),
                pk: pk.clone(),
                sm,
            });
            if out.len() >= n {
                break;
            }
        }
    }
    out
}

// ===========================================================================
// Field-element encodings (lvl1: p = 5·2^248 − 1; level-generic via FP_ENCODED_BYTES).
// ===========================================================================

/// Little-endian encoding of the prime p for the active level.
fn prime_p_le() -> [u8; FP_ENCODED_BYTES] {
    // p = c · 2^k − 1 with (c, k) = (5, 248) / (65, 376) / (27, 500).
    // Equivalently: low (FP_ENCODED_BYTES − 1) bytes are 0xFF; top byte is c − 1.
    let mut p = [0xFFu8; FP_ENCODED_BYTES];
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    {
        p[FP_ENCODED_BYTES - 1] = 0x04;
    }
    #[cfg(feature = "lvl3")]
    {
        p[FP_ENCODED_BYTES - 1] = 0x40;
    }
    #[cfg(feature = "lvl5")]
    {
        p[FP_ENCODED_BYTES - 1] = 0x1A;
    }
    p
}

fn fp_le(small: u64) -> [u8; FP_ENCODED_BYTES] {
    let mut out = [0u8; FP_ENCODED_BYTES];
    out[..8].copy_from_slice(&small.to_le_bytes());
    out
}

fn fp_le_neg(small: u64) -> [u8; FP_ENCODED_BYTES] {
    // p − small (assumes small > 0).
    let mut out = prime_p_le();
    let mut borrow = small;
    for b in out.iter_mut() {
        let (v, br) = (*b).overflowing_sub((borrow & 0xFF) as u8);
        *b = v;
        borrow = (borrow >> 8) + u64::from(br);
    }
    out
}

/// Write an Fp2 value (re, im) into `dst[off..off+FP2_ENCODED_BYTES]`.
fn put_fp2(dst: &mut [u8], off: usize, re: &[u8; FP_ENCODED_BYTES], im: &[u8; FP_ENCODED_BYTES]) {
    dst[off..off + FP_ENCODED_BYTES].copy_from_slice(re);
    dst[off + FP_ENCODED_BYTES..off + FP2_ENCODED_BYTES].copy_from_slice(im);
}

// ===========================================================================
// Test-case generation.
// ===========================================================================

fn build_cases() -> Vec<WycheproofCase> {
    assert_eq!(
        SIG_HINT_CHALL + 1,
        CRYPTO_BYTES,
        "signature layout mismatch"
    );
    assert_eq!(PK_HINT + 1, CRYPTO_PUBLICKEYBYTES, "pk layout mismatch");

    let kats = load_kats(5);
    let base = &kats[0];
    let mut cases = Vec::new();
    let mut tcid = 0u32;
    let mut push = |category: &'static str,
                    comment: String,
                    flags: Vec<&'static str>,
                    pk: Vec<u8>,
                    sm: Vec<u8>,
                    expected: Expect| {
        tcid += 1;
        cases.push(WycheproofCase {
            tcid,
            category,
            comment,
            flags,
            pk,
            sm,
            expected,
        });
    };

    // -----------------------------------------------------------------------
    // 1. Valid baseline.
    // -----------------------------------------------------------------------
    for (i, k) in kats.iter().enumerate() {
        push(
            "Baseline",
            format!("KAT vector {i} verifies"),
            vec![],
            k.pk.clone(),
            k.sm.clone(),
            Expect::Valid,
        );
    }
    // Empty message: sm = sig || "" with a fresh sign would be valid, but the
    // KAT message is non-empty. Stripping the message makes it Invalid.
    {
        let sm = base.sm[..CRYPTO_BYTES].to_vec();
        push(
            "Baseline",
            "valid signature, message stripped (empty msg)".into(),
            vec!["ModifiedMessage"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }

    // -----------------------------------------------------------------------
    // 2. Length manipulation.
    // -----------------------------------------------------------------------
    for &len in &[0usize, 1, CRYPTO_BYTES - 1] {
        push(
            "Length",
            format!("sm truncated to {len} bytes (< CRYPTO_BYTES)"),
            vec!["Truncated"],
            base.pk.clone(),
            base.sm[..len.min(base.sm.len())].to_vec(),
            Expect::InvalidFast,
        );
    }
    for &len in &[
        0usize,
        1,
        CRYPTO_PUBLICKEYBYTES - 1,
        CRYPTO_PUBLICKEYBYTES + 1,
    ] {
        let pk = if len <= base.pk.len() {
            base.pk[..len].to_vec()
        } else {
            let mut v = base.pk.clone();
            v.resize(len, 0);
            v
        };
        push(
            "Length",
            format!("pk length = {len} (≠ CRYPTO_PUBLICKEYBYTES)"),
            vec!["WrongLength"],
            pk,
            base.sm.clone(),
            Expect::InvalidFast,
        );
    }
    {
        // sm with one byte appended: message now has a trailing 0x00 → different msg.
        let mut sm = base.sm.clone();
        sm.push(0);
        push(
            "Length",
            "sm with one trailing zero byte appended (message extended)".into(),
            vec!["ModifiedMessage"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }

    // -----------------------------------------------------------------------
    // 3. Single-bit-flip sweep over the signature bytes.
    //
    // FINDING (matches C reference): some high bits of the basis-change-matrix
    // entries are malleable. The exact set is value- and level-dependent: a bit
    // at index `b` of entry `mat[i][j]` is malleable when both
    //   (a) flipping it leaves the entry < 2^(RL+2-bt) (canonical bound), and
    //   (b) the downstream isogeny chain produces the same commitment curve,
    //       which depends on point orders and is not a simple closed form.
    // We therefore mark the entire mat region as Acceptable (either outcome OK)
    // and characterise the actual malleable set in `wycheproof_malleability`.
    // Confirmed against the C reference (`tools/c_malleability_check.c`).
    // -----------------------------------------------------------------------
    for byte in 0..CRYPTO_BYTES {
        for bit in 0..8 {
            let mut sm = base.sm.clone();
            sm[byte] ^= 1 << bit;
            let field = sig_field_name(byte);
            let in_mat = (SIG_MAT..SIG_CHALL_COEFF).contains(&byte);
            let (expected, flags) = if in_mat {
                (
                    Expect::Acceptable,
                    vec!["BitFlip", field, "SignatureMalleabilityMatrix"],
                )
            } else {
                (Expect::Invalid, vec!["BitFlip", field])
            };
            push(
                "BitFlip",
                format!("flip sig byte {byte} bit {bit} ({field})"),
                flags,
                base.pk.clone(),
                sm,
                expected,
            );
        }
    }
    // Single-bit-flip sweep over the public key bytes.
    for byte in 0..CRYPTO_PUBLICKEYBYTES {
        for bit in 0..8 {
            let mut pk = base.pk.clone();
            pk[byte] ^= 1 << bit;
            let field = if byte < FP2_ENCODED_BYTES {
                "pk.curve_A"
            } else {
                "pk.hint"
            };
            push(
                "BitFlip",
                format!("flip pk byte {byte} bit {bit} ({field})"),
                vec!["BitFlip", field],
                pk,
                base.sm.clone(),
                Expect::Invalid,
            );
        }
    }
    // Single-bit-flip in the message portion (a few positions).
    for &byte in &[CRYPTO_BYTES, CRYPTO_BYTES + 1, base.sm.len() - 1] {
        let mut sm = base.sm.clone();
        sm[byte] ^= 0x01;
        push(
            "BitFlip",
            format!("flip message byte {} bit 0", byte - CRYPTO_BYTES),
            vec!["BitFlip", "Message"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }

    // -----------------------------------------------------------------------
    // 4. Field-element edge cases for sig.E_aux_A and pk.curve_A.
    // -----------------------------------------------------------------------
    let zero = fp_le(0);
    let one = fp_le(1);
    let two = fp_le(2);
    let neg_two = fp_le_neg(2);
    let p_minus_1 = fp_le_neg(1);
    let p = prime_p_le();
    let all_ff = [0xFFu8; FP_ENCODED_BYTES];

    let fe_cases: &[(
        &str,
        [u8; FP_ENCODED_BYTES],
        [u8; FP_ENCODED_BYTES],
        &[&'static str],
    )] = &[
        ("A = 0 (singular curve)", zero, zero, &["SingularCurve"]),
        ("A = 1", one, zero, &["EdgeCaseA"]),
        ("A = 2 (singular curve)", two, zero, &["SingularCurve"]),
        ("A = -2 (singular curve)", neg_two, zero, &["SingularCurve"]),
        ("A = p-1 (= -1)", p_minus_1, zero, &["EdgeCaseA"]),
        (
            "A = p (non-canonical encoding of 0)",
            p,
            zero,
            &["NonCanonicalEncoding"],
        ),
        (
            "A.re = all-0xFF (non-canonical, > p)",
            all_ff,
            zero,
            &["NonCanonicalEncoding"],
        ),
        (
            "A.im = p (non-canonical)",
            zero,
            p,
            &["NonCanonicalEncoding"],
        ),
        ("A = i (pure imaginary)", zero, one, &["EdgeCaseA"]),
    ];
    for (label, re, im, flags) in fe_cases {
        // sig.E_aux_A
        let mut sm = base.sm.clone();
        put_fp2(&mut sm, SIG_E_AUX_A, re, im);
        push(
            "FieldElement",
            format!("sig.E_aux_A: {label}"),
            flags.to_vec(),
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
        // pk.curve_A
        let mut pk = base.pk.clone();
        put_fp2(&mut pk, PK_CURVE_A, re, im);
        push(
            "FieldElement",
            format!("pk.curve_A: {label}"),
            flags.to_vec(),
            pk,
            base.sm.clone(),
            Expect::Invalid,
        );
    }
    // Non-canonical re-encoding of the *valid* pk.A: encode A.re as A.re + p.
    // fp_decode returns 0 (non-canonical) and zeros the output, so this should
    // be rejected. Tests strict canonical-input handling.
    {
        let mut pk = base.pk.clone();
        // Add p to the re component (mod 2^(8*FP_ENCODED_BYTES)).
        let p_bytes = prime_p_le();
        let mut carry = 0u16;
        for i in 0..FP_ENCODED_BYTES {
            let s = pk[PK_CURVE_A + i] as u16 + p_bytes[i] as u16 + carry;
            pk[PK_CURVE_A + i] = s as u8;
            carry = s >> 8;
        }
        push(
            "FieldElement",
            "pk.curve_A.re += p (non-canonical encoding of valid A)".into(),
            vec!["NonCanonicalEncoding", "Malleability"],
            pk,
            base.sm.clone(),
            Expect::Invalid,
        );
    }

    // -----------------------------------------------------------------------
    // 5. Signature-structure edge cases.
    // -----------------------------------------------------------------------
    // backtracking values (DoS / boundary).
    for &bt in &[1u8, 2, 126, 127, 128, 200, 255] {
        let mut sm = base.sm.clone();
        sm[SIG_BACKTRACKING] = bt;
        let exp = if bt as usize >= SQISIGN_RESPONSE_LENGTH + HD_EXTRA_TORSION {
            Expect::InvalidFast
        } else {
            Expect::Invalid
        };
        push(
            "Structure",
            format!("sig.backtracking = {bt}"),
            vec!["Backtracking"],
            base.pk.clone(),
            sm,
            exp,
        );
    }
    // two_resp_length values (skip the original value — that's a no-op).
    let orig_trl = base.sm[SIG_TWO_RESP_LEN];
    for &trl in &[0u8, 2, 126, 127, 128, 255] {
        if trl == orig_trl {
            continue;
        }
        let mut sm = base.sm.clone();
        sm[SIG_TWO_RESP_LEN] = trl;
        push(
            "Structure",
            format!("sig.two_resp_length = {trl}"),
            vec!["TwoRespLength"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }
    // pow_dim2_deg_resp == 1 special rejection (bt + trl = SQISIGN_RESPONSE_LENGTH − 1).
    {
        let mut sm = base.sm.clone();
        sm[SIG_BACKTRACKING] = 0;
        sm[SIG_TWO_RESP_LEN] = (SQISIGN_RESPONSE_LENGTH - 1) as u8;
        push(
            "Structure",
            "pow_dim2_deg_resp = 1 (explicit reject)".into(),
            vec!["Dim2Degree1"],
            base.pk.clone(),
            sm,
            Expect::InvalidFast,
        );
    }
    // hint bytes (skip values equal to the original — no-op mutations).
    for &h in &[0u8, 1, 0x7F, 0x80, 0xFF] {
        if h != base.sm[SIG_HINT_AUX] {
            let mut sm = base.sm.clone();
            sm[SIG_HINT_AUX] = h;
            push(
                "Structure",
                format!("sig.hint_aux = {h:#04x}"),
                vec!["Hint"],
                base.pk.clone(),
                sm,
                Expect::Invalid,
            );
        }
        if h != base.sm[SIG_HINT_CHALL] {
            let mut sm = base.sm.clone();
            sm[SIG_HINT_CHALL] = h;
            push(
                "Structure",
                format!("sig.hint_chall = {h:#04x}"),
                vec!["Hint"],
                base.pk.clone(),
                sm,
                Expect::Invalid,
            );
        }
        if h != base.pk[PK_HINT] {
            let mut pk = base.pk.clone();
            pk[PK_HINT] = h;
            push(
                "Structure",
                format!("pk.hint = {h:#04x}"),
                vec!["Hint"],
                pk,
                base.sm.clone(),
                Expect::Invalid,
            );
        }
    }
    // All-zero matrix.
    {
        let mut sm = base.sm.clone();
        for b in &mut sm[SIG_MAT..SIG_MAT + 4 * NBYTES_MAT] {
            *b = 0;
        }
        push(
            "Structure",
            "all-zero basis-change matrix".into(),
            vec!["ZeroMatrix"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }
    // Matrix entry at the canonical bound 2^(response_length+2−backtracking).
    // For KAT 0 (bt=0): bound = 2^(SQISIGN_RESPONSE_LENGTH+2). Set mat[0][0] = bound.
    {
        let bits = SQISIGN_RESPONSE_LENGTH + HD_EXTRA_TORSION; // bt=0 in KAT0
        let mut sm = base.sm.clone();
        for b in &mut sm[SIG_MAT..SIG_MAT + NBYTES_MAT] {
            *b = 0;
        }
        // Only representable if bits < 8*NBYTES_MAT; otherwise the encoding
        // width equals the bound and the bound check is vacuous (so the
        // rejection happens later in the isogeny chain, not fast).
        let exp = if bits < 8 * NBYTES_MAT {
            sm[SIG_MAT + bits / 8] = 1 << (bits % 8);
            Expect::InvalidFast
        } else {
            for b in &mut sm[SIG_MAT..SIG_MAT + NBYTES_MAT] {
                *b = 0xFF;
            }
            Expect::Invalid
        };
        push(
            "Structure",
            "mat[0][0] = 2^(response_length+2) (≥ canonical bound)".into(),
            vec!["MatrixBound"],
            base.pk.clone(),
            sm,
            exp,
        );
    }
    // chall_coeff = 0.
    {
        let mut sm = base.sm.clone();
        for b in &mut sm[SIG_CHALL_COEFF..SIG_CHALL_COEFF + NBYTES_CHALL] {
            *b = 0;
        }
        push(
            "Structure",
            "sig.chall_coeff = 0".into(),
            vec!["ZeroChallenge"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }
    // All-zero signature / all-zero pk. A=0 is the supersingular curve E0, so
    // ec_curve_verify_a does not early-reject; verification proceeds into the
    // isogeny chain and rejects later.
    {
        let mut sm = base.sm.clone();
        for b in &mut sm[..CRYPTO_BYTES] {
            *b = 0;
        }
        push(
            "Structure",
            "all-zero signature".into(),
            vec!["AllZero"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
        push(
            "Structure",
            "all-zero pk".into(),
            vec!["AllZero"],
            vec![0u8; CRYPTO_PUBLICKEYBYTES],
            base.sm.clone(),
            Expect::Invalid,
        );
    }

    // -----------------------------------------------------------------------
    // 6. Cross-vector mismatches.
    // -----------------------------------------------------------------------
    for j in 1..kats.len() {
        push(
            "CrossKey",
            format!("sig from KAT 0 with pk from KAT {j}"),
            vec!["WrongPublicKey"],
            kats[j].pk.clone(),
            base.sm.clone(),
            Expect::Invalid,
        );
        // Replace just the message (sig stays from KAT0, msg from KATj).
        let mut sm = base.sm[..CRYPTO_BYTES].to_vec();
        sm.extend_from_slice(&kats[j].msg);
        push(
            "CrossKey",
            format!("sig from KAT 0 with msg from KAT {j}"),
            vec!["WrongMessage"],
            base.pk.clone(),
            sm,
            Expect::Invalid,
        );
    }

    cases
}

fn sig_field_name(byte: usize) -> &'static str {
    if byte < SIG_BACKTRACKING {
        "E_aux_A"
    } else if byte == SIG_BACKTRACKING {
        "backtracking"
    } else if byte == SIG_TWO_RESP_LEN {
        "two_resp_length"
    } else if byte < SIG_CHALL_COEFF {
        "mat"
    } else if byte < SIG_HINT_AUX {
        "chall_coeff"
    } else if byte == SIG_HINT_AUX {
        "hint_aux"
    } else {
        "hint_chall"
    }
}

// ===========================================================================
// Runner.
// ===========================================================================

/// Upper bound for `InvalidFast` rejections. Generous to avoid CI flakiness;
/// the cases this guards are O(1) checks that should reject in microseconds.
/// The reference DoS case took ~18 ms in release before the fix.
const FAST_REJECT_MS: u128 = 5;

fn run_case(c: &WycheproofCase) -> Result<(), String> {
    let t0 = Instant::now();
    let result = crypto_sign_open(&c.sm, &c.pk);
    let elapsed = t0.elapsed();
    match c.expected {
        Expect::Valid => {
            let m = result.map_err(|()| "expected Valid, got reject".to_string())?;
            let expected_msg = &c.sm[CRYPTO_BYTES.min(c.sm.len())..];
            if m != expected_msg {
                return Err("recovered message mismatch".into());
            }
        }
        Expect::Invalid | Expect::InvalidFast => {
            if result.is_ok() {
                return Err("expected reject, signature was ACCEPTED".into());
            }
            if c.expected == Expect::InvalidFast && elapsed.as_millis() > FAST_REJECT_MS {
                return Err(format!(
                    "expected fast reject, took {elapsed:?} (> {FAST_REJECT_MS} ms)"
                ));
            }
        }
        Expect::Acceptable => {
            // Either outcome is permitted. Surface what actually happened so a
            // behaviour change is visible in CI logs without failing the test.
            return Ok(());
        }
    }
    Ok(())
}

#[test]
fn wycheproof_layout_sanity() {
    // Compile-time-ish check that our offset table matches the implementation.
    assert_eq!(SIG_HINT_CHALL + 1, CRYPTO_BYTES);
    assert_eq!(PK_HINT + 1, CRYPTO_PUBLICKEYBYTES);
}

#[test]
fn wycheproof_malleability() {
    // 7a. Re-encoding a valid signature must be byte-identical.
    let kats = load_kats(10);
    for (i, k) in kats.iter().enumerate() {
        let sig = Signature::try_from(&k.sm[..CRYPTO_BYTES]).unwrap();
        assert_eq!(
            &sig.to_bytes()[..],
            &k.sm[..CRYPTO_BYTES],
            "signature encoding is not canonical (KAT {i})"
        );
    }

    // 7b. Characterise basis-change-matrix malleability empirically.
    //
    // FINDING (matches C reference, confirmed via tools/c_malleability_check.c):
    // bit `SQISIGN_RESPONSE_LENGTH + 1 − backtracking` of every matrix entry is
    // malleable for EVERY honest signature, giving ≥ 2^4 valid encodings each.
    // This bit is the top bit consumed by `xdblmul` (kbits = RL+2−bt) but does
    // not affect the resulting point because the basis has order exactly 2^kbits.
    // The canonical-bound check uses < 2^(RL+2−bt) so this bit always passes.
    // Additional value-dependent lower bits may also be malleable.
    eprintln!(
        "matrix malleability ({CRYPTO_ALGNAME}, RL={SQISIGN_RESPONSE_LENGTH}, NBYTES_MAT={NBYTES_MAT}):"
    );
    for (i, k) in kats.iter().enumerate() {
        let bt = k.sm[SIG_BACKTRACKING];
        let trl = k.sm[SIG_TWO_RESP_LEN];
        let target_bit = SQISIGN_RESPONSE_LENGTH + 1 - bt as usize;
        assert!(target_bit < 8 * NBYTES_MAT);
        // Universal: flipping target_bit in each entry must still verify.
        for entry in 0..4 {
            let mut sm = k.sm.clone();
            sm[SIG_MAT + entry * NBYTES_MAT + target_bit / 8] ^= 1 << (target_bit % 8);
            assert!(
                crypto_sign_open(&sm, &k.pk).is_ok(),
                "KAT {i} bt={bt}: bit {target_bit} of mat entry {entry} not malleable (finding regressed?)"
            );
        }
        // The bit one above (RL+2−bt) must be rejected by the canonical bound,
        // when representable in the encoding.
        if target_bit + 1 < 8 * NBYTES_MAT {
            let mut sm = k.sm.clone();
            sm[SIG_MAT + (target_bit + 1) / 8] ^= 1 << ((target_bit + 1) % 8);
            assert!(
                crypto_sign_open(&sm, &k.pk).is_err(),
                "KAT {i}: bit {} not rejected by canonical bound",
                target_bit + 1
            );
        }
        eprintln!("  KAT {i}: bt={bt} trl={trl} -> bit {target_bit} malleable in all 4 entries");
    }

    // 7c. For KAT 0, also report any additional malleable bits in the top byte
    // of each entry (informational; value-dependent).
    let k = &kats[0];
    let bt = k.sm[SIG_BACKTRACKING];
    let mut extra: Vec<(usize, usize)> = Vec::new();
    for entry in 0..4 {
        let top = SIG_MAT + (entry + 1) * NBYTES_MAT - 1;
        for bit in 0..8 {
            let abs = 8 * (NBYTES_MAT - 1) + bit;
            if abs == SQISIGN_RESPONSE_LENGTH + 1 - bt as usize {
                continue;
            }
            let mut sm = k.sm.clone();
            sm[top] ^= 1 << bit;
            if crypto_sign_open(&sm, &k.pk).is_ok() {
                extra.push((entry, abs));
            }
        }
    }
    eprintln!("  KAT 0: additional value-dependent malleable bits in top byte: {extra:?}");
}

#[test]
fn wycheproof_generated() {
    let cases = build_cases();
    let mut failures = Vec::new();
    let mut by_category: std::collections::BTreeMap<&str, (usize, usize)> = Default::default();
    let mut acceptable_accepted = 0usize;
    let mut acceptable_total = 0usize;
    for c in &cases {
        let entry = by_category.entry(c.category).or_default();
        entry.0 += 1;
        if c.expected == Expect::Acceptable {
            acceptable_total += 1;
            if crypto_sign_open(&c.sm, &c.pk).is_ok() {
                acceptable_accepted += 1;
            }
        }
        if let Err(msg) = run_case(c) {
            entry.1 += 1;
            failures.push(format!(
                "  tcId {} [{}] '{}': {}",
                c.tcid, c.category, c.comment, msg
            ));
        }
    }
    println!("--- Wycheproof-style summary ({} cases) ---", cases.len());
    for (cat, (n, f)) in &by_category {
        println!("  {cat:14} {n:5} cases, {f} failures");
    }
    println!(
        "  acceptable: {acceptable_accepted}/{acceptable_total} were accepted by this implementation"
    );
    if !failures.is_empty() {
        println!("--- FAILURES ---");
        for f in &failures {
            println!("{f}");
        }
        panic!(
            "{} of {} wycheproof cases failed",
            failures.len(),
            cases.len()
        );
    }
}

// ===========================================================================
// JSON export (Wycheproof-ish schema). Gated on env var so normal `cargo test`
// runs don't write files.
// ===========================================================================

#[test]
fn wycheproof_export() {
    if std::env::var("WYCHEPROOF_EXPORT").is_err() {
        eprintln!("(WYCHEPROOF_EXPORT not set; skipping JSON export)");
        return;
    }
    let cases = build_cases();
    let kats = load_kats(1);
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(format!("wycheproof_{}.json", CRYPTO_ALGNAME.to_lowercase()));

    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"algorithm\": \"{CRYPTO_ALGNAME}\",\n"));
    out.push_str("  \"generatorVersion\": \"sqisign-rs/0.1\",\n");
    out.push_str(&format!("  \"numberOfTests\": {},\n", cases.len()));
    out.push_str("  \"notes\": {\n");
    out.push_str("    \"NonCanonicalEncoding\": \"Field element encoded with value >= p\",\n");
    out.push_str("    \"SingularCurve\": \"Montgomery A in {0, +-2} gives a singular curve\",\n");
    out.push_str("    \"BitFlip\": \"Single-bit mutation of an otherwise-valid input\",\n");
    out.push_str("    \"Malleability\": \"Alternative encoding of a valid input\",\n");
    out.push_str("    \"SignatureMalleabilityMatrixHighBit\": \"When backtracking==0, the high bit of each basis-change-matrix entry is unconstrained: flipping it yields a distinct accepted signature. Confirmed against the C reference implementation; not a Rust-port bug.\"\n");
    out.push_str("  },\n");
    out.push_str("  \"testGroups\": [\n    {\n");
    out.push_str(&format!(
        "      \"publicKey\": \"{}\",\n",
        hex::encode(&kats[0].pk)
    ));
    out.push_str(&format!(
        "      \"publicKeySize\": {CRYPTO_PUBLICKEYBYTES},\n"
    ));
    out.push_str(&format!("      \"signatureSize\": {CRYPTO_BYTES},\n"));
    out.push_str("      \"tests\": [\n");
    for (idx, c) in cases.iter().enumerate() {
        let result = match c.expected {
            Expect::Valid => "valid",
            Expect::Invalid | Expect::InvalidFast => "invalid",
            Expect::Acceptable => "acceptable",
        };
        let (sig_hex, msg_hex) = if c.sm.len() >= CRYPTO_BYTES {
            (
                hex::encode(&c.sm[..CRYPTO_BYTES]),
                hex::encode(&c.sm[CRYPTO_BYTES..]),
            )
        } else {
            (hex::encode(&c.sm), String::new())
        };
        let flags_json = c
            .flags
            .iter()
            .map(|f| format!("\"{f}\""))
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str("        {\n");
        out.push_str(&format!("          \"tcId\": {},\n", c.tcid));
        out.push_str(&format!(
            "          \"comment\": \"{}\",\n",
            json_escape(&c.comment)
        ));
        out.push_str(&format!("          \"flags\": [{flags_json}],\n"));
        out.push_str(&format!("          \"pk\": \"{}\",\n", hex::encode(&c.pk)));
        out.push_str(&format!("          \"msg\": \"{msg_hex}\",\n"));
        out.push_str(&format!("          \"sig\": \"{sig_hex}\",\n"));
        out.push_str(&format!("          \"result\": \"{result}\"\n"));
        out.push_str(if idx + 1 == cases.len() {
            "        }\n"
        } else {
            "        },\n"
        });
    }
    out.push_str("      ]\n    }\n  ]\n}\n");
    fs::write(&path, out).expect("write JSON");
    println!("wrote {} ({} cases)", path.display(), cases.len());
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
