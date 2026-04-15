//! NIST AES-256 CTR-DRBG used for deterministic KAT generation.
//!
//! Faithful port of `randombytes_ctrdrbg.c` from the NIST PQC submission
//! framework. This is **not** a secure RNG for production use; it exists
//! solely to reproduce Known Answer Tests byte-for-byte.

use aes::cipher::{generic_array::GenericArray, BlockEncrypt, KeyInit};
use aes::Aes256;
use std::sync::Mutex;

/// AES-256 CTR-DRBG state. Exposed so callers (e.g. parallel tests) can hold
/// independent instances; the global [`randombytes`] wrappers use a shared
/// instance to match the C API.
#[derive(Clone, zeroize::Zeroize, zeroize::ZeroizeOnDrop)]
pub struct DrbgState {
    key: [u8; 32],
    v: [u8; 16],
    reseed_counter: u32,
}

impl core::fmt::Debug for DrbgState {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("DrbgState").finish_non_exhaustive()
    }
}

#[inline]
fn aes256_ecb(cipher: &Aes256, ctr: &[u8; 16], out: &mut [u8; 16]) {
    let mut block = GenericArray::clone_from_slice(ctr);
    cipher.encrypt_block(&mut block);
    out.copy_from_slice(&block);
}

#[inline]
fn increment_v(v: &mut [u8; 16]) {
    for j in (0..16).rev() {
        if v[j] == 0xff {
            v[j] = 0x00;
        } else {
            v[j] += 1;
            break;
        }
    }
}

/// Derive next key/V from three blocks of `cipher`, mixing in `provided_data`.
fn ctr_drbg_update(
    cipher: &Aes256,
    provided_data: Option<&[u8; 48]>,
    key: &mut [u8; 32],
    v: &mut [u8; 16],
) {
    let mut temp = [0u8; 48];
    for i in 0..3 {
        increment_v(v);
        let mut block = [0u8; 16];
        aes256_ecb(cipher, v, &mut block);
        temp[16 * i..16 * (i + 1)].copy_from_slice(&block);
    }
    if let Some(pd) = provided_data {
        for i in 0..48 {
            temp[i] ^= pd[i];
        }
    }
    key.copy_from_slice(&temp[..32]);
    v.copy_from_slice(&temp[32..48]);
}

impl DrbgState {
    /// Create a new DRBG seeded with 48 bytes of entropy and an optional
    /// personalization string.
    pub fn new(entropy_input: &[u8; 48], personalization: Option<&[u8; 48]>) -> Self {
        let mut seed_material = *entropy_input;
        if let Some(ps) = personalization {
            for i in 0..48 {
                seed_material[i] ^= ps[i];
            }
        }
        let mut key = [0u8; 32];
        let mut v = [0u8; 16];
        let cipher = Aes256::new(GenericArray::from_slice(&key));
        ctr_drbg_update(&cipher, Some(&seed_material), &mut key, &mut v);
        DrbgState {
            key,
            v,
            reseed_counter: 1,
        }
    }

    /// Fill `x` with the next bytes from this DRBG instance.
    pub fn fill(&mut self, x: &mut [u8]) {
        // Expand the key schedule once per call, not per block.
        let cipher = Aes256::new(GenericArray::from_slice(&self.key));
        let mut block = [0u8; 16];
        for chunk in x.chunks_mut(16) {
            increment_v(&mut self.v);
            aes256_ecb(&cipher, &self.v, &mut block);
            chunk.copy_from_slice(&block[..chunk.len()]);
        }
        ctr_drbg_update(&cipher, None, &mut self.key, &mut self.v);
        self.reseed_counter += 1;
    }
}

static DRBG_CTX: Mutex<Option<DrbgState>> = Mutex::new(None);

/// Initialize the global DRBG with a 48-byte entropy seed.
/// Mirrors `randombytes_init(entropy_input, NULL, 256)` in the C code.
pub fn randombytes_init(entropy_input: &[u8]) {
    let entropy: &[u8; 48] = entropy_input
        .try_into()
        .expect("randombytes_init requires 48-byte entropy input");
    *DRBG_CTX.lock().unwrap() = Some(DrbgState::new(entropy, None));
}

/// Initialize with an optional personalization string (full NIST signature).
pub fn randombytes_init_full(
    entropy_input: &[u8; 48],
    personalization: Option<&[u8; 48]>,
    _security_strength: i32,
) {
    *DRBG_CTX.lock().unwrap() = Some(DrbgState::new(entropy_input, personalization));
}

/// Snapshot the global DRBG state for deterministic benchmark replays.
/// Gated because it exposes secret DRBG state (key + counter).
#[cfg(feature = "bench-internals")]
#[doc(hidden)]
pub fn snapshot() -> DrbgState {
    DRBG_CTX
        .lock()
        .unwrap()
        .clone()
        .expect("DRBG not initialized")
}

/// Restore the global DRBG state from a snapshot.
#[cfg(feature = "bench-internals")]
#[doc(hidden)]
pub fn restore(s: &DrbgState) {
    *DRBG_CTX.lock().unwrap() = Some(s.clone());
}

/// Fill `x` with deterministic pseudo-random bytes from the global DRBG.
/// Panics if `randombytes_init` has not been called.
pub fn randombytes(x: &mut [u8]) {
    DRBG_CTX
        .lock()
        .unwrap()
        .as_mut()
        .expect("randombytes called before randombytes_init")
        .fill(x);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test vectors extracted from the C reference (`randombytes_ctrdrbg.c`)
    /// seeded with KAT seed 0. Uses a local DrbgState to avoid racing with
    /// other tests on the global.
    #[test]
    fn drbg_matches_c_reference() {
        let seed: [u8; 48] = hex::decode(
            "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7\
             056A8C266F9EF97ED08541DBD2E1FFA1",
        )
        .unwrap()
        .try_into()
        .unwrap();
        let mut drbg = DrbgState::new(&seed, None);

        let mut out1 = [0u8; 32];
        drbg.fill(&mut out1);
        assert_eq!(
            hex::encode_upper(out1),
            "7C9935A0B07694AA0C6D10E4DB6B1ADD2FD81A25CCB148032DCD739936737F2D"
        );

        let mut out2 = [0u8; 16];
        drbg.fill(&mut out2);
        assert_eq!(hex::encode_upper(out2), "8626ED79D451140800E03B59B956F821");

        let mut out3 = [0u8; 33];
        drbg.fill(&mut out3);
        assert_eq!(
            hex::encode_upper(out3),
            "EFB3B24DA2BCF2C843FF1580EF5A1C1B25B59350EDFF47D56940692F0BB1B640FB"
        );
    }

    /// The PQCgenKAT outer loop seeds itself with entropy[i]=i, then derives
    /// the per-test seeds. The first derived seed must equal KAT seed 0.
    #[test]
    fn drbg_outer_kat_loop() {
        let ent: [u8; 48] = core::array::from_fn(|i| i as u8);
        let mut drbg = DrbgState::new(&ent, None);
        let mut s0 = [0u8; 48];
        drbg.fill(&mut s0);
        assert_eq!(
            hex::encode_upper(s0),
            "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7\
             056A8C266F9EF97ED08541DBD2E1FFA1"
        );
    }

    /// Smoke test for the global wrappers. Runs in isolation since other
    /// unit tests use local instances.
    #[test]
    fn global_wrappers() {
        let ent: [u8; 48] = core::array::from_fn(|i| i as u8);
        randombytes_init(&ent);
        let mut s0 = [0u8; 48];
        randombytes(&mut s0);
        assert_eq!(&s0[..4], &hex::decode("06155023").unwrap()[..],);
    }
}
