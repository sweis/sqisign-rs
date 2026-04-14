// SPDX-License-Identifier: Apache-2.0
//! NIST PQC API surface for SQIsign.

use crate::params::*;

/// Verify a signed message. Returns the recovered message on success.
pub fn crypto_sign_open(sm: &[u8], pk: &[u8]) -> Result<Vec<u8>, ()> {
    if sm.len() < CRYPTO_BYTES || pk.len() != CRYPTO_PUBLICKEYBYTES {
        return Err(());
    }
    let sig = &sm[..CRYPTO_BYTES];
    let m = &sm[CRYPTO_BYTES..];
    if crate::verification::sqisign_verify(m, sig, pk) {
        Ok(m.to_vec())
    } else {
        Err(())
    }
}

#[cfg(feature = "sign")]
pub fn crypto_sign_keypair() -> Result<([u8; CRYPTO_PUBLICKEYBYTES], [u8; CRYPTO_SECRETKEYBYTES]), ()> {
    crate::signature::sqisign_keypair()
}

#[cfg(feature = "sign")]
pub fn crypto_sign(m: &[u8], sk: &[u8]) -> Result<Vec<u8>, ()> {
    let sig = crate::signature::sqisign_sign(m, sk)?;
    let mut sm = Vec::with_capacity(CRYPTO_BYTES + m.len());
    sm.extend_from_slice(&sig);
    sm.extend_from_slice(m);
    Ok(sm)
}
