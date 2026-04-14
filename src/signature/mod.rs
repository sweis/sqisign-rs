// SPDX-License-Identifier: Apache-2.0
//! SQIsign key generation and signing.

use crate::params::*;

pub fn sqisign_keypair() -> Result<([u8; CRYPTO_PUBLICKEYBYTES], [u8; CRYPTO_SECRETKEYBYTES]), ()> {
    todo!("port from C signature/ref/")
}

pub fn sqisign_sign(_m: &[u8], _sk: &[u8]) -> Result<[u8; CRYPTO_BYTES], ()> {
    todo!("port from C signature/ref/")
}
