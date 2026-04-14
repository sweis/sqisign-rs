// SPDX-License-Identifier: Apache-2.0
//! Wire-format sizes (`encoded_sizes.h`, lvl1).

pub const SECURITY_BITS: usize = 128;
pub const SQISIGN_RESPONSE_LENGTH: usize = 126;
pub const HASH_ITERATIONS: usize = 64;
pub const FP_ENCODED_BYTES: usize = 32;
pub const FP2_ENCODED_BYTES: usize = 64;
pub const EC_CURVE_ENCODED_BYTES: usize = 64;
pub const EC_POINT_ENCODED_BYTES: usize = 64;
pub const EC_BASIS_ENCODED_BYTES: usize = 192;
pub const PUBLICKEY_BYTES: usize = 65;
pub const SECRETKEY_BYTES: usize = 353;
pub const SIGNATURE_BYTES: usize = 148;
