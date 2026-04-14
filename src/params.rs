// SPDX-License-Identifier: Apache-2.0
//! Compile-time parameter selection for SQIsign security levels.

#[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
pub mod active {
    pub const SECURITY_LEVEL: usize = 1;
    pub const CRYPTO_SECRETKEYBYTES: usize = 353;
    pub const CRYPTO_PUBLICKEYBYTES: usize = 65;
    pub const CRYPTO_BYTES: usize = 148;
    pub const CRYPTO_ALGNAME: &str = "SQIsign_lvl1";
}

#[cfg(feature = "lvl3")]
pub mod active {
    pub const SECURITY_LEVEL: usize = 3;
    pub const CRYPTO_SECRETKEYBYTES: usize = 529;
    pub const CRYPTO_PUBLICKEYBYTES: usize = 97;
    pub const CRYPTO_BYTES: usize = 224;
    pub const CRYPTO_ALGNAME: &str = "SQIsign_lvl3";
}

#[cfg(feature = "lvl5")]
pub mod active {
    pub const SECURITY_LEVEL: usize = 5;
    pub const CRYPTO_SECRETKEYBYTES: usize = 701;
    pub const CRYPTO_PUBLICKEYBYTES: usize = 129;
    pub const CRYPTO_BYTES: usize = 292;
    pub const CRYPTO_ALGNAME: &str = "SQIsign_lvl5";
}

pub use active::*;
