//! Rust implementation of SQIsign, a post-quantum digital signature scheme
//! based on isogenies of supersingular elliptic curves.
//!
//! Ported from the reference C implementation at <https://github.com/SQIsign/the-sqisign>.

#![cfg_attr(
    feature = "cryptobigint",
    allow(clippy::large_stack_arrays, clippy::large_stack_frames)
)]

pub mod common;
pub mod ec;
pub mod gf;
pub mod hd;
pub mod mp;
pub mod precomp;
pub mod verification;

#[cfg(feature = "sign")]
pub mod id2iso;
#[cfg(feature = "sign")]
pub mod quaternion;
#[cfg(feature = "sign")]
pub mod signature;

pub mod nistapi;

pub mod params;

#[cfg(test)]
pub(crate) mod test_util;
