// SPDX-License-Identifier: Apache-2.0
//! Rust implementation of SQIsign, a post-quantum digital signature scheme
//! based on isogenies of supersingular elliptic curves.
//!
//! Ported from the reference C implementation at <https://github.com/SQIsign/the-sqisign>.

#![allow(dead_code)] // TODO: remove once port is complete
#![allow(unused_imports)]

pub mod common;
pub mod mp;
pub mod gf;
pub mod ec;
pub mod precomp;
#[cfg(not(quaternion_isolation))]
pub mod hd;
#[cfg(not(quaternion_isolation))]
pub mod verification;

#[cfg(feature = "sign")]
pub mod quaternion;
#[cfg(feature = "sign")]
pub mod id2iso;
#[cfg(feature = "sign")]
pub mod signature;

#[cfg(not(quaternion_isolation))]
pub mod nistapi;

pub mod params;
