// SPDX-License-Identifier: Apache-2.0
//! Precomputed constants for SQIsign.
//!
//! Ported from `src/precomp/ref/lvl{1,3,5}/` in the C reference implementation.
//! Values are for RADIX=64, reference (non-broadwell) GF implementation.

pub use crate::gf::{Digit, Fp, Fp2};

pub mod hd_splitting_transforms;
pub mod level;

pub use hd_splitting_transforms::*;
pub use level::*;

#[cfg(all(
    feature = "sign",
    feature = "lvl1",
    not(feature = "lvl3"),
    not(feature = "lvl5")
))]
#[path = "sign_data_lvl1.rs"]
pub mod sign_data;
#[cfg(all(feature = "sign", feature = "lvl3", not(feature = "lvl5")))]
#[path = "sign_data_lvl3.rs"]
pub mod sign_data;
#[cfg(all(feature = "sign", feature = "lvl5"))]
#[path = "sign_data_lvl5.rs"]
pub mod sign_data;
#[cfg(feature = "sign")]
pub use sign_data::*;
