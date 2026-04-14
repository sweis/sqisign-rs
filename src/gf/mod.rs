// SPDX-License-Identifier: Apache-2.0
//! Finite-field arithmetic: GF(p) and GF(p²) for the active SQIsign prime.

#[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
#[path = "fp_p5248.rs"]
mod backend;

#[cfg(all(feature = "lvl3", not(feature = "lvl5")))]
#[path = "fp_p65376.rs"]
mod backend;

#[cfg(feature = "lvl5")]
#[path = "fp_p27500.rs"]
mod backend;

pub mod fp;
pub mod fp2;

pub use fp::*;
pub use fp2::*;
