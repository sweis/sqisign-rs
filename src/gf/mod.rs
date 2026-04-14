// SPDX-License-Identifier: Apache-2.0
//! Finite-field arithmetic: GF(p) and GF(p²) for the active SQIsign prime.

#[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
pub mod fp;
#[cfg(any(feature = "lvl3", feature = "lvl5"))]
compile_error!("gf: lvl3/lvl5 primes not yet ported");

pub mod fp2;

pub use fp::*;
pub use fp2::*;
