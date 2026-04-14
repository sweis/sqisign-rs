// SPDX-License-Identifier: Apache-2.0
//! Precomputed constants for SQIsign.
//!
//! Ported from `src/precomp/ref/lvl1/` in the C reference implementation.
//! Values are for RADIX=64, reference (non-broadwell) GF implementation.

// -----------------------------------------------------------------------------
// Local type aliases.
//
// These mirror the layout used by the C reference implementation and the
// `gf`/`ec` modules being ported in parallel. Reconcile by re-exporting the
// real types from `crate::gf` once that module lands; the data layout is
// identical.
// -----------------------------------------------------------------------------

// Field-element types come from `gf`; precomp data is stored in their
// internal representation.
pub use crate::gf::{Digit, Fp, Fp2};

pub mod fp_constants;
pub mod encoded_sizes;
pub mod ec_params;
pub mod e0_basis;
pub mod hd_splitting_transforms;
pub mod torsion_constants;

pub use fp_constants::*;
pub use encoded_sizes::*;
pub use ec_params::*;
pub use e0_basis::*;
pub use hd_splitting_transforms::*;
pub use torsion_constants::*;

#[cfg(feature = "sign")]
pub mod sign_data;
#[cfg(feature = "sign")]
pub use sign_data::*;
