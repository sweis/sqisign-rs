//! Finite-field arithmetic: GF(p) and GF(p²) for the active SQIsign prime.
//!
//! Two interchangeable GF(p) backends:
//! - default: in-tree modarith (unsaturated radix-2⁵¹/⁵⁵/⁵¹ Montgomery)
//! - `gf-fp2crate` feature: the `fp2` crate's saturated 4-limb Montgomery (lvl1 only)
//!
//! `fp2.rs` is backend-agnostic and builds GF(p²) on top of `Fp`.

#[cfg(all(
    not(feature = "gf-fp2crate"),
    feature = "lvl1",
    not(feature = "lvl3"),
    not(feature = "lvl5")
))]
#[path = "fp_p5248.rs"]
mod backend;

#[cfg(all(not(feature = "gf-fp2crate"), feature = "lvl3", not(feature = "lvl5")))]
#[path = "fp_p65376.rs"]
mod backend;

#[cfg(all(not(feature = "gf-fp2crate"), feature = "lvl5"))]
#[path = "fp_p27500.rs"]
mod backend;

#[cfg(not(feature = "gf-fp2crate"))]
pub mod fp;
#[cfg(feature = "gf-fp2crate")]
#[path = "fp_via_fp2crate.rs"]
pub mod fp;

pub mod fp2;

pub use fp::*;
pub use fp2::*;
