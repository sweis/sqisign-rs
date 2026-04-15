//! Finite-field arithmetic: GF(p) and GF(p²) for the active SQIsign prime.
//!
//! Two interchangeable GF(p) backends at lvl1:
//! - default: vendored saturated 4-limb Montgomery (`gf5_248_vendored`) with
//!   hand-written x86_64 BMI2/ADX inline asm for `mul`/`sqr` when the target
//!   supports it; pure-Rust intrinsics otherwise.
//! - `gf-portable` feature: in-tree modarith (unsaturated radix-2⁵¹) reference.
//!
//! lvl3/lvl5 always use the in-tree modarith reference.
//!
//! `fp2.rs` is backend-agnostic and builds GF(p²) on top of `Fp`.

#[cfg(all(
    feature = "lvl1",
    not(feature = "lvl3"),
    not(feature = "lvl5"),
    not(feature = "gf-portable")
))]
#[path = "fp_via_gf5248.rs"]
pub mod fp;

#[cfg(all(
    feature = "lvl1",
    not(feature = "lvl3"),
    not(feature = "lvl5"),
    feature = "gf-portable"
))]
#[path = "fp_p5248.rs"]
mod backend;

#[cfg(all(feature = "lvl3", not(feature = "lvl5")))]
#[path = "fp_p65376.rs"]
mod backend;

#[cfg(feature = "lvl5")]
#[path = "fp_p27500.rs"]
mod backend;

#[cfg(any(
    feature = "lvl3",
    feature = "lvl5",
    all(feature = "lvl1", feature = "gf-portable")
))]
pub mod fp;

pub mod fp2;

pub use fp::*;
pub use fp2::*;
