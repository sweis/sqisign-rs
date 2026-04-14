// SPDX-License-Identifier: Apache-2.0
//! Quaternion algebra arithmetic over big integers (`rug`-backed).
//!
//! Ported from `the-sqisign/src/quaternion/ref/generic/`. The function-level
//! API mirrors C (out-param first) so the rest of the signing path can be
//! ported mechanically.

pub mod intbig;
pub mod types;

pub mod algebra;
pub mod dim2;
pub mod dim4;
pub mod integers;

pub use algebra::*;
pub use dim2::*;
pub use dim4::*;
pub use intbig::*;
pub use integers::*;
pub use types::*;
