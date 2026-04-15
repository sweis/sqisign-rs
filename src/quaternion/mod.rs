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

pub mod dpe;
pub mod hnf;
pub mod ideal;
pub mod lat_ball;
pub mod lattice;
pub mod lll;
pub mod normeq;
pub mod rationals;

pub use algebra::*;
pub use dim2::*;
pub use dim4::*;
pub use hnf::*;
pub use ideal::*;
pub use intbig::*;
pub use integers::*;
pub use lat_ball::*;
pub use lattice::*;
pub use lll::*;
pub use normeq::*;
pub use rationals::*;
pub use types::*;
