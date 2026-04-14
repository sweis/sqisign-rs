// SPDX-License-Identifier: Apache-2.0
//! Common primitives: hashing, deterministic RNG, secure memory.
//!
//! Ports `src/common/` from the C reference. Timing/benchmark helpers
//! (`tools.c`, `bench.h`) are intentionally omitted.

pub mod ctrdrbg;
pub mod fips202;
pub mod mem;

pub use ctrdrbg::{randombytes, randombytes_init};
pub use fips202::{shake256, Shake256Inc};
pub use mem::secure_clear;
