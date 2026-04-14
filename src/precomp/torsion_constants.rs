// SPDX-License-Identifier: Apache-2.0
//! Big-integer torsion constants (`torsion_constants.c`, lvl1).
//!
//! In C these are `ibz_t` (GMP `mpz_t`) values. They are represented here as
//! little-endian `u64` limb arrays so they can be loaded into whatever bignum
//! type the `quaternion` module settles on (`rug::Integer`, etc.) without
//! pulling a GMP dependency into the verification path.

pub const TORSION_2POWER_BYTES: usize = 32;

/// 2^SECURITY_BITS = 2^128.
pub const TWO_TO_SECURITY_BITS: [u64; 3] = [0x0, 0x0, 0x1];

/// 2^TORSION_EVEN_POWER = 2^248.
pub const TORSION_PLUS_2POWER: [u64; 4] = [0x0, 0x0, 0x0, 0x100000000000000];

/// Secret-isogeny degree bound = 2^512 + 0x4b.
pub const SEC_DEGREE: [u64; 9] = [0x4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1];

/// Commitment-isogeny degree bound = 2^512 + 0x4b.
pub const COM_DEGREE: [u64; 9] = [0x4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1];
