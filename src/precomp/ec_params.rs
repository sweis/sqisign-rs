// SPDX-License-Identifier: Apache-2.0
//! Elliptic-curve parameters (`ec_params.h` / `ec_params.c`, lvl1).

use super::Digit;

/// Power of 2 dividing p+1.
pub const TORSION_EVEN_POWER: usize = 248;

/// (p+1) / 2^TORSION_EVEN_POWER.
pub const P_COFACTOR_FOR_2F: [Digit; 1] = [5];

/// Bit length of `P_COFACTOR_FOR_2F`.
pub const P_COFACTOR_FOR_2F_BITLENGTH: usize = 3;
