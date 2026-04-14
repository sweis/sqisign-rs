// SPDX-License-Identifier: Apache-2.0
//! Field-size constants (`fp_constants.h`, lvl1, RADIX=64, ref GF impl).

/// Number of limbs in a base-field element.
///
/// The reference (modarith) implementation uses an unsaturated 5-limb
/// representation; the Broadwell-optimised path uses 4 saturated limbs.
pub const NWORDS_FIELD: usize = 5;

/// Number of limbs in a group-order scalar.
pub const NWORDS_ORDER: usize = 4;

/// Bit length of the field modulus rounded up to a multiple of 64.
pub const BITS: usize = 256;

/// ⌈log₂ p⌉ mod 64, used by encoding routines.
pub const LOG2P: usize = 8;
