// SPDX-License-Identifier: Apache-2.0
//! Precomputed 2^e-torsion basis on the starting curve E₀ (`e0_basis.c`, lvl1).
//!
//! Values are stored in the internal (Montgomery-like) representation used by
//! the reference GF(p) implementation.

use super::{Fp, Fp2};

pub const BASIS_E0_PX: Fp2 = Fp2 {
    re: Fp([
        0x5bcab12000c08,
        0x452654b56d052,
        0x26f81b5190a0a,
        0x36cfd66a361eb,
        0x12726610d11b,
    ]),
    im: Fp([
        0x6b96065c83efc,
        0x29da1d4a82cd9,
        0x190797ab98bdf,
        0x6841aa6eeee05,
        0x1377c5431166,
    ]),
};

pub const BASIS_E0_QX: Fp2 = Fp2 {
    re: Fp([
        0x21dd55b97832f,
        0x210f2d30b26ad,
        0x680bcfcf6396,
        0x27b318ec126a7,
        0x4ffba5956012,
    ]),
    im: Fp([
        0x74590149117e3,
        0x4982edefcc606,
        0x2ae3db0cc6884,
        0x7d0384872f5ec,
        0x4fbb0fcb5a52,
    ]),
};
