// SPDX-License-Identifier: Apache-2.0
//! Precomputed splitting/normalisation transforms for the (2,2)-isogeny
//! theta-model code (`hd_splitting_transforms.c`, lvl1).

use super::{Fp, Fp2};

/// 4×4 matrix of indices into [`FP2_CONSTANTS`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrecompBasisChangeMatrix {
    pub m: [[u8; 4]; 4],
}

// Indices into FP2_CONSTANTS.
pub const FP2_ZERO: u8 = 0;
pub const FP2_ONE: u8 = 1;
pub const FP2_I: u8 = 2;
pub const FP2_MINUS_ONE: u8 = 3;
pub const FP2_MINUS_I: u8 = 4;

pub const EVEN_INDEX: [[i32; 2]; 10] = [
    [0, 0], [0, 1], [0, 2], [0, 3], [1, 0],
    [1, 2], [2, 0], [2, 1], [3, 0], [3, 3],
];

pub const CHI_EVAL: [[i32; 4]; 4] = [
    [1, 1, 1, 1],
    [1, -1, 1, -1],
    [1, 1, -1, -1],
    [1, -1, -1, 1],
];

const FP_ZERO: Fp = Fp([0, 0, 0, 0, 0]);
const FP_ONE_INTERNAL: Fp = Fp([0x19, 0x0, 0x0, 0x0, 0x300000000000]);
const FP_MINUS_ONE_INTERNAL: Fp = Fp([
    0x7ffffffffffe6,
    0x7ffffffffffff,
    0x7ffffffffffff,
    0x7ffffffffffff,
    0x1fffffffffff,
]);

/// {0, 1, i, -1, -i} in the internal field representation.
pub const FP2_CONSTANTS: [Fp2; 5] = [
    Fp2 { re: FP_ZERO, im: FP_ZERO },
    Fp2 { re: FP_ONE_INTERNAL, im: FP_ZERO },
    Fp2 { re: FP_ZERO, im: FP_ONE_INTERNAL },
    Fp2 { re: FP_MINUS_ONE_INTERNAL, im: FP_ZERO },
    Fp2 { re: FP_ZERO, im: FP_MINUS_ONE_INTERNAL },
];

pub const SPLITTING_TRANSFORMS: [PrecompBasisChangeMatrix; 10] = [
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_I, FP2_ONE, FP2_I],
            [FP2_ONE, FP2_MINUS_I, FP2_MINUS_ONE, FP2_I],
            [FP2_ONE, FP2_I, FP2_MINUS_ONE, FP2_MINUS_I],
            [FP2_MINUS_ONE, FP2_I, FP2_MINUS_ONE, FP2_I],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
            [FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO, FP2_ZERO],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
            [FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
            [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
            [FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
            [FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
            [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
            [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
            [FP2_MINUS_ONE, FP2_ONE, FP2_ONE, FP2_MINUS_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        ],
    },
];

pub const NORMALIZATION_TRANSFORMS: [PrecompBasisChangeMatrix; 6] = [
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE],
            [FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO],
            [FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO],
            [FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE],
            [FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE],
            [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE],
            [FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE],
            [FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE],
            [FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE],
            [FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I],
            [FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I],
            [FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE],
            [FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I],
            [FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I],
            [FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE],
        ],
    },
];
