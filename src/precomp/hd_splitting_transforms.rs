//! Precomputed splitting/normalisation transforms for the (2,2)-isogeny
//! theta-model code (`hd_splitting_transforms.c`, lvl1).

use super::Fp2;

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
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3],
];

pub const CHI_EVAL: [[i32; 4]; 4] = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]];

use crate::gf::{MINUS_ONE as N1, ONE as P1, ZERO as Z};

/// {0, 1, i, -1, -i} in the internal field representation.
pub const FP2_CONSTANTS: [Fp2; 5] = [
    Fp2 { re: Z, im: Z },
    Fp2 { re: P1, im: Z },
    Fp2 { re: Z, im: P1 },
    Fp2 { re: N1, im: Z },
    Fp2 { re: Z, im: N1 },
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
