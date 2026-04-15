//! Per-level precomputed constants. The active set is selected via cargo
//! feature; values are taken from `src/precomp/ref/lvl{1,3,5}/` in the
//! C reference (RADIX=64, non-broadwell GF impl).

#![allow(clippy::unreadable_literal)]

use super::{Digit, Fp, Fp2};

// Field-size constants are re-exported from the gf backend so they cannot drift.
pub use crate::gf::{NWORDS_FIELD, NWORDS_ORDER};
pub const BITS: usize = crate::gf::BITS as usize;
pub const LOG2P: usize = crate::gf::LOG2P as usize;
pub const FP_ENCODED_BYTES: usize = crate::gf::FP_ENCODED_BYTES;
pub const FP2_ENCODED_BYTES: usize = 2 * FP_ENCODED_BYTES;
pub const EC_CURVE_ENCODED_BYTES: usize = FP2_ENCODED_BYTES;
pub const EC_POINT_ENCODED_BYTES: usize = FP2_ENCODED_BYTES;
pub const EC_BASIS_ENCODED_BYTES: usize = 3 * FP2_ENCODED_BYTES;

// ===========================================================================
#[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
mod sel {
    use super::*;

    // encoded_sizes.h
    pub const SECURITY_BITS: usize = 128;
    pub const SQISIGN_RESPONSE_LENGTH: usize = 126;
    pub const HASH_ITERATIONS: usize = 64;
    pub const PUBLICKEY_BYTES: usize = 65;
    pub const SECRETKEY_BYTES: usize = 353;
    pub const SIGNATURE_BYTES: usize = 148;

    // ec_params.h
    pub const TORSION_EVEN_POWER: usize = 248;
    pub const P_COFACTOR_FOR_2F: [Digit; 1] = [5];
    pub const P_COFACTOR_FOR_2F_BITLENGTH: usize = 3;

    // torsion_constants.h
    pub const TORSION_2POWER_BYTES: usize = 32;
    pub const TWO_TO_SECURITY_BITS: &[u64] = &[0x0, 0x0, 0x1];
    pub const TORSION_PLUS_2POWER: &[u64] = &[0x0, 0x0, 0x0, 0x100000000000000];
    pub const SEC_DEGREE: &[u64] = &[0x4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1];
    pub const COM_DEGREE: &[u64] = &[0x4b, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1];

    // e0_basis.c
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
}

// ===========================================================================
#[cfg(all(feature = "lvl3", not(feature = "lvl5")))]
mod sel {
    use super::*;

    pub const SECURITY_BITS: usize = 192;
    pub const SQISIGN_RESPONSE_LENGTH: usize = 192;
    pub const HASH_ITERATIONS: usize = 256;
    pub const PUBLICKEY_BYTES: usize = 97;
    pub const SECRETKEY_BYTES: usize = 529;
    pub const SIGNATURE_BYTES: usize = 224;

    pub const TORSION_EVEN_POWER: usize = 376;
    pub const P_COFACTOR_FOR_2F: [Digit; 1] = [65];
    pub const P_COFACTOR_FOR_2F_BITLENGTH: usize = 7;

    pub const TORSION_2POWER_BYTES: usize = 48;
    pub const TWO_TO_SECURITY_BITS: &[u64] = &[0x0, 0x0, 0x0, 0x1];
    pub const TORSION_PLUS_2POWER: &[u64] = &[0x0, 0x0, 0x0, 0x0, 0x0, 0x100000000000000];
    pub const SEC_DEGREE: &[u64] = &[
        0xb7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1,
    ];
    pub const COM_DEGREE: &[u64] = &[
        0xb7, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1,
    ];

    pub const BASIS_E0_PX: Fp2 = Fp2 {
        re: Fp([
            0x94635b7b34b8c,
            0x431475975ec8c7,
            0x380f3b6b0f3d6c,
            0x2e90ddd88ba021,
            0x5eb0a59679b654,
            0x347706dc01cb41,
            0xb7765ed4a44a5,
        ]),
        im: Fp([
            0x412c2c3df0cc54,
            0x2338803450b7d0,
            0x206883ec0e5d2f,
            0x407a7e72205c5d,
            0x187f5a00661d99,
            0x5905b6352b7e4d,
            0x3032c0ad99418,
        ]),
    };
    pub const BASIS_E0_QX: Fp2 = Fp2 {
        re: Fp([
            0x7333ee4f4818b7,
            0x29c73aefc7681b,
            0x3db742e2128546,
            0x3f8774b65cc12a,
            0x332cf22a3425e2,
            0x4a219e343591d2,
            0x6d1dfdb6ea8ff,
        ]),
        im: Fp([
            0x1fdc82b11838c5,
            0x681f359137f9af,
            0x3eb05affc54924,
            0x509e310ef21e09,
            0x5a97b9d957fd56,
            0x6e7c043e0db389,
            0x4fbc3aab7429d,
        ]),
    };
}

// ===========================================================================
#[cfg(feature = "lvl5")]
mod sel {
    use super::*;

    pub const SECURITY_BITS: usize = 256;
    pub const SQISIGN_RESPONSE_LENGTH: usize = 253;
    pub const HASH_ITERATIONS: usize = 512;
    pub const PUBLICKEY_BYTES: usize = 129;
    pub const SECRETKEY_BYTES: usize = 701;
    pub const SIGNATURE_BYTES: usize = 292;

    pub const TORSION_EVEN_POWER: usize = 500;
    pub const P_COFACTOR_FOR_2F: [Digit; 1] = [27];
    pub const P_COFACTOR_FOR_2F_BITLENGTH: usize = 5;

    pub const TORSION_2POWER_BYTES: usize = 63;
    pub const TWO_TO_SECURITY_BITS: &[u64] = &[0x0, 0x0, 0x0, 0x0, 0x1];
    pub const TORSION_PLUS_2POWER: &[u64] = &[0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x10000000000000];
    pub const SEC_DEGREE: &[u64] = &[
        0x283, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1,
    ];
    pub const COM_DEGREE: &[u64] = &[
        0x283, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1,
    ];

    pub const BASIS_E0_PX: Fp2 = Fp2 {
        re: Fp([
            0xa6f4f854fc265d,
            0x4c8b84aa5fb427,
            0x1309a22f7f8bedd,
            0x12326c230bb7339,
            0x1177007f8e443b7,
            0x3227e897204471,
            0x173c12694021af7,
            0xd9af272f428697,
            0x523eda847ef5,
        ]),
        im: Fp([
            0xded07bbc792c63,
            0x11e1b26dec9cebc,
            0xc046644c2b6cd7,
            0x10fa781cd249b7d,
            0x1c100f6a2ab7eb,
            0x3268453a15b6a9,
            0x54d2827aa042c2,
            0x1976f2e8b7c96ec,
            0x16e01b2e8125f,
        ]),
    };
    pub const BASIS_E0_QX: Fp2 = Fp2 {
        re: Fp([
            0xf73af643285709,
            0xf149be2f088d45,
            0xcd261395ea3c0a,
            0x3a51f18f48bd2c,
            0x20878d18902069,
            0x1dde2b7d4cfad79,
            0x1cfc83af281db52,
            0xcb86b4138f7754,
            0xb4deb1e3f8a7,
        ]),
        im: Fp([
            0x924c7a12f1ab1e,
            0x37608c2f01a03,
            0x15ab8f95ccf5c3e,
            0x99e325091f7251,
            0xc375ef1b0a8b52,
            0x1e7185439fe829e,
            0x1393f18a069901e,
            0x1171a261ad16dd5,
            0x6573978c1c85,
        ]),
    };
}

pub use sel::*;
