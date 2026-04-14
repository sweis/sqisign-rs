// SPDX-License-Identifier: Apache-2.0
//! Precomputed signing-path constants (lvl1).
//! Auto-generated from the C reference by `tools/gen_precomp_sign.py`.
#![allow(clippy::all)]

use std::sync::OnceLock;
use crate::quaternion::{
    Ibz, IbzMat2x2, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal,
    QuatPExtremalMaximalOrder,
};
use crate::quaternion::intbig::{
    ibz_copy_digits, ibz_from_i64, ibz_mod, ibz_neg, ibz_pow, ibz_probab_prime,
    ibz_set_from_str,
};
use crate::ec::{EcCurve, EcBasis, EcPoint};
use crate::gf::{Fp, Fp2};

#[inline]
fn ibz_lit(neg: bool, limbs: &[u64]) -> Ibz {
    let mut z = Ibz::default();
    ibz_copy_digits(&mut z, limbs);
    if neg {
        let t = z.clone();
        ibz_neg(&mut z, &t);
    }
    z
}

// --- quaternion_constants.h ---
pub const QUAT_PRIMALITY_NUM_ITER: usize = 32;
pub const QUAT_REPRES_BOUND_INPUT: usize = 20;
pub const QUAT_EQUIV_BOUND_COEFF: usize = 64;
pub const FINDUV_BOX_SIZE: usize = 2;
pub const FINDUV_CUBE_SIZE: usize = 624;

// --- quaternion_data.c ---
pub fn quat_prime_cofactor() -> &'static Ibz {
    static V: OnceLock<Ibz> = OnceLock::new();
    V.get_or_init(|| ibz_lit(false, &[0x41,0x0,0x0,0x800000000000000]))
}
pub fn quatalg_pinfty() -> &'static QuatAlg {
    static V: OnceLock<QuatAlg> = OnceLock::new();
    V.get_or_init(|| QuatAlg { p: ibz_lit(false, &[0xffffffffffffffff,0xffffffffffffffff,0xffffffffffffffff,0x4ffffffffffffff]) })
}
pub const NUM_ALTERNATE_EXTREMAL_ORDERS: usize = 6;
pub const NUM_ALTERNATE_STARTING_CURVES: usize = 6;
pub fn extremal_orders() -> &'static [QuatPExtremalMaximalOrder; 7] {
    static V: OnceLock<[QuatPExtremalMaximalOrder; 7]> = OnceLock::new();
    V.get_or_init(|| [
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x2]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])],[Ibz::default(),ibz_lit(false, &[0x2]),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0x2]), coord: [Ibz::default(),ibz_lit(false, &[0x2]),Ibz::default(),Ibz::default()] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 1 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0x0,0x1000000000000000]), basis: [[ibz_lit(false, &[0x0,0x1000000000000000]),Ibz::default(),ibz_lit(false, &[0x0,0x800000000000000]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0x1]),Ibz::default(),ibz_lit(true, &[0x0,0x0,0x0,0x80000000000000])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x0,0x800000000000000]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0x1]),Ibz::default(),Ibz::default()]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0x0,0x1000000000000000]), coord: [Ibz::default(),ibz_lit(true, &[0x1]),Ibz::default(),ibz_lit(true, &[0x1])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 5 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0xf5f27a647b8578d4,0xb8746101369629b9]), basis: [[ibz_lit(false, &[0xf5f27a647b8578d4,0xb8746101369629b9]),Ibz::default(),ibz_lit(false, &[0xfaf93d323dc2bc6a,0x5c3a30809b4b14dc]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0x95ad2ad56fa47d47,0xc89877e749be8a4b,0x1]),Ibz::default(),ibz_lit(false, &[0x3e355e2970603f47,0x78dd10ae2a1bd950,0x0,0x280000000000000])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0xfaf93d323dc2bc6a,0x5c3a30809b4b14dc]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0x11]),Ibz::default(),ibz_lit(true, &[0xb19426e828ee3fe7,0xd6de568af586d7a])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0xf5f27a647b8578d4,0xb8746101369629b9]), coord: [Ibz::default(),ibz_lit(false, &[0x95ad2ad56fa47d47,0xc89877e749be8a4b,0x1]),Ibz::default(),ibz_lit(false, &[0x11])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 17 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0x3c6fa8e67715e5e2,0x17949bec872b9078]), basis: [[ibz_lit(false, &[0x3c6fa8e67715e5e2,0x17949bec872b9078]),Ibz::default(),ibz_lit(false, &[0x1e37d4733b8af2f1,0xbca4df64395c83c]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0xb034808274c8307a,0x9ab399ac43a4e8a]),Ibz::default(),ibz_lit(false, &[0x3d25ca466bc9954f,0x4f5822946ed431b,0xeb3e45306eb3e453,0x45306eb3e45306])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1e37d4733b8af2f1,0xbca4df64395c83c]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0x4]),Ibz::default(),ibz_lit(false, &[0xbd312454ca3a0e7f,0x2172f0cb4ce562])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0x3c6fa8e67715e5e2,0x17949bec872b9078]), coord: [Ibz::default(),ibz_lit(true, &[0xb034808274c8307a,0x9ab399ac43a4e8a]),Ibz::default(),ibz_lit(false, &[0x4])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 37 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0xde33c5116deeafa2,0x2df94f97c89ec8ce]), basis: [[ibz_lit(false, &[0xde33c5116deeafa2,0x2df94f97c89ec8ce]),Ibz::default(),ibz_lit(false, &[0x6f19e288b6f757d1,0x16fca7cbe44f6467]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0xd17aa943da6bdd36,0x44d44b0c564ce307]),Ibz::default(),ibz_lit(true, &[0xa0a2047cc4063a03,0x6cee07961df46dbc,0xc7ce0c7ce0c7ce0c,0x7ce0c7ce0c7ce0])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x6f19e288b6f757d1,0x16fca7cbe44f6467]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0x8]),Ibz::default(),ibz_lit(true, &[0xd9f82148a1e2188f,0xd6e1b21a072e79])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0xde33c5116deeafa2,0x2df94f97c89ec8ce]), coord: [Ibz::default(),ibz_lit(false, &[0xd17aa943da6bdd36,0x44d44b0c564ce307]),Ibz::default(),ibz_lit(true, &[0x8])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 41 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0x380014f2025b96a4,0x7bbeab7f79584e7c,0x1]), basis: [[ibz_lit(false, &[0x380014f2025b96a4,0x7bbeab7f79584e7c,0x1]),Ibz::default(),ibz_lit(false, &[0x1c000a79012dcb52,0xbddf55bfbcac273e]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0x4ba119e7333973e3,0xdbd0ee6227026ebc,0x7]),Ibz::default(),ibz_lit(false, &[0x9f01d923dd0ca33,0x83f7e395afe92f81,0xfffffffffffffffc,0x27fffffffffffff])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1c000a79012dcb52,0xbddf55bfbcac273e]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0x35]),Ibz::default(),ibz_lit(false, &[0x87f571c0f93ceb73,0x12fab9cbcb3c667a])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0x380014f2025b96a4,0x7bbeab7f79584e7c,0x1]), coord: [Ibz::default(),ibz_lit(true, &[0x4ba119e7333973e3,0xdbd0ee6227026ebc,0x7]),Ibz::default(),ibz_lit(false, &[0x35])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 53 },
            QuatPExtremalMaximalOrder { order: QuatLattice { denom: ibz_lit(false, &[0xe2b97b9e55af7ffa,0xc227f76b578ca7af,0xf]), basis: [[ibz_lit(false, &[0xe2b97b9e55af7ffa,0xc227f76b578ca7af,0xf]),Ibz::default(),ibz_lit(false, &[0xf15cbdcf2ad7bffd,0xe113fbb5abc653d7,0x7]),Ibz::default()],[Ibz::default(),ibz_lit(true, &[0xa2ef1ce7f02b0d16,0x66759632c56054b,0x6f]),Ibz::default(),ibz_lit(false, &[0x84ac06ea9d3bf0ab,0xd021882bdde962e5,0xffffffffffffffe2,0x13ffffffffffffff])],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0xf15cbdcf2ad7bffd,0xe113fbb5abc653d7,0x7]),Ibz::default()],[Ibz::default(),ibz_lit(false, &[0x308]),Ibz::default(),ibz_lit(false, &[0x77013f15c4a1f37,0x9281da3156007183])]] }, z: QuatAlgElem { denom: ibz_lit(false, &[0xe2b97b9e55af7ffa,0xc227f76b578ca7af,0xf]), coord: [Ibz::default(),ibz_lit(true, &[0xa2ef1ce7f02b0d16,0x66759632c56054b,0x6f]),Ibz::default(),ibz_lit(false, &[0x308])] }, t: QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()] }, q: 97 }
        ])
}
pub fn maxord_o0() -> &'static QuatLattice { &extremal_orders()[0].order }
pub fn standard_extremal_order() -> &'static QuatPExtremalMaximalOrder { &extremal_orders()[0] }
pub fn connecting_ideals() -> &'static [QuatLeftIdeal; 7] {
    static V: OnceLock<[QuatLeftIdeal; 7]> = OnceLock::new();
    V.get_or_init(|| [
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x2]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])],[Ibz::default(),ibz_lit(false, &[0x2]),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0x1]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x2,0x6000000000000000]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1,0x1000000000000000])],[Ibz::default(),ibz_lit(false, &[0x2,0x6000000000000000]),ibz_lit(false, &[0x1,0x5000000000000000]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0x1,0x3000000000000000]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x7f90157b8673f5fe,0x78f4a646d00bd2c5]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0xe65cd6d8002bfee5,0x5b1373de72d68a3])],[Ibz::default(),ibz_lit(false, &[0x7f90157b8673f5fe,0x78f4a646d00bd2c5]),ibz_lit(false, &[0x99333ea38647f719,0x73436f08e8de6a21]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0xbfc80abdc339faff,0x3c7a53236805e962]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x3c6fa8e67715e5e2,0x17949bec872b9078]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0xbb290a5a3af78597,0x84ff561d2d977c0])],[Ibz::default(),ibz_lit(false, &[0x3c6fa8e67715e5e2,0x17949bec872b9078]),ibz_lit(false, &[0x81469e8c3c1e604b,0xf44a68ab45218b7]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0x1e37d4733b8af2f1,0xbca4df64395c83c]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0xde33c5116deeafa2,0x2df94f97c89ec8ce]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0xd5f5cdcaa90b519b,0xe59b35483dd757a])],[Ibz::default(),ibz_lit(false, &[0xde33c5116deeafa2,0x2df94f97c89ec8ce]),ibz_lit(false, &[0x83df746c4e35e07,0x1f9f9c4344c15354]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0x6f19e288b6f757d1,0x16fca7cbe44f6467]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0x52a2ee77559419f2,0xb348218745c9f459]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1df48a96967adbd3,0x222419a0d707845])],[Ibz::default(),ibz_lit(false, &[0x52a2ee77559419f2,0xb348218745c9f459]),ibz_lit(false, &[0x34ae63e0bf193e1f,0xb125dfed38597c14]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0xa951773baaca0cf9,0x59a410c3a2e4fa2c]), parent_order: Some(maxord_o0()) },
            QuatLeftIdeal { lattice: QuatLattice { denom: ibz_lit(false, &[0x2]), basis: [[ibz_lit(false, &[0xd0316ad767cfaa3a,0x2996d852ebca0701]),Ibz::default(),Ibz::default(),ibz_lit(false, &[0xbc67edebd7ab0275,0x148ef2e5aeb5ad41])],[Ibz::default(),ibz_lit(false, &[0xd0316ad767cfaa3a,0x2996d852ebca0701]),ibz_lit(false, &[0x13c97ceb9024a7c5,0x1507e56d3d1459c0]),Ibz::default()],[Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1]),Ibz::default()],[Ibz::default(),Ibz::default(),Ibz::default(),ibz_lit(false, &[0x1])]] }, norm: ibz_lit(false, &[0xe818b56bb3e7d51d,0x14cb6c2975e50380]), parent_order: Some(maxord_o0()) }
        ])
}
pub fn conjugating_elements() -> &'static [QuatAlgElem; 7] {
    static V: OnceLock<[QuatAlgElem; 7]> = OnceLock::new();
    V.get_or_init(|| [
            QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [ibz_lit(false, &[0x1]),Ibz::default(),Ibz::default(),Ibz::default()] },
            QuatAlgElem { denom: ibz_lit(false, &[0x2]), coord: [ibz_lit(false, &[0x1,0x1000000000000000]),ibz_lit(true, &[0x1,0x1000000000000000]),ibz_lit(false, &[0x1]),ibz_lit(true, &[0x1])] },
            QuatAlgElem { denom: ibz_lit(false, &[0x2]), coord: [ibz_lit(false, &[0xcc7990f385eff94f,0x67e1008d1a8398d9]),ibz_lit(true, &[0xcc7990f385eff94f,0x67e1008d1a8398d9]),ibz_lit(true, &[0x3]),ibz_lit(false, &[0x3])] },
            QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [ibz_lit(false, &[0x1e37d4733b8af2f1,0xbca4df64395c83c]),ibz_lit(true, &[0x1e37d4733b8af2f1,0xbca4df64395c83c]),Ibz::default(),Ibz::default()] },
            QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [ibz_lit(false, &[0x6f19e288b6f757d1,0x16fca7cbe44f6467]),ibz_lit(true, &[0x6f19e288b6f757d1,0x16fca7cbe44f6467]),Ibz::default(),Ibz::default()] },
            QuatAlgElem { denom: ibz_lit(false, &[0x2]), coord: [ibz_lit(false, &[0xe869a36845fa6511,0xbdf3698988fc4db2]),ibz_lit(true, &[0xe869a36845fa6511,0xbdf3698988fc4db2]),ibz_lit(false, &[0x5]),ibz_lit(true, &[0x5])] },
            QuatAlgElem { denom: ibz_lit(false, &[0x1]), coord: [ibz_lit(false, &[0x77013f15c4a1e6b,0x9281da3156007183]),ibz_lit(true, &[0x77013f15c4a1e6b,0x9281da3156007183]),ibz_lit(true, &[0x4]),ibz_lit(false, &[0x4])] }
        ])
}

// --- endomorphism_action.c ---
/// Precomputed endomorphism rings applied to precomputed torsion bases.
#[derive(Clone)]
pub struct CurveWithEndomorphismRing {
    pub curve: EcCurve,
    pub basis_even: EcBasis,
    pub action_i: IbzMat2x2,
    pub action_j: IbzMat2x2,
    pub action_k: IbzMat2x2,
    pub action_gen2: IbzMat2x2,
    pub action_gen3: IbzMat2x2,
    pub action_gen4: IbzMat2x2,
}

pub fn curves_with_endomorphisms() -> &'static [CurveWithEndomorphismRing; 7] {
    static V: OnceLock<[CurveWithEndomorphismRing; 7]> = OnceLock::new();
    V.get_or_init(|| [
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x0,0x0,0x0,0x0,0x0]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0xc,0x0,0x0,0x0,0x400000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x5bcab12000c08,0x452654b56d052,0x26f81b5190a0a,0x36cfd66a361eb,0x12726610d11b]), im: Fp([0x6b96065c83efc,0x29da1d4a82cd9,0x190797ab98bdf,0x6841aa6eeee05,0x1377c5431166]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x21dd55b97832f,0x210f2d30b26ad,0x680bcfcf6396,0x27b318ec126a7,0x4ffba5956012]), im: Fp([0x74590149117e3,0x4982edefcc606,0x2ae3db0cc6884,0x7d0384872f5ec,0x4fbb0fcb5a52]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0xf6001dafb71a,0x75cb70989457f,0x5f2ab120f726c,0x7d12027e55817,0x6482fe24949]), im: Fp([0x63a39af1d2179,0x1c2884b0237f3,0x675979f836736,0x11de56ef443d1,0x462333fa18b7]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0xc5d3bda21b5456db,0x74759780861ddd06,0x7f9d34b241af33d1,0xcab471aa8c7f8c]),ibz_lit(false, &[0x7bfb7d32048b7d7a,0xa955918263d89bd3,0x76bf6861034403e1,0x574ae3eeb45cd0])],[ibz_lit(false, &[0x856fd6493698444f,0x189cafdf498f41db,0xf7e00bffe50bcb5b,0x1535daa88b47f9]),ibz_lit(false, &[0x3a2c425de4aba925,0x8b8a687f79e222f9,0x8062cb4dbe50cc2e,0x354b8e55738073])]], action_j: [[ibz_lit(false, &[0x36bad5fd54900abf,0xd14eea4a59da0f,0x914606f6a7aea3f0,0x7da2d2cde65004]),ibz_lit(false, &[0x611dbde3b7878680,0x819c9ec8b68a95f,0xbd7b5e31f73e2361,0x68240040d72b45])],[ibz_lit(false, &[0x1f0c9e126d204277,0x563f9d1cf854977f,0xe829af54c2ed00db,0xca7be80d8304fb]),ibz_lit(false, &[0xc9452a02ab6ff541,0xff2eb115b5a625f0,0x6eb9f90958515c0f,0x825d2d3219affb])]], action_k: [[ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6]),ibz_lit(false, &[0x55d2de69b685ad7a,0x925f591684e85675,0x83917c511cb68c0a,0xcd96ce11d1ffce])],[ibz_lit(false, &[0x959b1b9279bd3724,0x64a727d46f18b3ec,0x664bade78c7e9b4b,0x486a1da287a6d9]),ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59])]], action_gen2: [[ibz_lit(false, &[0xc5d3bda21b5456db,0x74759780861ddd06,0x7f9d34b241af33d1,0xcab471aa8c7f8c]),ibz_lit(false, &[0x7bfb7d32048b7d7a,0xa955918263d89bd3,0x76bf6861034403e1,0x574ae3eeb45cd0])],[ibz_lit(false, &[0x856fd6493698444f,0x189cafdf498f41db,0xf7e00bffe50bcb5b,0x1535daa88b47f9]),ibz_lit(false, &[0x3a2c425de4aba925,0x8b8a687f79e222f9,0x8062cb4dbe50cc2e,0x354b8e55738073])]], action_gen3: [[ibz_lit(false, &[0xfe4749cfb7f230cd,0xbaa37335683bdb8a,0x88719dd474aeebe0,0x242ba23c3967c8]),ibz_lit(false, &[0x6e8c9d8ade0981fd,0x58b7adb777a0a299,0x1a1d63497d4113a1,0xdfb77217c5c40b])],[ibz_lit(false, &[0x523e3a2dd1dc4363,0x376e267e20f1ecad,0xf004ddaa53fc661b,0x6fd8e15b07267a]),ibz_lit(false, &[0x1b8b630480dcf33,0x455c8cca97c42475,0x778e622b8b51141f,0xdbd45dc3c69837])]], action_gen4: [[ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3]),ibz_lit(false, &[0xaae96f34db42d6bd,0x492fac8b42742b3a,0x41c8be288e5b4605,0x66cb6708e8ffe7])],[ibz_lit(false, &[0x4acd8dc93cde9b92,0xb25393ea378c59f6,0xb325d6f3c63f4da5,0x24350ed143d36c]),ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x177f3bd3d98cf,0x568291dbf7092,0x755dcb3de2190,0x423388f314fe4,0x2a6f0241fb7]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x45dfcef4f6640,0x15a0a476fdc24,0x1d5772cf78864,0x708ce23cc53f9,0x2ca9bc0907ed]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x5f6259b797b43,0x157f63b3af2f9,0x7a3f4ea01dfa8,0x1dbb73e23680a,0x18914dc770b9]), im: Fp([0x8cb6e0ced492,0x5f20ac237154,0x7d25b71e8f3dd,0x4bf5fc15b1e6e,0x1dc3d80fa781]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x567ae7b4d67f3,0x5ccb6e9fa4f37,0x176489cb8f4ea,0x6a1c3c481062b,0x2c142d4feffe]), im: Fp([0x4c1bfcd30a39f,0x21b126ab96a61,0x60add76bd4a7,0x4a6a3d02240a9,0x1f52f1a6e758]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x23fad1a2013b,0x5e4194af99678,0x34468fab3bf1b,0x76e4e3f5b18c0,0x432503da9000]), im: Fp([0x34c912d2b3900,0x14d40850dcbe,0x672a3eab48ffe,0x2b790affecf8c,0x2ba92928eab]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0xe4058ceba8dcef13,0x3bbe28acfda5e2f5,0x5f5cb0ffee9141e5,0x95ef671e331920]),ibz_lit(false, &[0xb1b6fbce9e936b6e,0x6bcd20ae14b880bb,0xceb3c4a7feffb7f4,0xe9e00365bfd874])],[ibz_lit(false, &[0x523646b6c98847ff,0x7d56d563ec049694,0xe1958b0ac48f6833,0xdb58b2e957b64e]),ibz_lit(false, &[0x1bfa7314572310ed,0xc441d753025a1d0a,0xa0a34f00116ebe1a,0x6a1098e1cce6df])]], action_j: [[ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6]),ibz_lit(false, &[0x55d2de69b685ad7a,0x925f591684e85675,0x83917c511cb68c0a,0xcd96ce11d1ffce])],[ibz_lit(false, &[0x959b1b9279bd3724,0x64a727d46f18b3ec,0x664bade78c7e9b4b,0x486a1da287a6d9]),ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59])]], action_k: [[ibz_lit(false, &[0xe776f94c38f88d79,0x867742422d2e2bdf,0x8ee7a2e31736ddf0,0xa4bb554bb152ac]),ibz_lit(false, &[0xd49dc0b8e2806774,0x7a5dc53f25773b88,0x3ed5d6b24cfb3032,0xfc85b1584c27b8])],[ibz_lit(false, &[0xadb9d25b25cfc139,0x4e7a8867aa20bd39,0xacfc412aa81f8b24,0x201d50ab0cee2d]),ibz_lit(false, &[0x188906b3c7077287,0x7988bdbdd2d1d420,0x71185d1ce8c9220f,0x5b44aab44ead53])]], action_gen2: [[ibz_lit(false, &[0xe4058ceba8dcef13,0x3bbe28acfda5e2f5,0x5f5cb0ffee9141e5,0x95ef671e331920]),ibz_lit(false, &[0xb1b6fbce9e936b6e,0x6bcd20ae14b880bb,0xceb3c4a7feffb7f4,0xe9e00365bfd874])],[ibz_lit(false, &[0x523646b6c98847ff,0x7d56d563ec049694,0xe1958b0ac48f6833,0xdb58b2e957b64e]),ibz_lit(false, &[0x1bfa7314572310ed,0xc441d753025a1d0a,0xa0a34f00116ebe1a,0x6a1098e1cce6df])]], action_gen3: [[ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3]),ibz_lit(false, &[0xaae96f34db42d6bd,0x492fac8b42742b3a,0x41c8be288e5b4605,0x66cb6708e8ffe7])],[ibz_lit(false, &[0x4acd8dc93cde9b92,0xb25393ea378c59f6,0xb325d6f3c63f4da5,0x24350ed143d36c]),ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c])]], action_gen4: [[ibz_lit(false, &[0x994175298b307029,0x4553e3d77b3f2be8,0xc80bb49c7bef7065,0x181ece950cfa3e]),ibz_lit(false, &[0x7c8285e892ceb399,0x64f18924b18686eb,0x419631655e9a0d93,0x3155d501585e79])],[ibz_lit(false, &[0xaa0c720929f8da47,0x517c6e1939c9fc22,0x1edc20fccfa4c94e,0xdf85f0396de0d0]),ibz_lit(false, &[0x66be8ad674cf8fd7,0xbaac1c2884c0d417,0x37f44b6384108f9a,0xe7e1316af305c1])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x4d12b0e68b79f,0x337935267f3a8,0x380bf65840877,0x4bcc119304135,0x35da6e9613a8]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x1344ac39a2df4,0x6cde4d499fcea,0x2e02fd961021d,0x12f30464c104d,0x39769ba584ea]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x11012b71d2d54,0x76efaa195f3a3,0x6a89621403297,0xf05f07417877,0x58bafba5332]), im: Fp([0x3f3eaf5646a2d,0x6a0f369773854,0x15a15657d2442,0x667ba47d7dbf8,0x2d784590c43]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x3f45882691098,0x6a82534f3934f,0x6c6ead870b0ee,0x5669ed2bbb8da,0x2b9a1f281940]), im: Fp([0x41be7c586d896,0x22c68cb09ca5e,0x3c045adbe77b,0x506845058c043,0x2d2b7e8d71db]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x42b9c93c44402,0x461426db46e24,0x6d7aab066dc8c,0xbf26f540d0b8,0x4f6e2764cc0c]), im: Fp([0x72f03d7912cd,0x43aa6e7af9e21,0x679aa18a05871,0x14c0756affa95,0x2abcbd62f832]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0xe75d52b3a5945ff1,0xd9767d25d267dd09,0x10bf9aaec1a80bc5,0x70ae848de3e894]),ibz_lit(false, &[0xa6d796c0b9e011e6,0xf4c52f4404b6ee81,0xebb65b93e75d4597,0x163084c08e59c6])],[ibz_lit(false, &[0x479f60313463f41d,0x404e1e6b159f6fe7,0xad84c1f788a8e302,0xab3d6631758b50]),ibz_lit(false, &[0x18a2ad4c5a6ba00f,0x268982da2d9822f6,0xef4065513e57f43a,0x8f517b721c176b])]], action_j: [[ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6]),ibz_lit(false, &[0x55d2de69b685ad7a,0x925f591684e85675,0x83917c511cb68c0a,0xcd96ce11d1ffce])],[ibz_lit(false, &[0x959b1b9279bd3724,0x64a727d46f18b3ec,0x664bade78c7e9b4b,0x486a1da287a6d9]),ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59])]], action_k: [[ibz_lit(false, &[0x56c4c7ef1fbeffc3,0x1ced36aafa5c2834,0xe528890a31d9076,0xc298bcef6887d2]),ibz_lit(false, &[0x8109b5c6d7404098,0xc0d081c23a8c0299,0x656a89969243e848,0x4cb9f56c998c87])],[ibz_lit(false, &[0xc5fe55b5feed712b,0x7814177577a9e867,0x397386b173b14780,0xb001fa7f0b797a]),ibz_lit(false, &[0xa93b3810e041003d,0xe312c95505a3d7cb,0xf1ad776f5ce26f89,0x3d67431097782d])]], action_gen2: [[ibz_lit(false, &[0xe75d52b3a5945ff1,0xd9767d25d267dd09,0x10bf9aaec1a80bc5,0x70ae848de3e894]),ibz_lit(false, &[0xa6d796c0b9e011e6,0xf4c52f4404b6ee81,0xebb65b93e75d4597,0x163084c08e59c6])],[ibz_lit(false, &[0x479f60313463f41d,0x404e1e6b159f6fe7,0xad84c1f788a8e302,0xab3d6631758b50]),ibz_lit(false, &[0x18a2ad4c5a6ba00f,0x268982da2d9822f6,0xef4065513e57f43a,0x8f517b721c176b])]], action_gen3: [[ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3]),ibz_lit(false, &[0xaae96f34db42d6bd,0x492fac8b42742b3a,0x41c8be288e5b4605,0x66cb6708e8ffe7])],[ibz_lit(false, &[0x4acd8dc93cde9b92,0xb25393ea378c59f6,0xb325d6f3c63f4da5,0x24350ed143d36c]),ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c])]], action_gen4: [[ibz_lit(false, &[0x7e74cc3f1bd65d2b,0xd6d49f84fba04fea,0x4f4e68b188d142a1,0x2eb13ba60c13ec]),ibz_lit(false, &[0x649aaa1694487b4f,0x8df1d3fd3bc2e4b4,0xe8968caa4078931e,0x7c166ebed5deec])],[ibz_lit(false, &[0x2b40d32d51aa101d,0xb323257ac5f807ba,0x2c3ddc8f20c59bde,0x9917b954a64e7b]),ibz_lit(false, &[0x818b33c0e429a2d5,0x292b607b045fb015,0xb0b1974e772ebd5e,0xd14ec459f3ec13])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0xc17103986f53,0x6268ee5a8a215,0x11304cb0efe57,0x3846c2af6c518,0x2f57c43f40f7]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x2305c40e61be1,0x789a3b96a2885,0x44c132c3bf95,0x6e11b0abdb146,0x37d5f10fd03d]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x5b79ca4d5d6e0,0x39395e18e3349,0x75887ba6eb031,0x7d3b20412639b,0x13cf1bccb9dd]), im: Fp([0x76561c962386e,0x6f0884ce0b2e6,0x20dd8220aacb5,0x19375e2d543a7,0x4da1583c8553]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x4854b149c6d0c,0x7904efa1d89aa,0x343394a9e5c0f,0x68d9d640ad69d,0x2d711f0af96f]), im: Fp([0x3bcab7a6e1d94,0x6c35a91df0293,0x1b51f6ef1b777,0x6e9f0bb3d284,0x464e4d547390]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x376c342849596,0x657b69dced4b6,0x44b159aeb5eca,0x54b8abf1bdbfe,0x202393a746e4]), im: Fp([0x260478ad25e9b,0x21652ecc55014,0x728048f1594da,0x6b5eb728d6d3,0x3f305db59a7f]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0x415c44557ed2323f,0xcc1176ef42825876,0x340547291142bdab,0xc57c1f17791155]),ibz_lit(false, &[0x8a694faca958c9ce,0x8c191a17999731e1,0x8113c0eb68c7d118,0x3fc94ef8c862fd])],[ibz_lit(false, &[0x127c8ea0ed4741cb,0x3826ce8cf74c69d5,0xe695056b6bf33da2,0xd0581784acc45c]),ibz_lit(false, &[0xbea3bbaa812dcdc1,0x33ee8910bd7da789,0xcbfab8d6eebd4254,0x3a83e0e886eeaa])]], action_j: [[ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59]),ibz_lit(false, &[0xaa2d2196497a5286,0x6da0a6e97b17a98a,0x7c6e83aee34973f5,0x326931ee2e0031])],[ibz_lit(false, &[0x6a64e46d8642c8dc,0x9b58d82b90e74c13,0x99b45218738164b4,0xb795e25d785926]),ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6])]], action_k: [[ibz_lit(false, &[0xb6e0c901be7a7363,0xd659bc42779d6a56,0x923a12f438683476,0x3d5f954126fa8]),ibz_lit(false, &[0x14b5c144fd4edb4,0x4bb9c11d6bef702,0x5085b159de259f10,0x1dc4d36f42b0a9])],[ibz_lit(false, &[0x6b01c7ed4974e873,0x20d6c641f3d4affb,0x451e9f012d69ca22,0xdeb43bef65fa05]),ibz_lit(false, &[0x491f36fe41858c9d,0x29a643bd886295a9,0x6dc5ed0bc797cb89,0xfc2a06abed9057])]], action_gen2: [[ibz_lit(false, &[0x415c44557ed2323f,0xcc1176ef42825876,0x340547291142bdab,0xc57c1f17791155]),ibz_lit(false, &[0x8a694faca958c9ce,0x8c191a17999731e1,0x8113c0eb68c7d118,0x3fc94ef8c862fd])],[ibz_lit(false, &[0x127c8ea0ed4741cb,0x3826ce8cf74c69d5,0xe695056b6bf33da2,0xd0581784acc45c]),ibz_lit(false, &[0xbea3bbaa812dcdc1,0x33ee8910bd7da789,0xcbfab8d6eebd4254,0x3a83e0e886eeaa])]], action_gen3: [[ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c]),ibz_lit(false, &[0x551690cb24bd2943,0xb6d05374bd8bd4c5,0xbe3741d771a4b9fa,0x993498f7170018])],[ibz_lit(false, &[0xb5327236c321646e,0x4dac6c15c873a609,0x4cda290c39c0b25a,0xdbcaf12ebc2c93]),ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3])]], action_gen4: [[ibz_lit(false, &[0x490c34a61036ca5b,0x1c171590771bc0ed,0xdb988054977e4855,0xfe77221175b26f]),ibz_lit(false, &[0xb0cbae0a821cf543,0x2340d2bf5a80642d,0xdace4e38e1ce0c8f,0x59bc4807e3445])],[ibz_lit(false, &[0x39a70fd4c3dc6e85,0xb0fe0b83777b5158,0x53dc103f45b355de,0x12b18b6cb27e2]),ibz_lit(false, &[0xb6f3cb59efc935a5,0xe3e8ea6f88e43f12,0x24677fab6881b7aa,0x188ddee8a4d90])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x14612b0c4c481,0x7219e19939ca1,0x2bc69d2a0a8bd,0x5f4b0bcbad964,0x25664a8d484e]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x25184ac31312d,0x3c8678664e728,0xaf1a74a82a2f,0x57d2c2f2eb659,0xd5992a35213]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x58c095baf6ada,0x741ce646cd96,0x5007b4e8336a8,0x5010ebbfe93f9,0x1b2013c1eb92]), im: Fp([0x75c0724e94e91,0x77664d380f258,0xfb261c9ef941,0x749554a3cd77c,0x1b77c23de11f]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x71850cee2e1ca,0x1826b78a3cc19,0xddebf5154aa,0x696aeeba62d78,0x8953ba03b47]), im: Fp([0x2dc44634da928,0x4ea539513e1b6,0x5728c1bb241c3,0x3686f2152057e,0x2f6351277b8b]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x4c38023ba1341,0xee167e7a402b,0x7cbae09cd7aee,0x442bf312e4537,0x658d9f7ab76]), im: Fp([0x370f1db4d5016,0x4e773feecb28a,0x427c305ffbe2,0x687ab9f2e04cb,0x1feaa39f031c]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0x206ab453d052900d,0xfb21c57931f2e61d,0xf9c1f38f02bbc870,0xeb58d147f183aa]),ibz_lit(false, &[0x9e04ebc3a5e8727e,0x8ea968e038d7f1eb,0x82c048eb83318f77,0x54f2213583b0a3])],[ibz_lit(false, &[0x34e9c9b349e4dba9,0xcfb0cae0d8767ab9,0x1e302c9826b36177,0x713bdc53cc4a38]),ibz_lit(false, &[0xdf954bac2fad6ff3,0x4de3a86ce0d19e2,0x63e0c70fd44378f,0x14a72eb80e7c55])]], action_j: [[ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6]),ibz_lit(false, &[0x55d2de69b685ad7a,0x925f591684e85675,0x83917c511cb68c0a,0xcd96ce11d1ffce])],[ibz_lit(false, &[0x959b1b9279bd3724,0x64a727d46f18b3ec,0x664bade78c7e9b4b,0x486a1da287a6d9]),ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59])]], action_k: [[ibz_lit(false, &[0x89600cfe1b002417,0x222cc00f42d2662e,0xbcac863f278b7671,0x1b4c5a6e5edb9c]),ibz_lit(false, &[0xef865a2dd92b21e8,0x6b378ae01483f492,0xf4ec69c57b907f78,0xf8829616602fb9])],[ibz_lit(false, &[0xbc0e9aea5dc538ff,0xc311447b775dbea5,0x162a15fdb63af01c,0xc52a0d9defab76]),ibz_lit(false, &[0x769ff301e4ffdbe9,0xddd33ff0bd2d99d1,0x435379c0d874898e,0xe4b3a591a12463])]], action_gen2: [[ibz_lit(false, &[0x206ab453d052900d,0xfb21c57931f2e61d,0xf9c1f38f02bbc870,0xeb58d147f183aa]),ibz_lit(false, &[0x9e04ebc3a5e8727e,0x8ea968e038d7f1eb,0x82c048eb83318f77,0x54f2213583b0a3])],[ibz_lit(false, &[0x34e9c9b349e4dba9,0xcfb0cae0d8767ab9,0x1e302c9826b36177,0x713bdc53cc4a38]),ibz_lit(false, &[0xdf954bac2fad6ff3,0x4de3a86ce0d19e2,0x63e0c70fd44378f,0x14a72eb80e7c55])]], action_gen3: [[ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3]),ibz_lit(false, &[0xaae96f34db42d6bd,0x492fac8b42742b3a,0x41c8be288e5b4605,0x66cb6708e8ffe7])],[ibz_lit(false, &[0x4acd8dc93cde9b92,0xb25393ea378c59f6,0xb325d6f3c63f4da5,0x24350ed143d36c]),ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c])]], action_gen4: [[ibz_lit(false, &[0x731aacbf269320f0,0xf8361bcdd8ceb0f3,0xd3aad60444d58469,0x3f9086cdc34aa8]),ibz_lit(false, &[0x9534932fe55acc11,0x614f2956af432895,0x25560c6e4b24e84,0x13c61ed70dbb14])],[ibz_lit(false, &[0x54e4d24f450d28d6,0xadbe9b71081d6e67,0x4b684ebf61481088,0xe5e1a665c8829e]),ibz_lit(false, &[0x8ce55340d96cdf10,0x7c9e43227314f0c,0x2c5529fbbb2a7b96,0xc06f79323cb557])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x27e67b1ad4c35,0x4c9b9707ea7be,0x54e830f39a013,0x2661741eb0d4,0x40d297b19c53]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x49f99ec6b531a,0x7326e5c1fa9ef,0x153a0c3ce6804,0x609985d07ac35,0x1434a5ec6714]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x6292649ab6ec5,0x514c3aa63eaa8,0x42b95b0dce14a,0x5617e6b3d022,0x262a0b6ad948]), im: Fp([0x296936f8959c,0x7829b486d8303,0x51e4d11693064,0x3559dbc9d0dae,0x282ba45c8a46]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0xbd0e9751b3df,0x29bd7a6842bbd,0x61480930054f6,0x7c90f1cdb870a,0x10fc8988a92c]), im: Fp([0x6ee415f437e26,0x2244aa9d1a613,0x437f0b45ef3a9,0x749d8893337b5,0xe5a6eeb752b]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x85579e1b8722,0x632525e90080b,0x35539378e8d10,0x47389416f49d3,0x11c1e7bbb047]), im: Fp([0x367f0f2e5527c,0x763bfb94f5016,0x70df1a057bfdc,0x42460f20b8757,0x4a07f8d23dd8]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0xcd7513e0493127cb,0x9ff95a913de76846,0xb97226eca6d6a270,0x3f52fce4b80b44]),ibz_lit(false, &[0x1d16ca745a382d7e,0xafc28c2916742547,0x79572c7348562349,0xad04d33c3e67e1])],[ibz_lit(false, &[0xcf15321e27fbd8d7,0x7ed75fbd6f8efbd3,0xb73c593758d6f394,0x2264ace0270bfb]),ibz_lit(false, &[0x328aec1fb6ced835,0x6006a56ec21897b9,0x468dd91359295d8f,0xc0ad031b47f4bb])]], action_j: [[ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59]),ibz_lit(false, &[0xaa2d2196497a5286,0x6da0a6e97b17a98a,0x7c6e83aee34973f5,0x326931ee2e0031])],[ibz_lit(false, &[0x6a64e46d8642c8dc,0x9b58d82b90e74c13,0x99b45218738164b4,0xb795e25d785926]),ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6])]], action_k: [[ibz_lit(false, &[0xebb2cd7f6dc794df,0xc882825811db290,0xc8c37d64959ad514,0x321eea106a16a9]),ibz_lit(false, &[0x19fe1464e778e08c,0x4af0a98f1d24ef25,0x95ea2b828e4d70d3,0xf1695c2673277])],[ibz_lit(false, &[0x4bc9e2a4b6e1f1df,0x9a383fca6365dc85,0x1984ca7fed030ee2,0x8bc1731efee709]),ibz_lit(false, &[0x144d328092386b21,0xf377d7da7ee24d6f,0x373c829b6a652aeb,0xcde115ef95e956])]], action_gen2: [[ibz_lit(false, &[0xcd7513e0493127cb,0x9ff95a913de76846,0xb97226eca6d6a270,0x3f52fce4b80b44]),ibz_lit(false, &[0x1d16ca745a382d7e,0xafc28c2916742547,0x79572c7348562349,0xad04d33c3e67e1])],[ibz_lit(false, &[0xcf15321e27fbd8d7,0x7ed75fbd6f8efbd3,0xb73c593758d6f394,0x2264ace0270bfb]),ibz_lit(false, &[0x328aec1fb6ced835,0x6006a56ec21897b9,0x468dd91359295d8f,0xc0ad031b47f4bb])]], action_gen3: [[ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c]),ibz_lit(false, &[0x551690cb24bd2943,0xb6d05374bd8bd4c5,0xbe3741d771a4b9fa,0x993498f7170018])],[ibz_lit(false, &[0xb5327236c321646e,0x4dac6c15c873a609,0x4cda290c39c0b25a,0xdbcaf12ebc2c93]),ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3])]], action_gen4: [[ibz_lit(false, &[0x1ee9a31f187d647,0x9418bfedeec2193b,0xae854c740e9a15b6,0x423f1226773ebe]),ibz_lit(false, &[0x7bbf416c99be68ff,0xa8ff7682609fd44f,0xc0768e77ec03a6bb,0xf5a15296767873])],[ibz_lit(false, &[0xa38ebd23f7da3739,0x1a76b0e76908cf9,0x51015ac7a2bd77f0,0x952d3aa9223aae]),ibz_lit(false, &[0xfe1165ce0e7829b9,0x6be74012113de6c4,0x517ab38bf165ea49,0xbdc0edd988c141])]] },
            CurveWithEndomorphismRing { curve: EcCurve { a: Fp2 { re: Fp([0x1fd635b4f2c83,0x3ddd0240b9934,0x53881afe8d4a1,0x723f462627973,0x147962843332]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, c: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, a24: EcPoint { x: Fp2 { re: Fp([0x7f58d6d3cb2d,0x2f7740902e64d,0x74e206bfa3528,0x5c8fd18989e5c,0x311e58a10ccc]), im: Fp([0x0,0x0,0x0,0x0,0x0]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, is_a24_computed_and_normalized: true }, basis_even: EcBasis { p: EcPoint { x: Fp2 { re: Fp([0x472a0432a50a5,0x2584ec65ccf85,0x5a5586ba27eff,0x248f2f0f9bd37,0x42892709fd53]), im: Fp([0x3727bdaab80d,0x229e05a5546f4,0x4bad4d3212000,0x79e6087aee2df,0x42f9bfaf2bc8]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, q: EcPoint { x: Fp2 { re: Fp([0x140e00d2ad002,0x3235e1c701b8d,0x272d7237bc84d,0x44426d7ad2303,0x459a7fa89b08]), im: Fp([0x4246142cac789,0x1a160f97cc85d,0x43707cb72dff1,0x30e5aa57a2936,0x2c228ad830fe]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } }, pmq: EcPoint { x: Fp2 { re: Fp([0x519b1a003883d,0x356e25ed579a9,0x6b2a143d80555,0x1039d06c01ead,0xa3c331e0448]), im: Fp([0x45ddc052cdef3,0x20a40813439ef,0x52630baf0e697,0x4b49649819137,0x14d0e0cfb056]) }, z: Fp2 { re: Fp([0x19,0x0,0x0,0x0,0x300000000000]), im: Fp([0x0,0x0,0x0,0x0,0x0]) } } }, action_i: [[ibz_lit(false, &[0xc57273deb1867177,0xfe177031c0ee9802,0xed41e2a741c5bc2e,0x1ef5bc9ff91cbf]),ibz_lit(false, &[0x75c232b6bae3726a,0x382ad1726e79e003,0x6a39a56379628a51,0xa0f6f0c9109cdd])],[ibz_lit(false, &[0xbd8754e69fa9246b,0x25fa64701c7015b5,0x7eb5a6e989403f5c,0x16a8df54a16109]),ibz_lit(false, &[0x3a8d8c214e798e89,0x1e88fce3f1167fd,0x12be1d58be3a43d1,0xe10a436006e340])]], action_j: [[ibz_lit(false, &[0xb19c16401af2231b,0xf39a683ee470f713,0x904ec26e7a543289,0x4455fc6a0cd5a6]),ibz_lit(false, &[0x55d2de69b685ad7a,0x925f591684e85675,0x83917c511cb68c0a,0xcd96ce11d1ffce])],[ibz_lit(false, &[0x959b1b9279bd3724,0x64a727d46f18b3ec,0x664bade78c7e9b4b,0x486a1da287a6d9]),ibz_lit(false, &[0x4e63e9bfe50ddce5,0xc6597c11b8f08ec,0x6fb13d9185abcd76,0xbbaa0395f32a59])]], action_k: [[ibz_lit(false, &[0x9cbe086c2b021975,0x737ed9a7b1c37576,0xf9bf7652a2454de1,0x8eaa1dc2c4bf8]),ibz_lit(false, &[0x337c717746bcee88,0x3366b65740dc92b6,0x114640eb2b986c8a,0xe3a22fb00ae116])],[ibz_lit(false, &[0x40b9e24864d3f28d,0xcf3582ea82bb5141,0x6e88d71f0003faf0,0xcc6b9ef4c97ac9]),ibz_lit(false, &[0x6341f793d4fde68b,0x8c8126584e3c8a89,0x64089ad5dbab21e,0xf7155e23d3b407])]], action_gen2: [[ibz_lit(false, &[0xc57273deb1867177,0xfe177031c0ee9802,0xed41e2a741c5bc2e,0x1ef5bc9ff91cbf]),ibz_lit(false, &[0x75c232b6bae3726a,0x382ad1726e79e003,0x6a39a56379628a51,0xa0f6f0c9109cdd])],[ibz_lit(false, &[0xbd8754e69fa9246b,0x25fa64701c7015b5,0x7eb5a6e989403f5c,0x16a8df54a16109]),ibz_lit(false, &[0x3a8d8c214e798e89,0x1e88fce3f1167fd,0x12be1d58be3a43d1,0xe10a436006e340])]], action_gen3: [[ibz_lit(false, &[0xd8ce0b200d79118e,0xf9cd341f72387b89,0x482761373d2a1944,0x222afe35066ad3]),ibz_lit(false, &[0xaae96f34db42d6bd,0x492fac8b42742b3a,0x41c8be288e5b4605,0x66cb6708e8ffe7])],[ibz_lit(false, &[0x4acd8dc93cde9b92,0xb25393ea378c59f6,0xb325d6f3c63f4da5,0x24350ed143d36c]),ibz_lit(false, &[0x2731f4dff286ee73,0x632cbe08dc78476,0xb7d89ec8c2d5e6bb,0xddd501caf9952c])]], action_gen4: [[ibz_lit(false, &[0xc440bcc48ad184a0,0x784e15c646ea94e1,0x2bee0630d26f0190,0x3ce06193ce74b1]),ibz_lit(false, &[0xa69bfa45dafe1e2b,0xa0f9927df670b77e,0x5f229607c897ccb5,0xd9f781086747bf])],[ibz_lit(false, &[0x19bd02adbd3f2aa2,0xd17e3fff3b95e6f0,0x3b21da467888f3a6,0xc43e505301cc57]),ibz_lit(false, &[0x3bbf433b752e7b60,0x87b1ea39b9156b1e,0xd411f9cf2d90fe6f,0xc31f9e6c318b4e])]] }
        ])
}
pub fn curve_e0() -> &'static EcCurve { &curves_with_endomorphisms()[0].curve }
pub fn basis_even() -> &'static EcBasis { &curves_with_endomorphisms()[0].basis_even }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precomp::{P_COFACTOR_FOR_2F, TORSION_EVEN_POWER};

    #[test]
    fn pinfty_is_prime_p() {
        let mut expected = ibz_from_i64(P_COFACTOR_FOR_2F[0] as i64);
        let mut t = Ibz::default();
        ibz_pow(&mut t, crate::quaternion::ibz_const_two(), TORSION_EVEN_POWER as u32);
        let e = expected.clone();
        crate::quaternion::ibz_mul(&mut expected, &e, &t);
        let e = expected.clone();
        crate::quaternion::ibz_sub(&mut expected, &e, crate::quaternion::ibz_const_one());
        assert_eq!(quatalg_pinfty().p, expected);
        assert!(ibz_probab_prime(&quatalg_pinfty().p, 40) > 0);
    }

    #[test]
    fn extremal_order_0_is_o0() {
        let o0 = &extremal_orders()[0];
        assert_eq!(o0.order.denom, ibz_from_i64(2));
        assert_eq!(o0.q, 1);
        let exp: [[i64; 4]; 4] = [[2, 0, 0, 1], [0, 2, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(o0.order.basis[i][j], ibz_from_i64(exp[i][j]), "basis[{i}][{j}]");
            }
        }
        assert_eq!(o0.z.denom, ibz_from_i64(2));
        assert_eq!(o0.z.coord[1], ibz_from_i64(2));
        assert_eq!(o0.t.denom, ibz_from_i64(1));
        assert_eq!(o0.t.coord[2], ibz_from_i64(1));
    }

    #[test]
    fn connecting_ideal_0_is_unit() {
        let ci = &connecting_ideals()[0];
        assert_eq!(ci.norm, ibz_from_i64(1));
        assert!(ci.parent_order.is_some());
        assert!(std::ptr::eq(ci.parent_order.unwrap(), maxord_o0()));
    }

    #[test]
    fn curve_e0_is_a_zero() {
        let e0 = curve_e0();
        assert!(e0.is_a24_computed_and_normalized);
        assert_eq!(crate::gf::fp2_is_zero(&e0.a), 0xFFFF_FFFF);
        let p = &curves_with_endomorphisms()[0].basis_even.p;
        assert_eq!(p.x.re.0, crate::precomp::BASIS_E0_PX.re.0);
        assert_eq!(p.x.im.0, crate::precomp::BASIS_E0_PX.im.0);
    }

    #[test]
    fn endo_action_i_squares_to_neg_id() {
        let cwe = &curves_with_endomorphisms()[0];
        let m = &cwe.action_i;
        let mut modn = Ibz::default();
        ibz_pow(&mut modn, crate::quaternion::ibz_const_two(), TORSION_EVEN_POWER as u32);
        let r = |x: &Ibz| -> Ibz {
            let mut o = Ibz::default();
            ibz_mod(&mut o, x, &modn);
            o
        };
        let m2_00 = r(&(m[0][0].clone() * &m[0][0] + m[0][1].clone() * &m[1][0]));
        let m2_01 = r(&(m[0][0].clone() * &m[0][1] + m[0][1].clone() * &m[1][1]));
        let m2_10 = r(&(m[1][0].clone() * &m[0][0] + m[1][1].clone() * &m[1][0]));
        let m2_11 = r(&(m[1][0].clone() * &m[0][1] + m[1][1].clone() * &m[1][1]));
        let neg1 = &modn - ibz_from_i64(1);
        assert_eq!(m2_00, neg1);
        assert_eq!(m2_11, neg1);
        assert_eq!(m2_01, ibz_from_i64(0));
        assert_eq!(m2_10, ibz_from_i64(0));
    }

    #[test]
    fn all_populated() {
        let n = extremal_orders().len();
        assert_eq!(connecting_ideals().len(), n);
        assert_eq!(conjugating_elements().len(), n);
        assert_eq!(curves_with_endomorphisms().len(), n);
        assert!(extremal_orders()[n - 1].q > 1);
        assert!(connecting_ideals()[n - 1].norm > ibz_from_i64(1));
    }
}

#[cfg(test)]
mod tests_lvl1 {
    use super::*;

    #[test]
    fn prime_cofactor_value() {
        let mut expected = Ibz::default();
        ibz_set_from_str(
            &mut expected,
            "3618502788666131106986593281521497120414687020801267626233049500247285301313",
            10,
        );
        assert_eq!(*quat_prime_cofactor(), expected);
    }

    #[test]
    fn c_golden_cross_check() {
        assert_eq!(extremal_orders()[6].q, 97);
        assert_eq!(
            connecting_ideals()[6].norm.to_string(),
            "27640789963059351638707589454869550365"
        );
        assert_eq!(conjugating_elements()[6].coord[3], ibz_from_i64(4));
        assert_eq!(
            curves_with_endomorphisms()[0].action_gen4[1][1].to_string(),
            "391943321623591284286034922686343541417807403958897095695530545493876076147"
        );
        assert_eq!(
            curves_with_endomorphisms()[0].basis_even.q.x.re.0[0],
            0x21dd55b97832f
        );
    }
}
