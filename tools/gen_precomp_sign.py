#!/usr/bin/env python3
"""Parse the lvl1 quaternion_data.c and endomorphism_action.c precomputed
constant tables and emit Rust source.

We extract only the GMP_LIMB_BITS==64 and RADIX==64 non-BROADWELL branches,
flatten the C aggregate initializers into a token stream, then walk the known
struct layouts to produce Rust OnceLock-backed constructors.
"""
import re, sys, pathlib

C_DIR = pathlib.Path("/root/src/personal-hacking/the-sqisign/src/precomp/ref/lvl1")
OUT_DIR = pathlib.Path("/root/src/personal-hacking/sqisign-rs/src/precomp")

# ---------------------------------------------------------------------------
# Step 1: preprocess — keep only the 64-bit / non-broadwell #elif branches.

def preprocess(src: str) -> str:
    """Remove preprocessor conditionals, keeping only the 64-bit non-broadwell branch.
    Handles both `GMP_LIMB_BITS == 64` and `RADIX == 64` gating, and the nested
    `#if defined(SQISIGN_GF_IMPL_BROADWELL)` blocks."""
    lines = src.splitlines()
    out = []
    stack = []  # each entry: [emitting?, ever_emitted?]
    emitting = lambda: all(s[0] for s in stack)
    for ln in lines:
        s = ln.strip()
        if s.startswith("#if 0"):
            stack.append([False, False])
        elif s.startswith("#if defined(SQISIGN_GF_IMPL_BROADWELL)"):
            stack.append([False, False])
        elif s.startswith("#elif GMP_LIMB_BITS == 64") or s.startswith("#elif RADIX == 64"):
            stack[-1] = [not stack[-1][1], True]
        elif s.startswith("#elif"):
            stack[-1] = [False, stack[-1][1]]
        elif s.startswith("#else"):
            stack[-1] = [not stack[-1][1], True]
        elif s.startswith("#endif"):
            stack.pop()
        elif s.startswith("#include") or s.startswith("#ifndef") or s.startswith("#define"):
            pass
        else:
            if emitting():
                out.append(ln)
    return "\n".join(out)

# ---------------------------------------------------------------------------
# Step 2: tokenize the preprocessed C into a stream of leaf values.
#
# After preprocessing, the source is a sequence of:
#   const TYPE NAME = INIT ;
#   const TYPE NAME[N] = { INIT, INIT, ... } ;
# where INIT is nested {…} aggregates, ibz literals, hex limb arrays,
# integers, or `true`/`false`.

IBZ_RE = re.compile(
    r"\{\{\._mp_alloc\s*=\s*\d+,\s*\._mp_size\s*=\s*(-?\d+),\s*"
    r"\._mp_d\s*=\s*\(mp_limb_t\[\]\)\s*\{([^}]*)\}\}\}"
)

class Stream:
    def __init__(self, toks):
        self.toks = toks
        self.i = 0
    def next(self):
        t = self.toks[self.i]; self.i += 1; return t
    def expect(self, t):
        got = self.next()
        assert got == t, f"expected {t!r}, got {got!r} at idx {self.i-1}; context: {self.toks[max(0,self.i-5):self.i+5]}"
    def peek(self):
        return self.toks[self.i]

def tokenize(src: str):
    # First, lift ibz literals into atomic tokens so brace nesting inside
    # them doesn't confuse the structural tokenizer.
    ibzs = []
    def repl(m):
        sign = int(m.group(1))
        limbs = [x.strip() for x in m.group(2).split(",") if x.strip()]
        ibzs.append((sign, limbs))
        return f" @IBZ{len(ibzs)-1}@ "
    src = IBZ_RE.sub(repl, src)
    # Now tokenize: { } , ; = [ ] identifiers numbers @IBZk@
    toks = re.findall(r"@IBZ\d+@|[{},;=\[\]]|true|false|0x[0-9a-fA-F]+|\d+|[A-Za-z_][A-Za-z0-9_]*", src)
    return toks, ibzs

# ---------------------------------------------------------------------------
# Step 3: structure-aware emitters.

def emit_ibz(ibzs, tok):
    idx = int(tok[4:-1])
    sign, limbs = ibzs[idx]
    if sign == 0:
        return "Ibz::new()"
    neg = "true" if sign < 0 else "false"
    return f"ibz_lit({neg}, &[{','.join(limbs)}])"

def take_ibz(s, ibzs):
    return emit_ibz(ibzs, s.next())

def take_fp(s):
    # { 0x.., 0x.., 0x.., 0x.., 0x.. }
    s.expect("{")
    limbs = []
    while s.peek() != "}":
        limbs.append(s.next())
        if s.peek() == ",": s.next()
    s.expect("}")
    assert len(limbs) == 5, f"expected 5 fp limbs, got {len(limbs)}"
    return f"Fp([{','.join(limbs)}])"

def take_fp2(s):
    # { fp_re , fp_im }
    s.expect("{")
    re_ = take_fp(s); s.expect(",")
    im_ = take_fp(s)
    s.expect("}")
    return f"Fp2 {{ re: {re_}, im: {im_} }}"

def take_ec_point(s):
    # { fp2_x , fp2_z }
    s.expect("{")
    x = take_fp2(s); s.expect(",")
    z = take_fp2(s)
    s.expect("}")
    return f"EcPoint {{ x: {x}, z: {z} }}"

def take_ec_curve(s):
    # { fp2_A , fp2_C , ec_point_A24 , bool }
    s.expect("{")
    a = take_fp2(s); s.expect(",")
    c = take_fp2(s); s.expect(",")
    a24 = take_ec_point(s); s.expect(",")
    norm = s.next()
    s.expect("}")
    return (f"EcCurve {{ a: {a}, c: {c}, a24: {a24}, "
            f"is_a24_computed_and_normalized: {norm} }}")

def take_ec_basis(s):
    s.expect("{")
    p = take_ec_point(s); s.expect(",")
    q = take_ec_point(s); s.expect(",")
    pmq = take_ec_point(s)
    s.expect("}")
    return f"EcBasis {{ p: {p}, q: {q}, pmq: {pmq} }}"

def take_ibz_mat_2x2(s, ibzs):
    s.expect("{")
    rows = []
    for _ in range(2):
        s.expect("{")
        a = take_ibz(s, ibzs); s.expect(",")
        b = take_ibz(s, ibzs)
        s.expect("}")
        rows.append(f"[{a},{b}]")
        if s.peek() == ",": s.next()
    s.expect("}")
    return f"[{rows[0]},{rows[1]}]"

def take_ibz_mat_4x4(s, ibzs):
    s.expect("{")
    rows = []
    for _ in range(4):
        s.expect("{")
        cols = []
        for _ in range(4):
            cols.append(take_ibz(s, ibzs))
            if s.peek() == ",": s.next()
        s.expect("}")
        rows.append(f"[{','.join(cols)}]")
        if s.peek() == ",": s.next()
    s.expect("}")
    return f"[{','.join(rows)}]"

def take_ibz_vec_4(s, ibzs):
    s.expect("{")
    cols = []
    for _ in range(4):
        cols.append(take_ibz(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}")
    return f"[{','.join(cols)}]"

def take_quat_lattice(s, ibzs):
    s.expect("{")
    denom = take_ibz(s, ibzs); s.expect(",")
    basis = take_ibz_mat_4x4(s, ibzs)
    s.expect("}")
    return f"QuatLattice {{ denom: {denom}, basis: {basis} }}"

def take_quat_alg_elem(s, ibzs):
    s.expect("{")
    denom = take_ibz(s, ibzs); s.expect(",")
    coord = take_ibz_vec_4(s, ibzs)
    s.expect("}")
    return f"QuatAlgElem {{ denom: {denom}, coord: {coord} }}"

def take_quat_p_extremal_maximal_order(s, ibzs):
    s.expect("{")
    order = take_quat_lattice(s, ibzs); s.expect(",")
    z = take_quat_alg_elem(s, ibzs); s.expect(",")
    t = take_quat_alg_elem(s, ibzs); s.expect(",")
    q = s.next()
    s.expect("}")
    return f"QuatPExtremalMaximalOrder {{ order: {order}, z: {z}, t: {t}, q: {q} }}"

def take_quat_left_ideal(s, ibzs, parent):
    s.expect("{")
    lat = take_quat_lattice(s, ibzs); s.expect(",")
    norm = take_ibz(s, ibzs); s.expect(",")
    # parent_order is &EXTREMAL_ORDERS->order — emit as Some(ref)
    # consume `& IDENT -> order` or NULL
    tok = s.next()
    if tok == "NULL":
        po = "None"
    else:
        # token stream: & EXTREMAL_ORDERS -> order  (-> is 2 chars but our regex
        # split it; check)
        # Actually `&EXTREMAL_ORDERS->order` tokenizes as nothing standard.
        # The regex splits on word boundaries; `&` and `->` aren't in our set.
        # We need to handle them. Easiest: in tokenize(), also capture `&` and `->`.
        # We'll fix tokenize() and re-handle here.
        raise RuntimeError(f"unexpected parent_order token {tok!r}")
    s.expect("}")
    return f"QuatLeftIdeal {{ lattice: {lat}, norm: {norm}, parent_order: {po} }}"

def take_curve_with_endo(s, ibzs):
    s.expect("{")
    curve = take_ec_curve(s); s.expect(",")
    basis_even = take_ec_basis(s); s.expect(",")
    mats = []
    for _ in range(6):
        mats.append(take_ibz_mat_2x2(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}")
    return (f"CurveWithEndomorphismRing {{ curve: {curve}, basis_even: {basis_even}, "
            f"action_i: {mats[0]}, action_j: {mats[1]}, action_k: {mats[2]}, "
            f"action_gen2: {mats[3]}, action_gen3: {mats[4]}, action_gen4: {mats[5]} }}")

# ---------------------------------------------------------------------------
# Step 4: drivers per file.

HEADER = """\
// SPDX-License-Identifier: Apache-2.0
//! Precomputed signing-path constants (lvl1).
//! Auto-generated from the C reference by `tools/gen_precomp_sign.py`.
#![allow(clippy::all)]

use std::sync::OnceLock;
use rug::Integer;
use rug::integer::Order;
use crate::quaternion::{
    Ibz, IbzMat2x2, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal,
    QuatPExtremalMaximalOrder,
};
use crate::ec::{EcCurve, EcBasis, EcPoint};
use crate::gf::{Fp, Fp2};

#[inline]
fn ibz_lit(neg: bool, limbs: &[u64]) -> Ibz {
    let mut z = Integer::new();
    z.assign_digits(limbs, Order::Lsf);
    if neg { z = -z; }
    z
}

"""

def gen_quaternion_data():
    src = (C_DIR / "quaternion_data.c").read_text()
    src = preprocess(src)
    toks, ibzs = tokenize(src)
    s = Stream(toks)
    out = []

    # const ibz_t QUAT_prime_cofactor = IBZ ;
    s.expect("const"); s.expect("ibz_t"); s.expect("QUAT_prime_cofactor"); s.expect("=")
    pc = take_ibz(s, ibzs); s.expect(";")
    out.append(f"pub fn quat_prime_cofactor() -> &'static Ibz {{\n"
               f"    static V: OnceLock<Ibz> = OnceLock::new();\n"
               f"    V.get_or_init(|| {pc})\n}}\n")

    # const quat_alg_t QUATALG_PINFTY = { IBZ } ;
    s.expect("const"); s.expect("quat_alg_t"); s.expect("QUATALG_PINFTY"); s.expect("=")
    s.expect("{"); p = take_ibz(s, ibzs); s.expect("}"); s.expect(";")
    out.append(f"pub fn quatalg_pinfty() -> &'static QuatAlg {{\n"
               f"    static V: OnceLock<QuatAlg> = OnceLock::new();\n"
               f"    V.get_or_init(|| QuatAlg {{ p: {p} }})\n}}\n")

    # const quat_p_extremal_maximal_order_t EXTREMAL_ORDERS[7] = { ... } ;
    s.expect("const"); s.expect("quat_p_extremal_maximal_order_t")
    s.expect("EXTREMAL_ORDERS"); s.expect("["); s.expect("7"); s.expect("]"); s.expect("=")
    s.expect("{")
    eords = []
    for _ in range(7):
        eords.append(take_quat_p_extremal_maximal_order(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}"); s.expect(";")
    body = ",\n            ".join(eords)
    out.append(f"pub fn extremal_orders() -> &'static [QuatPExtremalMaximalOrder; 7] {{\n"
               f"    static V: OnceLock<[QuatPExtremalMaximalOrder; 7]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")
    out.append("pub fn maxord_o0() -> &'static QuatLattice { &extremal_orders()[0].order }\n")
    out.append("pub fn standard_extremal_order() -> &'static QuatPExtremalMaximalOrder { &extremal_orders()[0] }\n")

    # const quat_left_ideal_t CONNECTING_IDEALS[7] = { ... } ;
    s.expect("const"); s.expect("quat_left_ideal_t")
    s.expect("CONNECTING_IDEALS"); s.expect("["); s.expect("7"); s.expect("]"); s.expect("=")
    s.expect("{")
    cids = []
    for _ in range(7):
        cids.append(take_quat_left_ideal_nostrict(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}"); s.expect(";")
    body = ",\n            ".join(cids)
    out.append(f"pub fn connecting_ideals() -> &'static [QuatLeftIdeal; 7] {{\n"
               f"    static V: OnceLock<[QuatLeftIdeal; 7]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")

    # const quat_alg_elem_t CONJUGATING_ELEMENTS[7] = { ... } ;
    s.expect("const"); s.expect("quat_alg_elem_t")
    s.expect("CONJUGATING_ELEMENTS"); s.expect("["); s.expect("7"); s.expect("]"); s.expect("=")
    s.expect("{")
    cels = []
    for _ in range(7):
        cels.append(take_quat_alg_elem(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}"); s.expect(";")
    body = ",\n            ".join(cels)
    out.append(f"pub fn conjugating_elements() -> &'static [QuatAlgElem; 7] {{\n"
               f"    static V: OnceLock<[QuatAlgElem; 7]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")

    return "".join(out)

def take_quat_left_ideal_nostrict(s, ibzs):
    """parent_order in C is `&EXTREMAL_ORDERS->order`; our tokenizer drops the
    `&` and `->` punctuation, so the stream after norm is:
    EXTREMAL_ORDERS order } ...  — consume those two identifiers."""
    s.expect("{")
    lat = take_quat_lattice(s, ibzs); s.expect(",")
    norm = take_ibz(s, ibzs); s.expect(",")
    # Consume identifiers until '}'
    while s.peek() != "}":
        s.next()
    s.expect("}")
    return (f"QuatLeftIdeal {{ lattice: {lat}, norm: {norm}, "
            f"parent_order: Some(maxord_o0()) }}")

def gen_endomorphism_action():
    src = (C_DIR / "endomorphism_action.c").read_text()
    src = preprocess(src)
    toks, ibzs = tokenize(src)
    s = Stream(toks)
    out = []

    out.append("""\
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

""")

    s.expect("const"); s.expect("curve_with_endomorphism_ring_t")
    s.expect("CURVES_WITH_ENDOMORPHISMS"); s.expect("["); s.expect("7"); s.expect("]"); s.expect("=")
    s.expect("{")
    cwe = []
    for _ in range(7):
        cwe.append(take_curve_with_endo(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}"); s.expect(";")
    body = ",\n            ".join(cwe)
    out.append(f"pub fn curves_with_endomorphisms() -> &'static [CurveWithEndomorphismRing; 7] {{\n"
               f"    static V: OnceLock<[CurveWithEndomorphismRing; 7]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")
    out.append("pub fn curve_e0() -> &'static EcCurve { &curves_with_endomorphisms()[0].curve }\n")
    out.append("pub fn basis_even() -> &'static EcBasis { &curves_with_endomorphisms()[0].basis_even }\n")

    return "".join(out)

def gen_quaternion_constants():
    src = (C_DIR / "include/quaternion_constants.h").read_text()
    out = []
    for m in re.finditer(r"#define\s+(\w+)\s+(\d+)", src):
        name, val = m.group(1), m.group(2)
        out.append(f"pub const {name.upper()}: usize = {val};\n")
    return "".join(out)

def main():
    body = HEADER
    body += "// --- quaternion_constants.h ---\n"
    body += gen_quaternion_constants()
    body += "\n// --- quaternion_data.c ---\n"
    body += gen_quaternion_data()
    body += "\n// --- endomorphism_action.c ---\n"
    body += gen_endomorphism_action()
    body += "\npub const NUM_ALTERNATE_EXTREMAL_ORDERS: usize = 6;\n"
    body += "pub const NUM_ALTERNATE_STARTING_CURVES: usize = 6;\n"
    body += TESTS
    (OUT_DIR / "sign_data.rs").write_text(body)
    print(f"wrote {OUT_DIR/'sign_data.rs'} ({len(body)} bytes)")

TESTS = r"""
#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;
    use rug::ops::RemRounding;

    #[test]
    fn pinfty_is_prime_p() {
        // p = 5 * 2^248 - 1
        let mut expected = Integer::from(5);
        expected <<= 248;
        expected -= 1;
        assert_eq!(quatalg_pinfty().p, expected);
    }

    #[test]
    fn prime_cofactor_value() {
        // 2^251 + 65 (from limbs [0x41, 0, 0, 0x800000000000000])
        let mut expected = Integer::from(1);
        expected <<= 251;
        expected += 65;
        assert_eq!(*quat_prime_cofactor(), expected);
    }

    #[test]
    fn extremal_order_0_is_o0() {
        let o0 = &extremal_orders()[0];
        // Standard maximal order O₀: denom=2, basis=[[2,0,0,1],[0,2,1,0],[0,0,1,0],[0,0,0,1]], q=1.
        assert_eq!(o0.order.denom, 2);
        assert_eq!(o0.q, 1);
        let exp: [[i32; 4]; 4] = [[2,0,0,1],[0,2,1,0],[0,0,1,0],[0,0,0,1]];
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(o0.order.basis[i][j], exp[i][j], "basis[{i}][{j}]");
            }
        }
        // z = i (stored as 2i/2), t = j.
        assert_eq!(o0.z.denom, 2);
        assert_eq!(o0.z.coord[0], 0); assert_eq!(o0.z.coord[1], 2);
        assert_eq!(o0.z.coord[2], 0); assert_eq!(o0.z.coord[3], 0);
        assert_eq!(o0.t.denom, 1);
        assert_eq!(o0.t.coord[2], 1);
    }

    #[test]
    fn connecting_ideal_0_is_unit() {
        let ci = &connecting_ideals()[0];
        assert_eq!(ci.norm, 1);
        assert!(ci.parent_order.is_some());
        assert!(std::ptr::eq(ci.parent_order.unwrap(), maxord_o0()));
    }

    #[test]
    fn curve_e0_is_a_zero() {
        // E₀ has Montgomery A=0; its A24 is precomputed and normalized.
        let e0 = curve_e0();
        assert!(e0.is_a24_computed_and_normalized);
        assert_eq!(crate::gf::fp2_is_zero(&e0.a), 0xFFFF_FFFF);
        // basis_even.P.x must equal precomp::BASIS_E0_PX (cross-check between
        // the two precomp sources).
        let p = &curves_with_endomorphisms()[0].basis_even.p;
        assert_eq!(p.x.re.0, crate::precomp::BASIS_E0_PX.re.0);
        assert_eq!(p.x.im.0, crate::precomp::BASIS_E0_PX.im.0);
    }

    #[test]
    fn endo_action_i_squares_to_neg_id() {
        // The 2×2 action of i on the even-torsion basis must square to -1 mod 2^248.
        let cwe = &curves_with_endomorphisms()[0];
        let m = &cwe.action_i;
        let mut modn = Integer::from(1); modn <<= 248;
        let r = |x: &Integer| -> Integer { x.clone().rem_euc(modn.clone()) };
        let m2_00 = r(&(m[0][0].clone()*&m[0][0] + m[0][1].clone()*&m[1][0]));
        let m2_01 = r(&(m[0][0].clone()*&m[0][1] + m[0][1].clone()*&m[1][1]));
        let m2_10 = r(&(m[1][0].clone()*&m[0][0] + m[1][1].clone()*&m[1][0]));
        let m2_11 = r(&(m[1][0].clone()*&m[0][1] + m[1][1].clone()*&m[1][1]));
        let neg1 = modn.clone() - 1;
        assert_eq!(m2_00, neg1);
        assert_eq!(m2_11, neg1);
        assert_eq!(m2_01, 0);
        assert_eq!(m2_10, 0);
    }

    #[test]
    fn c_golden_cross_check() {
        // Values dumped from the C build (tools: /tmp/dump_precomp.c).
        assert_eq!(extremal_orders()[6].q, 97);
        assert_eq!(
            connecting_ideals()[6].norm.to_string(),
            "27640789963059351638707589454869550365"
        );
        assert_eq!(conjugating_elements()[6].coord[3], 4);
        assert_eq!(
            curves_with_endomorphisms()[0].action_gen4[1][1].to_string(),
            "391943321623591284286034922686343541417807403958897095695530545493876076147"
        );
        assert_eq!(
            curves_with_endomorphisms()[0].basis_even.q.x.re.0[0],
            0x21dd55b97832f
        );
    }

    #[test]
    fn all_seven_populated() {
        assert_eq!(extremal_orders().len(), 7);
        assert_eq!(connecting_ideals().len(), 7);
        assert_eq!(conjugating_elements().len(), 7);
        assert_eq!(curves_with_endomorphisms().len(), 7);
        // Spot-check that later entries have non-trivial data.
        assert!(extremal_orders()[6].q > 1);
        assert!(connecting_ideals()[6].norm > 1);
    }
}
"""

if __name__ == "__main__":
    main()
