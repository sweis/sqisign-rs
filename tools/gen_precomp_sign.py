#!/usr/bin/env python3
"""Parse the lvl1 quaternion_data.c and endomorphism_action.c precomputed
constant tables and emit Rust source.

We extract only the GMP_LIMB_BITS==64 and RADIX==64 non-BROADWELL branches,
flatten the C aggregate initializers into a token stream, then walk the known
struct layouts to produce Rust OnceLock-backed constructors.
"""
import re, sys, pathlib

LEVEL = sys.argv[1] if len(sys.argv) > 1 else "lvl1"
assert LEVEL in ("lvl1", "lvl3", "lvl5")
C_DIR = pathlib.Path(f"/root/src/personal-hacking/the-sqisign/src/precomp/ref/{LEVEL}")
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
        return "Ibz::default()"
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
    nw = {"lvl1": 5, "lvl3": 7, "lvl5": 9}[LEVEL]
    assert len(limbs) == nw, f"expected {nw} fp limbs, got {len(limbs)}"
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
//! Precomputed signing-path constants (""" + LEVEL + """).
//! Auto-generated from the C reference by `tools/gen_precomp_sign.py`.
#![allow(clippy::all)]
#![cfg_attr(rustfmt, rustfmt::skip)]

use std::sync::OnceLock;
use crate::quaternion::{
    Ibz, IbzMat2x2, QuatAlg, QuatAlgElem, QuatLattice, QuatLeftIdeal,
    QuatPExtremalMaximalOrder,
};
use crate::quaternion::intbig::{ibz_copy_digits, ibz_neg};
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

    def take_array(typename, name, taker):
        s.expect("const"); s.expect(typename)
        s.expect(name); s.expect("["); n = int(s.next()); s.expect("]"); s.expect("=")
        s.expect("{")
        items = []
        for _ in range(n):
            items.append(taker(s, ibzs))
            if s.peek() == ",": s.next()
        s.expect("}"); s.expect(";")
        return n, items

    n, eords = take_array("quat_p_extremal_maximal_order_t", "EXTREMAL_ORDERS",
                          take_quat_p_extremal_maximal_order)
    out.append(f"pub const NUM_ALTERNATE_EXTREMAL_ORDERS: usize = {n - 1};\n")
    out.append(f"pub const NUM_ALTERNATE_STARTING_CURVES: usize = {n - 1};\n")
    body = ",\n            ".join(eords)
    out.append(f"pub fn extremal_orders() -> &'static [QuatPExtremalMaximalOrder; {n}] {{\n"
               f"    static V: OnceLock<[QuatPExtremalMaximalOrder; {n}]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")
    out.append("pub fn maxord_o0() -> &'static QuatLattice { &extremal_orders()[0].order }\n")
    out.append("pub fn standard_extremal_order() -> &'static QuatPExtremalMaximalOrder { &extremal_orders()[0] }\n")

    n, cids = take_array("quat_left_ideal_t", "CONNECTING_IDEALS", take_quat_left_ideal_nostrict)
    body = ",\n            ".join(cids)
    out.append(f"pub fn connecting_ideals() -> &'static [QuatLeftIdeal; {n}] {{\n"
               f"    static V: OnceLock<[QuatLeftIdeal; {n}]> = OnceLock::new();\n"
               f"    V.get_or_init(|| [\n            {body}\n        ])\n}}\n")

    n, cels = take_array("quat_alg_elem_t", "CONJUGATING_ELEMENTS", take_quat_alg_elem)
    body = ",\n            ".join(cels)
    out.append(f"pub fn conjugating_elements() -> &'static [QuatAlgElem; {n}] {{\n"
               f"    static V: OnceLock<[QuatAlgElem; {n}]> = OnceLock::new();\n"
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
#[derive(Clone, Debug)]
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
    s.expect("CURVES_WITH_ENDOMORPHISMS"); s.expect("["); n = int(s.next()); s.expect("]"); s.expect("=")
    s.expect("{")
    cwe = []
    for _ in range(n):
        cwe.append(take_curve_with_endo(s, ibzs))
        if s.peek() == ",": s.next()
    s.expect("}"); s.expect(";")
    body = ",\n            ".join(cwe)
    out.append(f"pub fn curves_with_endomorphisms() -> &'static [CurveWithEndomorphismRing; {n}] {{\n"
               f"    static V: OnceLock<[CurveWithEndomorphismRing; {n}]> = OnceLock::new();\n"
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
    body += TESTS_GENERIC
    if LEVEL == "lvl1":
        body += TESTS_LVL1
    out = OUT_DIR / f"sign_data_{LEVEL}.rs"
    out.write_text(body)
    print(f"wrote {out} ({len(body)} bytes)")

TESTS_GENERIC = r"""
#[cfg(test)]
mod tests {
    use crate::quaternion::intbig::*;
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
        let neg1 = modn.clone() - ibz_from_i64(1);
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
"""

TESTS_LVL1 = r"""
#[cfg(test)]
mod tests_lvl1 {
    use crate::quaternion::intbig::*;
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
            ibz_convert_to_str(&connecting_ideals()[6].norm, 10).unwrap(),
            "27640789963059351638707589454869550365"
        );
        assert_eq!(conjugating_elements()[6].coord[3], ibz_from_i64(4));
        assert_eq!(
            ibz_convert_to_str(&curves_with_endomorphisms()[0].action_gen4[1][1], 10).unwrap(),
            "391943321623591284286034922686343541417807403958897095695530545493876076147"
        );
        assert_eq!(
            curves_with_endomorphisms()[0].basis_even.q.x.re.0[0],
            0x21dd55b97832f
        );
    }
}
"""

if __name__ == "__main__":
    main()
