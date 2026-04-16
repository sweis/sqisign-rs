//! Prototype AVX-512 IFMA backend for GF(p), p = 5·2²⁴⁸ − 1 (lvl1).
//!
//! Strategy A: single field element packed into one zmm as 5 unsaturated
//! radix-2⁵² limbs (lanes 0–4; lanes 5–7 zero). Product accumulated with
//! `vpmadd52lo/hi`, then column-serial Montgomery reduction (R = 2²⁶⁰).
//! This is a standalone prototype: it has its own Montgomery domain and
//! converts at the byte boundary, so correctness is checked against the
//! production `Fp` via `encode`/`decode`.
//!
//! Gated on `target_feature = "avx512ifma"` (enable with
//! `RUSTFLAGS="-C target-cpu=native"` on Ice Lake / Sapphire Rapids).
#![cfg(all(
    feature = "lvl1",
    target_arch = "x86_64",
    target_feature = "avx512ifma",
    target_feature = "avx512f"
))]
#![allow(unsafe_code, clippy::cast_possible_truncation)]

use core::arch::x86_64::*;

pub const LIMB_BITS: u32 = 52;
pub const NLIMBS: usize = 5;
pub const MASK52: u64 = (1u64 << LIMB_BITS) - 1;
/// Top limb of (p+1) in radix-2⁵²: (5·2²⁴⁸)/2²⁰⁸ = 5·2⁴⁰.
pub const P_TOP_52: u64 = 5u64 << 40;
/// 2²⁶⁰ mod p = 2²⁴⁸ + 819 (used to fold limb-4 overflow back).
pub const TWO_260_MOD_P_LO: u64 = 819;
pub const TWO_260_MOD_P_LIMB4: u64 = 1u64 << 40;

/// One field element, radix-2⁵², Montgomery domain R = 2²⁶⁰.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(align(64))]
pub struct FpIfma(pub [u64; 8]); // lanes 0-4 used, 5-7 zero

impl FpIfma {
    pub const ZERO: Self = Self([0; 8]);

    #[inline(always)]
    fn load(&self) -> __m512i {
        unsafe { _mm512_load_si512(self.0.as_ptr().cast()) }
    }

    /// Propagate carries so every limb < 2⁵² (required before IFMA mul).
    #[inline(always)]
    pub fn normalize(&mut self) {
        let mut c = 0u64;
        for i in 0..5 {
            let t = self.0[i] + c;
            self.0[i] = t & MASK52;
            c = t >> LIMB_BITS;
        }
        // Fold limb-5 overflow back via 2²⁶⁰ ≡ 2²⁴⁸ + 819.
        // c is at most a few bits; one extra pass is enough.
        self.0[0] += c * TWO_260_MOD_P_LO;
        self.0[4] += c * TWO_260_MOD_P_LIMB4;
        let t = self.0[0];
        self.0[0] = t & MASK52;
        self.0[1] += t >> LIMB_BITS;
    }

    /// Canonical 32-byte LE → radix-2⁵² limbs (NOT Montgomery).
    pub fn from_le_bytes_raw(b: &[u8; 32]) -> Self {
        let mut w = [0u64; 4];
        for (i, c) in b.chunks_exact(8).enumerate() {
            w[i] = u64::from_le_bytes(c.try_into().unwrap());
        }
        let mut r = Self::ZERO;
        r.0[0] = w[0] & MASK52;
        r.0[1] = (w[0] >> 52 | w[1] << 12) & MASK52;
        r.0[2] = (w[1] >> 40 | w[2] << 24) & MASK52;
        r.0[3] = (w[2] >> 28 | w[3] << 36) & MASK52;
        r.0[4] = w[3] >> 16; // ≤ 48 bits since input < 2²⁵⁶
        r
    }

    /// Radix-2⁵² limbs (lazy, integer value may be < ~4p) → fully reduced 32-byte LE.
    pub fn to_le_bytes_canonical(mut self) -> [u8; 32] {
        // Full carry-propagate then a couple of conditional p-subtractions.
        self.normalize();
        // After normalize, integer < ~5·2²⁵² ≈ 4p; subtract p while ≥ p.
        // p limbs: [M, M, M, M, 5·2⁴⁰-1].
        let p4 = P_TOP_52 - 1;
        for _ in 0..3 {
            // trial subtract p
            let mut borrow = 0i128;
            let mut t = [0u64; 5];
            let p = [MASK52, MASK52, MASK52, MASK52, p4];
            for i in 0..5 {
                let v = self.0[i] as i128 - p[i] as i128 - borrow;
                t[i] = (v as u64) & MASK52;
                borrow = i128::from(v < 0);
            }
            if borrow == 0 {
                self.0[..5].copy_from_slice(&t);
            }
        }
        // pack 5×52 → 4×64
        let mut w = [0u64; 4];
        w[0] = self.0[0] | self.0[1] << 52;
        w[1] = self.0[1] >> 12 | self.0[2] << 40;
        w[2] = self.0[2] >> 24 | self.0[3] << 28;
        w[3] = self.0[3] >> 36 | self.0[4] << 16;
        let mut out = [0u8; 32];
        for (i, c) in out.chunks_exact_mut(8).enumerate() {
            c.copy_from_slice(&w[i].to_le_bytes());
        }
        out
    }

    /// Canonical bytes → Montgomery (×R via montmul with R²).
    pub fn from_le_bytes_montgomery(b: &[u8; 32]) -> Self {
        let raw = Self::from_le_bytes_raw(b);
        fp_mul_ifma(&raw, &R2_52)
    }

    /// Montgomery → canonical bytes (montmul with 1, i.e., one MontRed).
    pub fn encode(self) -> [u8; 32] {
        let one = {
            let mut o = Self::ZERO;
            o.0[0] = 1;
            o
        };
        fp_mul_ifma(&self, &one).to_le_bytes_canonical()
    }
}

// ---------------------------------------------------------------------------
// Reference: scalar radix-2⁵² Montgomery mul (for correctness oracle).
// Same column-serial algorithm as modarith but with LIMB_BITS=52.
// ---------------------------------------------------------------------------

pub fn fp_mul_ref52(a: &FpIfma, b: &FpIfma) -> FpIfma {
    let p4 = P_TOP_52 as u128;
    let a = &a.0;
    let b = &b.0;
    let mut t: u128 = 0;
    let mut v = [0u64; 5];
    let mut c = [0u64; 8];
    macro_rules! col {
        ($k:expr; $($i:expr,$j:expr);*) => {{
            $( t += (a[$i] as u128) * (b[$j] as u128); )*
            if $k >= 4 { t += (v[$k-4] as u128) * p4; }
            let lo = (t as u64) & MASK52;
            t >>= LIMB_BITS;
            lo
        }};
    }
    v[0] = col!(0; 0,0);
    v[1] = col!(1; 0,1; 1,0);
    v[2] = col!(2; 0,2; 1,1; 2,0);
    v[3] = col!(3; 0,3; 1,2; 2,1; 3,0);
    v[4] = col!(4; 0,4; 1,3; 2,2; 3,1; 4,0);
    c[0] = col!(5; 1,4; 2,3; 3,2; 4,1);
    c[1] = col!(6; 2,4; 3,3; 4,2);
    c[2] = col!(7; 3,4; 4,3);
    c[3] = col!(8; 4,4);
    c[4] = t as u64;
    let mut r = FpIfma(c);
    r.normalize();
    r
}

// ---------------------------------------------------------------------------
// IFMA: row-major schoolbook into 8-lane accumulator + 2 scalar overflow limbs,
// then scalar column-serial Montgomery reduction.
// ---------------------------------------------------------------------------

/// `a · b · R⁻¹ mod p` (R = 2²⁶⁰), inputs must have limbs < 2⁵².
#[inline]
pub fn fp_mul_ifma(a: &FpIfma, b: &FpIfma) -> FpIfma {
    unsafe {
        let zero = _mm512_setzero_si512();
        let av = a.load();
        // Six left-shifted copies of a: a in lanes k..k+4 (a₄ dropped for k≥4).
        // valignq(hi, lo, n) = (hi:lo) >> (n·64); to shift a left by k lanes
        // (zeros in low lanes), use valignq(av, zero, 8-k).
        let a0 = av; // lanes 0-4
        let a1 = _mm512_alignr_epi64(av, zero, 7); // lanes 1-5
        let a2 = _mm512_alignr_epi64(av, zero, 6); // lanes 2-6
        let a3 = _mm512_alignr_epi64(av, zero, 5); // lanes 3-7
        let a4 = _mm512_alignr_epi64(av, zero, 4); // lanes 4-7 (a₀..a₃), a₄ dropped
        let a5 = _mm512_alignr_epi64(av, zero, 3); // lanes 5-7 (a₀..a₂)

        let bb = b.0;
        let b0 = _mm512_set1_epi64(bb[0] as i64);
        let b1 = _mm512_set1_epi64(bb[1] as i64);
        let b2 = _mm512_set1_epi64(bb[2] as i64);
        let b3 = _mm512_set1_epi64(bb[3] as i64);
        let b4 = _mm512_set1_epi64(bb[4] as i64);

        // acc lanes 0-7 accumulate product limbs 0-7.
        let mut acc = _mm512_madd52lo_epu64(zero, a0, b0); // limbs 0-4
        acc = _mm512_madd52hi_epu64(acc, a1, b0); // limbs 1-5
        acc = _mm512_madd52lo_epu64(acc, a1, b1); // limbs 1-5
        acc = _mm512_madd52hi_epu64(acc, a2, b1); // limbs 2-6
        acc = _mm512_madd52lo_epu64(acc, a2, b2); // limbs 2-6
        acc = _mm512_madd52hi_epu64(acc, a3, b2); // limbs 3-7
        acc = _mm512_madd52lo_epu64(acc, a3, b3); // limbs 3-7
        acc = _mm512_madd52hi_epu64(acc, a4, b3); // limbs 4-7 (a₄·b₃ hi → limb 8: scalar)
        acc = _mm512_madd52lo_epu64(acc, a4, b4); // limbs 4-7 (a₄·b₄ lo → limb 8: scalar)
        acc = _mm512_madd52hi_epu64(acc, a5, b4); // limbs 5-7 (a₃,a₄·b₄ hi → limb 8,9: scalar)

        // Scalar overflow contributions to limbs 8-9.
        let aa4 = a.0[4];
        let p_a4b3 = (aa4 as u128) * (bb[3] as u128);
        let p_a4b4 = (aa4 as u128) * (bb[4] as u128);
        let p_a3b4 = (a.0[3] as u128) * (bb[4] as u128);
        let z8 = (p_a4b3 >> 52) as u64 + (p_a4b4 as u64 & MASK52) + (p_a3b4 >> 52) as u64;
        let z9 = (p_a4b4 >> 52) as u64;

        // Spill acc to z[0..7] for the (serial) Montgomery reduction.
        let mut z = [0u64; 10];
        _mm512_storeu_si512(z.as_mut_ptr().cast(), acc);
        z[8] = z8;
        z[9] = z9;

        montred_5x52(&z)
    }
}

/// Column-serial Montgomery reduction of a 10-limb product to 5 limbs.
#[inline(always)]
fn montred_5x52(z: &[u64; 10]) -> FpIfma {
    let p4 = P_TOP_52 as u128;
    let mut t: u128 = 0;
    let mut v = [0u64; 5];
    let mut c = FpIfma::ZERO;
    for k in 0..9 {
        t += z[k] as u128;
        if k >= 4 {
            t += (v[k - 4] as u128) * p4;
        }
        let lo = (t as u64) & MASK52;
        if k < 5 {
            v[k] = lo;
        } else {
            c.0[k - 5] = lo;
        }
        t >>= LIMB_BITS;
    }
    t += z[9] as u128;
    c.0[4] = t as u64;
    c.normalize();
    c
}

// ---------------------------------------------------------------------------
// Strategy B: batched 8-way mul. Eight independent Fp elements in SoA layout:
// limb i of all 8 elements packs into one __m512i. No permutes; the full
// SIMD width is used. This is where IFMA actually wins.
// ---------------------------------------------------------------------------

/// Eight Fp elements in radix-2⁵² SoA: `limbs[i]` lane k = element k's limb i.
#[allow(missing_debug_implementations)]
#[derive(Clone, Copy)]
#[repr(align(64))]
pub struct FpIfma8 {
    pub limbs: [__m512i; 5],
}

impl Default for FpIfma8 {
    fn default() -> Self {
        unsafe { core::mem::zeroed() }
    }
}

// ---------------------------------------------------------------------------
// Radix-64 ↔ radix-52 SoA bridge.
//
// The production `Fp` (GF5_248) holds y in [0, 2²⁵¹) with y ≡ x·2²⁵⁶ mod p.
// Bit-splitting y to radix-52 gives the same integer; an IFMA mul then yields
// `a·b·2⁻²⁶⁰` instead of the desired `a·b·2⁻²⁵⁶`. We correct by left-shifting
// one operand by 4 bits *before* the bit-split (so the product gains a 2⁴
// factor). With both inputs < 2²⁵² (the `add_noreduce` bound), the shifted
// operand is < 2²⁵⁶ so its limb-4 stays < 2⁴⁸, and `fp_mul_ifma8`'s output is
// then < `2²⁵²·2²⁵⁶/2²⁶⁰ + p` < 2²⁵¹ — back in `GF5_248`'s representable
// range after bit-join.
// ---------------------------------------------------------------------------

/// 8×4 u64 AoS → 4× zmm SoA (lane k of `out[i]` = `e[k][i]`).
#[inline(always)]
unsafe fn transpose_8x4(e: &[[u64; 4]; 8]) -> [__m512i; 4] {
    let r0 = _mm512_loadu_si512(e[0].as_ptr().cast()); // e0 e1
    let r1 = _mm512_loadu_si512(e[2].as_ptr().cast()); // e2 e3
    let r2 = _mm512_loadu_si512(e[4].as_ptr().cast()); // e4 e5
    let r3 = _mm512_loadu_si512(e[6].as_ptr().cast()); // e6 e7
                                                       // Interleave 64-bit lanes within 128-bit pairs.
    let a0 = _mm512_unpacklo_epi64(r0, r1); // e0.0 e2.0 e0.2 e2.2 e1.0 e3.0 e1.2 e3.2
    let a1 = _mm512_unpackhi_epi64(r0, r1); // e0.1 e2.1 e0.3 e2.3 e1.1 e3.1 e1.3 e3.3
    let a2 = _mm512_unpacklo_epi64(r2, r3);
    let a3 = _mm512_unpackhi_epi64(r2, r3);
    // Gather all limb-0 lanes etc. with explicit 16-way permutes.
    static IDX_LO: [u64; 8] = [0, 4, 1, 5, 8, 12, 9, 13];
    static IDX_HI: [u64; 8] = [2, 6, 3, 7, 10, 14, 11, 15];
    let lo = _mm512_loadu_si512(IDX_LO.as_ptr().cast());
    let hi = _mm512_loadu_si512(IDX_HI.as_ptr().cast());
    [
        _mm512_permutex2var_epi64(a0, lo, a2),
        _mm512_permutex2var_epi64(a1, lo, a3),
        _mm512_permutex2var_epi64(a0, hi, a2),
        _mm512_permutex2var_epi64(a1, hi, a3),
    ]
}

/// Bit-split eight 4×64-bit integers into SoA radix-2⁵². No arithmetic.
#[inline(always)]
fn pack_u64x4(e: &[[u64; 4]; 8]) -> FpIfma8 {
    unsafe {
        let m = _mm512_set1_epi64(MASK52 as i64);
        let [w0, w1, w2, w3] = transpose_8x4(e);
        FpIfma8 {
            limbs: [
                _mm512_and_si512(w0, m),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<52>(w0), _mm512_slli_epi64::<12>(w1)),
                    m,
                ),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<40>(w1), _mm512_slli_epi64::<24>(w2)),
                    m,
                ),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<28>(w2), _mm512_slli_epi64::<36>(w3)),
                    m,
                ),
                _mm512_srli_epi64::<16>(w3),
            ],
        }
    }
}

/// As `pack_u64x4`, but pre-shift each integer left by 4 bits (the R₆₄/R₅₂
/// domain-correction factor). Inputs must be < 2²⁵².
#[inline(always)]
fn pack_u64x4_x16(e: &[[u64; 4]; 8]) -> FpIfma8 {
    unsafe {
        let m = _mm512_set1_epi64(MASK52 as i64);
        let [w0, w1, w2, w3] = transpose_8x4(e);
        FpIfma8 {
            limbs: [
                _mm512_and_si512(_mm512_slli_epi64::<4>(w0), m),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<48>(w0), _mm512_slli_epi64::<16>(w1)),
                    m,
                ),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<36>(w1), _mm512_slli_epi64::<28>(w2)),
                    m,
                ),
                _mm512_and_si512(
                    _mm512_or_si512(_mm512_srli_epi64::<24>(w2), _mm512_slli_epi64::<40>(w3)),
                    m,
                ),
                _mm512_srli_epi64::<12>(w3),
            ],
        }
    }
}

/// Carry-propagate, fold to < 2²⁵¹, bit-join SoA radix-2⁵² → eight 4×64-bit
/// integers (each a valid `GF5_248` representative).
#[inline]
fn unpack_u64x4(mut a: FpIfma8) -> [[u64; 4]; 8] {
    unsafe {
        let mask = _mm512_set1_epi64(MASK52 as i64);
        let mut c = _mm512_setzero_si512();
        for i in 0..5 {
            let t = _mm512_add_epi64(a.limbs[i], c);
            a.limbs[i] = _mm512_and_si512(t, mask);
            c = _mm512_srli_epi64::<52>(t);
        }
        debug_assert_eq!(_mm512_reduce_or_epi64(c), 0);
    }
    let mut buf = [[0u64; 8]; 5];
    for i in 0..5 {
        unsafe { _mm512_storeu_si512(buf[i].as_mut_ptr().cast(), a.limbs[i]) };
    }
    let mut out = [[0u64; 4]; 8];
    for k in 0..8 {
        out[k] = join_and_reduce(buf[0][k], buf[1][k], buf[2][k], buf[3][k], buf[4][k]);
    }
    out
}

/// Bit-join five radix-2⁵² limbs (each < 2⁵²) into 4×64, partial-reduced
/// to < 2²⁵¹ via 5·2²⁴⁸ ≡ 1. Input integer is Σ bᵢ·2⁵²ⁱ < 2²⁶⁰.
#[inline(always)]
fn join_and_reduce(b0: u64, mut b1: u64, mut b2: u64, mut b3: u64, b4: u64) -> [u64; 4] {
    // h = bits ≥248 (limb-4 bits ≥40). h < 2¹²; ⌊h/5⌋ < 820 (10 bits).
    let h = b4 >> 40;
    let quo = (h * 0xCCCD) >> 18;
    let rem = h - 5 * quo;
    let l0 = b0 + quo; // < 2⁵² + 2¹⁰
    b1 += l0 >> 52;
    b2 += b1 >> 52;
    b3 += b2 >> 52;
    // b3>>52 (≤ 1) carries into limb 4's place (bit 208).
    let l4 = (b4 & 0xFF_FFFF_FFFF) + (rem << 40) + (b3 >> 52); // < 6·2⁴⁰ < 2⁴³
    let w = [
        (l0 & MASK52) | (b1 & MASK52) << 52,
        (b1 & MASK52) >> 12 | (b2 & MASK52) << 40,
        (b2 & MASK52) >> 24 | (b3 & MASK52) << 28,
        (b3 & MASK52) >> 36 | l4 << 16,
    ];
    debug_assert!(w[3] >> 59 == 0, "ifma8 unpack: still ≥ 2²⁵¹");
    w
}

/// `r[k] ← a[k]·b[k]·2⁻²⁵⁶ mod p` for k=0..8, where `a[k]`, `b[k]` are
/// radix-64 integers in [0, 2²⁵²). Each output is in [0, 2²⁵¹), i.e. a valid
/// `GF5_248` representative, so this composes with the production `Fp` ops
/// without further reduction.
#[inline(never)]
pub fn fp_mul8_u64x4(a: &[[u64; 4]; 8], b: &[[u64; 4]; 8]) -> [[u64; 4]; 8] {
    let r = fp_mul_ifma8(&pack_u64x4(a), &pack_u64x4_x16(b));
    unpack_u64x4(r)
}

// ---------------------------------------------------------------------------
// SoA-resident pipeline: stay in radix-52 across many muls.
//
// `fp_mul_ifma8` is Montgomery wrt R₅₂ = 2²⁶⁰; the asm `Fp` is wrt R₆₄ = 2²⁵⁶.
// `fp_mul8_r64(a,b) := 16·(a·b·R₅₂⁻¹) = a·b·R₆₄⁻¹` keeps SoA values in the
// R₆₄ domain so they compose freely and convert losslessly to `GF5_248` at
// the boundary.
// ---------------------------------------------------------------------------

/// `a·b·2⁻²⁵⁶ mod p` per lane. Inputs must have limbs < 2⁵²; output limbs
/// are *lazy* (limb 0 < 2⁵² + 2¹⁴, limb 4 < 2⁵² + 2⁴⁴, others < 2⁵²) and
/// must be normalised before being fed back into another `fp_mul_ifma8`.
#[inline]
pub fn fp_mul8_r64(a: &FpIfma8, b: &FpIfma8) -> FpIfma8 {
    let mut r = fp_mul_ifma8(a, b);
    unsafe {
        let m = _mm512_set1_epi64(MASK52 as i64);
        let mut prev_hi = _mm512_setzero_si512();
        for i in 0..5 {
            let lo = _mm512_and_si512(_mm512_slli_epi64::<4>(r.limbs[i]), m);
            let hi = _mm512_srli_epi64::<48>(r.limbs[i]);
            r.limbs[i] = _mm512_or_si512(lo, prev_hi);
            prev_hi = hi;
        }
        let fold_lo = _mm512_set1_epi64(TWO_260_MOD_P_LO as i64);
        let fold_l4 = _mm512_set1_epi64(TWO_260_MOD_P_LIMB4 as i64);
        r.limbs[0] = _mm512_madd52lo_epu64(r.limbs[0], prev_hi, fold_lo);
        r.limbs[4] = _mm512_madd52lo_epu64(r.limbs[4], prev_hi, fold_l4);
    }
    r
}

/// Signed-limb carry-propagate so every limb ends in [0, 2⁵²). Accepts
/// limbs in roughly [−2⁵⁶, 2⁵⁶] with non-negative integer value per lane.
#[inline(always)]
pub fn normalize_signed8(a: &mut FpIfma8) {
    unsafe {
        let m = _mm512_set1_epi64(MASK52 as i64);
        let fold_lo = _mm512_set1_epi64(TWO_260_MOD_P_LO as i64);
        let fold_l4 = _mm512_set1_epi64(TWO_260_MOD_P_LIMB4 as i64);
        // Pass 1: |c| ≲ 16; pass 2: |c| ≤ 1; pass 3 settles limb0/4 (since
        // the integer value is ≥ 0 the third fold's c = 0).
        for _ in 0..3 {
            let mut c = _mm512_setzero_si512();
            for i in 0..5 {
                let t = _mm512_add_epi64(a.limbs[i], c);
                a.limbs[i] = _mm512_and_si512(t, m);
                c = _mm512_srai_epi64::<52>(t);
            }
            a.limbs[0] = _mm512_add_epi64(a.limbs[0], _mm512_mullo_epi64(c, fold_lo));
            a.limbs[4] = _mm512_add_epi64(a.limbs[4], _mm512_mullo_epi64(c, fold_l4));
        }
        let t = a.limbs[0];
        a.limbs[0] = _mm512_and_si512(t, m);
        a.limbs[1] = _mm512_add_epi64(a.limbs[1], _mm512_srli_epi64::<52>(t));
        #[cfg(debug_assertions)]
        for i in 0..5 {
            let mut s = [0u64; 8];
            _mm512_storeu_si512(s.as_mut_ptr().cast(), a.limbs[i]);
            for &v in &s {
                assert!(v < (1u64 << 52), "normalize_signed8: limb {i} = {v:#x}");
            }
        }
    }
}

/// Permute lanes by a fixed index table.
#[inline(always)]
pub fn permute8(a: &FpIfma8, idx: &[u64; 8]) -> FpIfma8 {
    unsafe {
        let iv = _mm512_loadu_si512(idx.as_ptr().cast());
        FpIfma8 {
            limbs: a.limbs.map(|l| _mm512_permutexvar_epi64(iv, l)),
        }
    }
}

/// Pack eight `GF5_248`-layout limbs (4×64 each) into SoA radix-2⁵².
#[inline(always)]
pub fn pack8(e: &[[u64; 4]; 8]) -> FpIfma8 {
    pack_u64x4(e)
}

/// Inverse of `pack8` (with partial reduce to < 2²⁵¹).
#[inline(always)]
pub fn unpack8(a: FpIfma8) -> [[u64; 4]; 8] {
    unpack_u64x4(a)
}

/// As `fp_mul8_u64x4` but with a single shared SoA operand `a` against two
/// different `b`'s — saves one pack when computing `a∘b₁` and `a∘b₂`.
#[inline]
pub fn fp_mul8_u64x4_shared_a(
    a: &[[u64; 4]; 8],
    b1: &[[u64; 4]; 8],
    b2: &[[u64; 4]; 8],
) -> ([[u64; 4]; 8], [[u64; 4]; 8]) {
    let pa = pack_u64x4(a);
    let r1 = fp_mul_ifma8(&pa, &pack_u64x4_x16(b1));
    let r2 = fp_mul_ifma8(&pa, &pack_u64x4_x16(b2));
    (unpack_u64x4(r1), unpack_u64x4(r2))
}

impl FpIfma8 {
    /// Build from 8 single elements (test/bench only; the real win needs callers
    /// to keep data in SoA form across operations).
    pub fn from_singles(s: &[FpIfma; 8]) -> Self {
        let mut buf = [[0u64; 8]; 5];
        for k in 0..8 {
            for i in 0..5 {
                buf[i][k] = s[k].0[i];
            }
        }
        let mut r = Self::default();
        for i in 0..5 {
            r.limbs[i] = unsafe { _mm512_loadu_si512(buf[i].as_ptr().cast()) };
        }
        r
    }
    pub fn to_singles(self) -> [FpIfma; 8] {
        let mut buf = [[0u64; 8]; 5];
        for i in 0..5 {
            unsafe { _mm512_storeu_si512(buf[i].as_mut_ptr().cast(), self.limbs[i]) };
        }
        let mut out = [FpIfma::ZERO; 8];
        for k in 0..8 {
            for i in 0..5 {
                out[k].0[i] = buf[i][k];
            }
        }
        out
    }
}

/// `r[k] ← a[k]·b[k]·R⁻¹ mod p` for k=0..8. Inputs must have limbs < 2⁵².
#[inline]
pub fn fp_mul_ifma8(a: &FpIfma8, b: &FpIfma8) -> FpIfma8 {
    unsafe {
        let zero = _mm512_setzero_si512();
        let mask = _mm512_set1_epi64(MASK52 as i64);
        let ptop = _mm512_set1_epi64(P_TOP_52 as i64);
        let fold_lo = _mm512_set1_epi64(TWO_260_MOD_P_LO as i64);
        let fold_l4 = _mm512_set1_epi64(TWO_260_MOD_P_LIMB4 as i64);
        let al = a.limbs;
        let bl = b.limbs;

        // Schoolbook 5×5: z[k] = Σ_{i+j=k} lo(a_i b_j) + Σ_{i+j=k-1} hi(a_i b_j).
        // 25 lo + 25 hi = 50 madd52, no permutes.
        let mut z = [zero; 10];
        macro_rules! madd {
            ($dst:expr, lo, $i:expr, $j:expr) => {
                $dst = _mm512_madd52lo_epu64($dst, al[$i], bl[$j]);
            };
            ($dst:expr, hi, $i:expr, $j:expr) => {
                $dst = _mm512_madd52hi_epu64($dst, al[$i], bl[$j]);
            };
        }
        madd!(z[0], lo, 0, 0);
        madd!(z[1], lo, 0, 1);
        madd!(z[1], lo, 1, 0);
        madd!(z[1], hi, 0, 0);
        madd!(z[2], lo, 0, 2);
        madd!(z[2], lo, 1, 1);
        madd!(z[2], lo, 2, 0);
        madd!(z[2], hi, 0, 1);
        madd!(z[2], hi, 1, 0);
        madd!(z[3], lo, 0, 3);
        madd!(z[3], lo, 1, 2);
        madd!(z[3], lo, 2, 1);
        madd!(z[3], lo, 3, 0);
        madd!(z[3], hi, 0, 2);
        madd!(z[3], hi, 1, 1);
        madd!(z[3], hi, 2, 0);
        madd!(z[4], lo, 0, 4);
        madd!(z[4], lo, 1, 3);
        madd!(z[4], lo, 2, 2);
        madd!(z[4], lo, 3, 1);
        madd!(z[4], lo, 4, 0);
        madd!(z[4], hi, 0, 3);
        madd!(z[4], hi, 1, 2);
        madd!(z[4], hi, 2, 1);
        madd!(z[4], hi, 3, 0);
        madd!(z[5], lo, 1, 4);
        madd!(z[5], lo, 2, 3);
        madd!(z[5], lo, 3, 2);
        madd!(z[5], lo, 4, 1);
        madd!(z[5], hi, 0, 4);
        madd!(z[5], hi, 1, 3);
        madd!(z[5], hi, 2, 2);
        madd!(z[5], hi, 3, 1);
        madd!(z[5], hi, 4, 0);
        madd!(z[6], lo, 2, 4);
        madd!(z[6], lo, 3, 3);
        madd!(z[6], lo, 4, 2);
        madd!(z[6], hi, 1, 4);
        madd!(z[6], hi, 2, 3);
        madd!(z[6], hi, 3, 2);
        madd!(z[6], hi, 4, 1);
        madd!(z[7], lo, 3, 4);
        madd!(z[7], lo, 4, 3);
        madd!(z[7], hi, 2, 4);
        madd!(z[7], hi, 3, 3);
        madd!(z[7], hi, 4, 2);
        madd!(z[8], lo, 4, 4);
        madd!(z[8], hi, 3, 4);
        madd!(z[8], hi, 4, 3);
        madd!(z[9], hi, 4, 4);

        // Montgomery reduction. v[0..3] depend only on z[0..3]+carries (no folds),
        // so extract them first, then fold v[0..3]·P_TOP into z[4..8] in parallel,
        // then resume the carry chain for v[4], c[0..3], c[4].
        let mut carry = zero;
        let mut v = [zero; 5];
        for k in 0..4 {
            let t = _mm512_add_epi64(z[k], carry);
            v[k] = _mm512_and_si512(t, mask);
            carry = _mm512_srli_epi64::<52>(t);
        }
        // Parallel fold of v[0..3]·P_TOP into z[4..8] (lo to k+4, hi to k+5).
        z[4] = _mm512_madd52lo_epu64(z[4], v[0], ptop);
        z[5] = _mm512_madd52hi_epu64(z[5], v[0], ptop);
        z[5] = _mm512_madd52lo_epu64(z[5], v[1], ptop);
        z[6] = _mm512_madd52hi_epu64(z[6], v[1], ptop);
        z[6] = _mm512_madd52lo_epu64(z[6], v[2], ptop);
        z[7] = _mm512_madd52hi_epu64(z[7], v[2], ptop);
        z[7] = _mm512_madd52lo_epu64(z[7], v[3], ptop);
        z[8] = _mm512_madd52hi_epu64(z[8], v[3], ptop);
        // Resume carry chain at column 4.
        let t = _mm512_add_epi64(z[4], carry);
        v[4] = _mm512_and_si512(t, mask);
        carry = _mm512_srli_epi64::<52>(t);
        z[8] = _mm512_madd52lo_epu64(z[8], v[4], ptop);
        z[9] = _mm512_madd52hi_epu64(z[9], v[4], ptop);
        let mut c = [zero; 5];
        for k in 5..9 {
            let t = _mm512_add_epi64(z[k], carry);
            c[k - 5] = _mm512_and_si512(t, mask);
            carry = _mm512_srli_epi64::<52>(t);
        }
        c[4] = _mm512_add_epi64(z[9], carry);

        // Normalize: c[4] may exceed 2⁵²; fold via 2²⁶⁰ ≡ 2²⁴⁸+819.
        let mut carry = zero;
        for i in 0..5 {
            let t = _mm512_add_epi64(c[i], carry);
            c[i] = _mm512_and_si512(t, mask);
            carry = _mm512_srli_epi64::<52>(t);
        }
        c[0] = _mm512_madd52lo_epu64(c[0], carry, fold_lo);
        c[4] = _mm512_madd52lo_epu64(c[4], carry, fold_l4);
        // One more limb-0 carry (carry ≤ ~8, 819·8 = 6552 < 2¹³, always fits).
        let t = c[0];
        c[0] = _mm512_and_si512(t, mask);
        c[1] = _mm512_add_epi64(c[1], _mm512_srli_epi64::<52>(t));

        FpIfma8 { limbs: c }
    }
}

// ---------------------------------------------------------------------------
// R² mod p in radix-2⁵² (computed offline; see test `r2_const`).
// R = 2²⁶⁰; R² mod p = 2⁵²⁰ mod p.
// ---------------------------------------------------------------------------

pub const R2_52: FpIfma = FpIfma([
    0x33333333d70a3,
    0x3333333333333,
    0x3333333333333,
    0x3333333333333,
    0x43333333333,
    0,
    0,
    0,
]);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf::{Fp, FpInner};
    use crate::test_util::Prng;

    #[test]
    fn r2_const() {
        // Recompute R² mod p and check the hardcoded constant.
        // p = 5·2²⁴⁸ − 1; R = 2²⁶⁰. Compute via 520 doublings of 1 mod p.
        // Do it in u128 chunks with the existing 4-limb Fp for reference.
        let one = Fp::from_small(1);
        let mut r = one;
        for _ in 0..520 {
            r = r.dbl();
        }
        let bytes = r.encode();
        let want = FpIfma::from_le_bytes_raw(&bytes);
        assert_eq!(&R2_52.0[..5], &want.0[..5], "R2_52 constant mismatch");
    }

    #[test]
    fn ref52_matches_asm() {
        let mut prng = Prng(0x51FA);
        for _ in 0..10_000 {
            let a = prng.fp();
            let b = prng.fp();
            let want = (a * b).encode();
            let ai = FpIfma::from_le_bytes_montgomery(&a.encode());
            let bi = FpIfma::from_le_bytes_montgomery(&b.encode());
            let r = fp_mul_ref52(&ai, &bi);
            assert_eq!(r.encode(), want, "scalar radix-52 mul mismatch");
        }
    }

    #[test]
    fn ifma_matches_asm() {
        let mut prng = Prng(0x1F4A);
        // Edge cases first.
        let edges: [Fp; 4] = [Fp::ZERO, Fp::ONE, Fp::MINUS_ONE, Fp::from_small(u64::MAX)];
        for &a in &edges {
            for &b in &edges {
                check(a, b);
            }
        }
        for _ in 0..100_000 {
            check(prng.fp(), prng.fp());
        }
        fn check(a: Fp, b: Fp) {
            let want = (a * b).encode();
            let ai = FpIfma::from_le_bytes_montgomery(&a.encode());
            let bi = FpIfma::from_le_bytes_montgomery(&b.encode());
            let r = fp_mul_ifma(&ai, &bi);
            assert_eq!(r.encode(), want, "IFMA mul mismatch a={a:?} b={b:?}");
        }
    }

    #[test]
    fn mul8_u64x4_matches_fp() {
        let mut prng = Prng(0x764);
        let edges = [Fp::ZERO, Fp::ONE, Fp::MINUS_ONE, Fp::from_small(u64::MAX)];
        for _ in 0..50_000 {
            let mut a = [Fp::ZERO; 8];
            let mut b = [Fp::ZERO; 8];
            for k in 0..8 {
                a[k] = if k < 4 { edges[k] } else { prng.fp() };
                b[k] = if k < 4 { edges[3 - k] } else { prng.fp() };
            }
            let aw = a.map(|x| x.0 .0);
            let bw = b.map(|x| x.0 .0);
            let r = fp_mul8_u64x4(&aw, &bw);
            for k in 0..8 {
                assert_eq!(Fp(FpInner(r[k])), a[k] * b[k], "lane {k}");
            }
        }
    }

    #[test]
    fn mul8_u64x4_noreduce_inputs() {
        // Inputs from add_noreduce / sub_2p_noreduce can be up to ~2²⁵².
        let mut prng = Prng(0x765);
        for _ in 0..20_000 {
            let mut aw = [[0u64; 4]; 8];
            let mut bw = [[0u64; 4]; 8];
            let mut a = [Fp::ZERO; 8];
            let mut b = [Fp::ZERO; 8];
            for k in 0..8 {
                let (x, y, z, w) = (prng.fp(), prng.fp(), prng.fp(), prng.fp());
                let an = Fp::add_noreduce(&x, &y);
                let bn = Fp::sub_2p_noreduce(&z, &w);
                aw[k] = an.0 .0;
                bw[k] = bn.0 .0;
                a[k] = x + y;
                b[k] = z - w;
            }
            let r = fp_mul8_u64x4(&aw, &bw);
            for k in 0..8 {
                assert_eq!(Fp(FpInner(r[k])), a[k] * b[k], "lane {k}");
            }
        }
    }

    #[test]
    fn ifma8_matches_asm() {
        let mut prng = Prng(0xB8);
        for _ in 0..10_000 {
            let mut a = [Fp::ZERO; 8];
            let mut b = [Fp::ZERO; 8];
            for k in 0..8 {
                a[k] = prng.fp();
                b[k] = prng.fp();
            }
            let ai =
                FpIfma8::from_singles(&a.map(|x| FpIfma::from_le_bytes_montgomery(&x.encode())));
            let bi =
                FpIfma8::from_singles(&b.map(|x| FpIfma::from_le_bytes_montgomery(&x.encode())));
            let ri = fp_mul_ifma8(&ai, &bi).to_singles();
            for k in 0..8 {
                let want = (a[k] * b[k]).encode();
                assert_eq!(ri[k].encode(), want, "lane {k}");
            }
        }
    }

    #[test]
    fn ifma_matches_ref52() {
        let mut prng = Prng(0xC0FFEE);
        for _ in 0..100_000 {
            let mut a = FpIfma::ZERO;
            let mut b = FpIfma::ZERO;
            for i in 0..5 {
                a.0[i] = prng.next() & MASK52;
                b.0[i] = prng.next() & MASK52;
            }
            assert_eq!(fp_mul_ifma(&a, &b).0, fp_mul_ref52(&a, &b).0);
        }
    }
}
