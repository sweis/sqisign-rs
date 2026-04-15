//! Vendored from <https://github.com/GiacomoPope/gf5_248_experiments>
//! (4-limb Montgomery GF(p), p = 5·2²⁴⁸ - 1, with x86_64 BMI2/ADX asm).
//!
//! Local modifications:
//!   - utils64 vendored alongside; fp2-crate trait impls and rand_core removed.
//!   - asm gated on `cfg(gf5_248_asm)` (set by build.rs).
//!   - `set_mul`/`set_square` use inline `asm!()` (`inline_asm.rs`) so they
//!     can be inlined by LLVM; the extern-C `gf5_248.S` versions are kept for
//!     reference and for the `*_via_extern` test/bench oracles.
//!   - `set_add`/`set_sub` always use the portable intrinsics path (call
//!     overhead on the asm versions outweighs the ~10-instruction body).
#![allow(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(unused)]
#![allow(unsafe_code)]
#![allow(missing_debug_implementations)]

mod utils64;
use utils64::{addcarry_u64, lzcnt, sgnw, subborrow_u64, umull, umull_add, umull_x2, umull_x2_add};

#[cfg(gf5_248_asm)]
mod inline_asm;

use core::convert::TryFrom;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[cfg(gf5_248_asm)]
unsafe extern "C" {
    unsafe fn fp_add_asm(dst: *mut GF5_248, a: *const GF5_248, b: *const GF5_248);
    unsafe fn fp_sub_asm(dst: *mut GF5_248, a: *const GF5_248, b: *const GF5_248);
    unsafe fn fp_mul_asm(dst: *mut GF5_248, a: *const GF5_248, b: *const GF5_248);
    unsafe fn fp_sqr_asm(dst: *mut GF5_248, a: *const GF5_248);
    unsafe fn fp_sop_asm(dst: *mut GF5_248, a: *const GF5_248, b: *const GF5_248);
    unsafe fn fp_dop_asm(dst: *mut GF5_248, a: *const GF5_248, b: *const GF5_248);
    unsafe fn fp2_sqr_re_asm(dst: *mut GF5_248, a: *const GF5_248);
    unsafe fn fp2_sqr_im_asm(dst: *mut GF5_248, a: *const GF5_248);
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct GF5_248([u64; 4]);

impl zeroize::DefaultIsZeroes for GF5_248 {}

impl GF5_248 {
    // IMPLEMENTATION NOTES
    // ====================
    //
    // Modulus is: q = 5*2^248 - 1
    // We represent the value x (integer in Z_q) with the integer y,
    // with the following rules:
    //    y mod q = x*2^256 mod q
    //    0 <= y < 2^251

    // Element encoding length (in bytes); always 32 bytes.
    pub const ENCODED_LENGTH: usize = 32;
    pub const N: usize = 4;
    pub const BIT_LENGTH: usize = 251;
    pub const SUM_OF_PRODUCTS_ADDITIONAL_SUB: bool = false;

    // Modulus q in base 2^64 (low-to-high order).
    pub const MODULUS: [u64; 4] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x04FFFFFFFFFFFFFF,
    ];

    // 2*q
    const MOD_X2: [u64; 4] = [
        0xFFFFFFFFFFFFFFFE,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x09FFFFFFFFFFFFFF,
    ];

    pub const ZERO: Self = Self([0, 0, 0, 0]);
    pub const ONE: Self = Self([
        0x0000000000000033,
        0x0000000000000000,
        0x0000000000000000,
        0x0100000000000000,
    ]);
    pub const TWO: Self = Self([
        0x0000000000000066,
        0x0000000000000000,
        0x0000000000000000,
        0x0200000000000000,
    ]);
    pub const THREE: Self = Self([
        0x0000000000000099,
        0x0000000000000000,
        0x0000000000000000,
        0x0300000000000000,
    ]);
    pub const FOUR: Self = Self([
        0x00000000000000cc,
        0x0000000000000000,
        0x0000000000000000,
        0x0400000000000000,
    ]);
    pub const MINUS_ONE: Self = Self([
        0xFFFFFFFFFFFFFFCC,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x03FFFFFFFFFFFFFF,
    ]);

    // 1/2^244 in the field, in Montgomery representation.
    const INVT244: Self = Self([
        0x0000000000001000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ]);

    // Montgomery representation of 2^256 (i.e. 2^512 mod q).
    const R2: Self = Self([
        0x3333333333333D70,
        0x3333333333333333,
        0x3333333333333333,
        0x0333333333333333,
    ]);

    /// Get a new instance containing the provided 256-bit integer,
    /// which is implicitly reduced modulo q. All 256-bit values are
    /// accepted. This function is meant for constant (compile-time)
    /// evaluation; at runtime, consider using `from_w64le()`, which
    /// is faster.
    /// This function expects the four 64-bit limbs in little-endian order.
    pub const fn w64le(x0: u64, x1: u64, x2: u64, x3: u64) -> Self {
        Self::const_mmul(Self([x0, x1, x2, x3]), Self::R2)
    }

    /// Get a new instance containing the provided 256-bit integer,
    /// which is implicitly reduced modulo q. All 256-bit values are
    /// accepted. This function is meant for constant (compile-time)
    /// evaluation; at runtime, consider using `from_w64be()`, which
    /// is faster.
    /// This function expects the four 64-bit limbs in big-endian order.
    pub const fn w64be(x3: u64, x2: u64, x1: u64, x0: u64) -> Self {
        Self::const_mmul(Self([x0, x1, x2, x3]), Self::R2)
    }

    /// Get a new instance containing the provided 256-bit integer,
    /// which is implicitly reduced modulo q. All 256-bit values are
    /// accepted.
    #[inline(always)]
    pub fn from_w64le(x0: u64, x1: u64, x2: u64, x3: u64) -> Self {
        let mut r = Self([x0, x1, x2, x3]);
        r.set_partial_reduce();
        r.set_mul(&Self::R2);
        r
    }

    /// Get a new instance containing the provided 256-bit integer,
    /// which is implicitly reduced modulo q. All 256-bit values are
    /// accepted.
    #[inline(always)]
    pub fn from_w64be(x3: u64, x2: u64, x1: u64, x0: u64) -> Self {
        let mut r = Self([x0, x1, x2, x3]);
        r.set_partial_reduce();
        r.set_mul(&Self::R2);
        r
    }

    // Note: extern-C asm for add/sub is *slower* than the portable
    // intrinsics version when composed into larger algorithms (call
    // overhead defeats LLVM scheduling). Always use portable add/sub.
    #[inline(always)]
    fn set_add(&mut self, rhs: &Self) {
        // Raw addition.
        let (d0, cc) = addcarry_u64(self.0[0], rhs.0[0], 0);
        let (d1, cc) = addcarry_u64(self.0[1], rhs.0[1], cc);
        let (d2, cc) = addcarry_u64(self.0[2], rhs.0[2], cc);
        let (d3, _) = addcarry_u64(self.0[3], rhs.0[3], cc);

        // Inputs were up to 2^251 - 1; sum can be up to 2^252 - 2.
        // We subtract q if the value is not lower than 2^251. Subtraction
        // of q is done by adding -q (modulo 2^256) (because -q has two
        // limbs equal to zero, which saves a few operations).
        let f = d3 >> 59;
        let (d0, cc) = addcarry_u64(d0, f, 0);
        let (d1, cc) = addcarry_u64(d1, 0, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, (0xFBu64 << 56) & f.wrapping_neg(), cc);

        // If the value is still not lower than 2^251 then we subtract q again.
        let f = d3 >> 59;
        let (d0, cc) = addcarry_u64(d0, f, 0);
        let (d1, cc) = addcarry_u64(d1, 0, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, (0xFBu64 << 56) & f.wrapping_neg(), cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    #[inline(always)]
    fn set_sub(&mut self, rhs: &Self) {
        // Raw subtraction.
        let (d0, cc) = subborrow_u64(self.0[0], rhs.0[0], 0);
        let (d1, cc) = subborrow_u64(self.0[1], rhs.0[1], cc);
        let (d2, cc) = subborrow_u64(self.0[2], rhs.0[2], cc);
        let (d3, cc) = subborrow_u64(self.0[3], rhs.0[3], cc);

        // If the result is negative, we add 2*q (by subtracting -2*q).
        let (m, _) = subborrow_u64(0, 0, cc);
        let (d0, cc) = subborrow_u64(d0, m & 2, 0);
        let (d1, cc) = subborrow_u64(d1, 0, cc);
        let (d2, cc) = subborrow_u64(d2, 0, cc);
        let (d3, _) = subborrow_u64(d3, (0xF6u64 << 56) & m, cc);

        // If the value is not lower than 2^251 then we subtract q (i.e.
        // we add -q).
        let f = d3 >> 59;
        let (d0, cc) = addcarry_u64(d0, f, 0);
        let (d1, cc) = addcarry_u64(d1, 0, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, (0xFBu64 << 56) & f.wrapping_neg(), cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    // Negate this value (in place).
    #[inline(always)]
    pub fn set_neg(&mut self) {
        // 2*q - self.
        let (d0, cc) = subborrow_u64(Self::MOD_X2[0], self.0[0], 0);
        let (d1, cc) = subborrow_u64(Self::MOD_X2[1], self.0[1], cc);
        let (d2, cc) = subborrow_u64(Self::MOD_X2[2], self.0[2], cc);
        let (d3, _) = subborrow_u64(Self::MOD_X2[3], self.0[3], cc);

        // If the value is not lower than 2^251 then we subtract q (i.e.
        // we add -q).
        let f = d3 >> 59;
        let (d0, cc) = addcarry_u64(d0, f, 0);
        let (d1, cc) = addcarry_u64(d1, 0, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, (0xFBu64 << 56) & f.wrapping_neg(), cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    // Conditionally copy the provided value ('a') into self:
    //  - If ctl == 0xFFFFFFFF, then the value of 'a' is copied into self.
    //  - If ctl == 0, then the value of self is unchanged.
    // ctl MUST be equal to 0 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn set_cond(&mut self, a: &Self, ctl: u32) {
        let cw = ((ctl as i32) as i64) as u64;
        self.0[0] ^= cw & (self.0[0] ^ a.0[0]);
        self.0[1] ^= cw & (self.0[1] ^ a.0[1]);
        self.0[2] ^= cw & (self.0[2] ^ a.0[2]);
        self.0[3] ^= cw & (self.0[3] ^ a.0[3]);
    }

    /// Negate this value if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline]
    pub fn set_cond_neg(&mut self, ctl: u32) {
        let v = -(self as &Self);
        self.set_cond(&v, ctl);
    }

    /// Set this value to either a or b, depending on whether the control
    /// word ctl is 0x00000000 or 0xFFFFFFFF, respectively.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline]
    pub fn set_select(&mut self, a: &Self, b: &Self, ctl: u32) {
        let c = ((ctl as i32) as i64) as u64;
        self.0[0] = a.0[0] ^ (c & (a.0[0] ^ b.0[0]));
        self.0[1] = a.0[1] ^ (c & (a.0[1] ^ b.0[1]));
        self.0[2] = a.0[2] ^ (c & (a.0[2] ^ b.0[2]));
        self.0[3] = a.0[3] ^ (c & (a.0[3] ^ b.0[3]));
    }

    // Return a value equal to either a0 (if ctl == 0) or a1 (if
    // ctl == 0xFFFFFFFF). Value ctl MUST be either 0 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn select(a0: &Self, a1: &Self, ctl: u32) -> Self {
        let mut r = *a0;
        r.set_cond(a1, ctl);
        r
    }

    // Conditionally swap two elements: values a and b are exchanged if
    // ctl == 0xFFFFFFFF, or not exchanged if ctl == 0x00000000. Value
    // ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn cond_swap(a: &mut Self, b: &mut Self, ctl: u32) {
        let cw = ((ctl as i32) as i64) as u64;
        let t = cw & (a.0[0] ^ b.0[0]);
        a.0[0] ^= t;
        b.0[0] ^= t;
        let t = cw & (a.0[1] ^ b.0[1]);
        a.0[1] ^= t;
        b.0[1] ^= t;
        let t = cw & (a.0[2] ^ b.0[2]);
        a.0[2] ^= t;
        b.0[2] ^= t;
        let t = cw & (a.0[3] ^ b.0[3]);
        a.0[3] ^= t;
        b.0[3] ^= t;
    }

    #[inline(always)]
    pub fn set_half(&mut self) {
        // Right-shift the value by 1 bit.
        let d0 = (self.0[0] >> 1) | (self.0[1] << 63);
        let d1 = (self.0[1] >> 1) | (self.0[2] << 63);
        let d2 = (self.0[2] >> 1) | (self.0[3] << 63);
        let d3 = self.0[3] >> 1;

        // If the dropped bit was 1, then we got (x - 1)/2 and we need
        // to add 1/2 mod q = (q + 1)/2 = 5*2^247.
        let d3 = d3 + ((self.0[0] & 1).wrapping_neg() & (5u64 << 55));

        // Value is necessarily in range: max value is
        // (2^251 - 1)/2 + 5*2^247, which is lower than 2^251.

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    #[inline(always)]
    pub fn half(self) -> Self {
        let mut r = self;
        r.set_half();
        r
    }

    // Multiply this value by 2.
    #[inline(always)]
    pub fn set_mul2(&mut self) {
        let r = *self;
        self.set_add(&r);
    }

    #[inline(always)]
    pub fn mul2(self) -> Self {
        let mut r = self;
        r.set_mul2();
        r
    }

    // Multiply this value by 3.
    #[inline(always)]
    pub fn set_mul3(&mut self) {
        // Tripling as an integer.
        let a0 = self.0[0];
        let a1 = self.0[1];
        let a2 = self.0[2];
        let a3 = self.0[3];
        let b0 = a0 << 1;
        let b1 = (a1 << 1) | (a0 >> 63);
        let b2 = (a2 << 1) | (a1 >> 63);
        let b3 = (a3 << 1) | (a2 >> 63);
        let (d0, cc) = addcarry_u64(a0, b0, 0);
        let (d1, cc) = addcarry_u64(a1, b1, cc);
        let (d2, cc) = addcarry_u64(a2, b2, cc);
        let (d3, _) = addcarry_u64(a3, b3, cc);
        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
        self.set_partial_reduce();
    }

    #[inline(always)]
    pub fn mul3(self) -> Self {
        let mut r = self;
        r.set_mul3();
        r
    }

    // Multiply this value by 4.
    #[inline(always)]
    pub fn set_mul4(&mut self) {
        self.0[3] = (self.0[3] << 2) | (self.0[2] >> 62);
        self.0[2] = (self.0[2] << 2) | (self.0[1] >> 62);
        self.0[1] = (self.0[1] << 2) | (self.0[0] >> 62);
        self.0[0] <<= 2;
        self.set_partial_reduce();
    }

    #[inline(always)]
    pub fn mul4(self) -> Self {
        let mut r = self;
        r.set_mul4();
        r
    }

    // Multiply this value by 8.
    #[inline(always)]
    pub fn set_mul8(&mut self) {
        self.0[3] = (self.0[3] << 3) | (self.0[2] >> 61);
        self.0[2] = (self.0[2] << 3) | (self.0[1] >> 61);
        self.0[1] = (self.0[1] << 3) | (self.0[0] >> 61);
        self.0[0] <<= 3;
        self.set_partial_reduce();
    }

    #[inline(always)]
    pub fn mul8(self) -> Self {
        let mut r = self;
        r.set_mul8();
        r
    }

    // Multiply this value by 16.
    #[inline(always)]
    pub fn set_mul16(&mut self) {
        self.0[3] = (self.0[3] << 4) | (self.0[2] >> 60);
        self.0[2] = (self.0[2] << 4) | (self.0[1] >> 60);
        self.0[1] = (self.0[1] << 4) | (self.0[0] >> 60);
        self.0[0] <<= 4;
        self.set_partial_reduce();
    }

    #[inline(always)]
    pub fn mul16(self) -> Self {
        let mut r = self;
        r.set_mul16();
        r
    }

    // Multiply this value by 32.
    #[inline(always)]
    pub fn set_mul32(&mut self) {
        self.0[3] = (self.0[3] << 5) | (self.0[2] >> 59);
        self.0[2] = (self.0[2] << 5) | (self.0[1] >> 59);
        self.0[1] = (self.0[1] << 5) | (self.0[0] >> 59);
        self.0[0] <<= 5;
        self.set_partial_reduce();
    }

    #[inline(always)]
    pub fn mul32(self) -> Self {
        let mut r = self;
        r.set_mul32();
        r
    }

    // Multiply this value by a small integer.
    #[inline(always)]
    pub fn set_mul_small(&mut self, k: i32) {
        // Get the absolute value of the multiplier (but remember the sign).
        let sign = (k >> 31) as u32;
        let b = ((k as u32) ^ sign).wrapping_sub(sign) as u64;

        // Compute the product as an integer over five words.
        // Max value is (2^32 - 1)*(2^251 - 1), so the top word (d4) is
        // at most 2^27 - 2.
        let (d0, d1) = umull(self.0[0], b);
        let (d2, d3) = umull(self.0[2], b);
        let (lo, hi) = umull(self.0[1], b);
        let (d1, cc) = addcarry_u64(d1, lo, 0);
        let (d2, cc) = addcarry_u64(d2, hi, cc);
        let (lo, d4) = umull(self.0[3], b);
        let (d3, cc) = addcarry_u64(d3, lo, cc);
        let (d4, _) = addcarry_u64(d4, 0, cc);

        // Extract low 248-bit part, and the high word (35 bits).
        let h = (d4 << 8) | (d3 >> 56);
        let d3 = d3 & 0x00FFFFFFFFFFFFFF;

        // Fold value h by dividing it by 5. Since 5*2^248 = 1 mod q,
        // We add floor(h/5) + (h mod 5)*2^248 to the value.
        let (_, z) = umull(h, 0xCCCCCCCCCCCCCCCD);
        let quo = z >> 2;
        let rem = h - (5 * quo);
        let (d0, cc) = addcarry_u64(d0, quo, 0);
        let (d1, cc) = addcarry_u64(d1, 0, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, rem << 56, cc);

        // The value is at most 5*2^248 + 6871947672, which is in
        // range.
        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;

        // We computed self*|k|; we must adjust for the sign of k.
        self.set_cond_neg(sign);
    }

    #[inline(always)]
    pub fn mul_small(self, x: i32) -> Self {
        let mut r = self;
        r.set_mul_small(x);
        r
    }

    // Montgomery reduction, i.e. division by 2^256. Output is also
    // normalized in the 0 to q-1 range.
    #[inline(always)]
    fn set_montgomery_reduce(&mut self) {
        // Let m = -1/q mod 2^256 = 5*2^248 + 1.
        // For input x, we compute f = x*m mod 2^256, then
        // h = x + f*q, which is a multiple of 2^256. The output
        // is then h/2^256. Note that if x < 2^256, then:
        //   h <= 2^256 - 1 + (2^256 - 1)*q
        //   h <= q*2^256 + 2^256 - q - 1
        // Since h = 0 mod 2^256 and, this implies that h <= q*2^256.
        // The output h/2^256 is thus in the 0 to q range (inclusive).

        // Input x.
        let x0 = self.0[0];
        let x1 = self.0[1];
        let x2 = self.0[2];
        let x3 = self.0[3];

        // f = x*(-1/q) mod 2^256
        let f0 = x0;
        let f1 = x1;
        let f2 = x2;
        let f3 = x3.wrapping_add(x0.wrapping_mul(5) << 56);

        // g = f*q
        let (g3, hi) = umull(f0, 5u64 << 56);
        let (g4, hi) = umull_add(f1, 5u64 << 56, hi);
        let (g5, hi) = umull_add(f2, 5u64 << 56, hi);
        let (g6, g7) = umull_add(f3, 5u64 << 56, hi);
        let (g0, cc) = subborrow_u64(0, f0, 0);
        let (g1, cc) = subborrow_u64(0, f1, cc);
        let (g2, cc) = subborrow_u64(0, f2, cc);
        let (g3, cc) = subborrow_u64(g3, f3, cc);
        let (g4, cc) = subborrow_u64(g4, 0, cc);
        let (g5, cc) = subborrow_u64(g5, 0, cc);
        let (g6, cc) = subborrow_u64(g6, 0, cc);
        let (g7, _) = subborrow_u64(g7, 0, cc);

        // h = x + f*q
        // We drop the low 256 bits.
        let (_, cc) = addcarry_u64(g0, x0, 0);
        let (_, cc) = addcarry_u64(g1, x1, cc);
        let (_, cc) = addcarry_u64(g2, x2, cc);
        let (_, cc) = addcarry_u64(g3, x3, cc);
        let (d0, cc) = addcarry_u64(g4, 0, cc);
        let (d1, cc) = addcarry_u64(g5, 0, cc);
        let (d2, cc) = addcarry_u64(g6, 0, cc);
        let (d3, _) = addcarry_u64(g7, 0, cc);

        // h \in [0..q]
        // Normalize: if output is q, replace with zero.
        let t = d0 & d1 & d2 & (d3 ^ !Self::MODULUS[3]);
        let (_, cc) = addcarry_u64(t, 1, 0);
        let (w, _) = subborrow_u64(0, 0, cc);
        let w = !w;

        self.0[0] = d0 & w;
        self.0[1] = d1 & w;
        self.0[2] = d2 & w;
        self.0[3] = d3 & w;
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn set_mul(&mut self, rhs: &Self) {
        let a = self.0;
        inline_asm::fp_mul(&mut self.0, &a, &rhs.0);
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    fn set_mul(&mut self, rhs: &Self) {
        let (a0, a1, a2, a3) = (self.0[0], self.0[1], self.0[2], self.0[3]);
        let (b0, b1, b2, b3) = (rhs.0[0], rhs.0[1], rhs.0[2], rhs.0[3]);

        // Product -> 502 bits
        let (e0, e1) = umull(a0, b0);
        let (e2, e3) = umull(a1, b1);
        let (e4, e5) = umull(a2, b2);
        let (e6, e7) = umull(a3, b3);

        let (lo, hi) = umull(a0, b1);
        let (e1, cc) = addcarry_u64(e1, lo, 0);
        let (e2, cc) = addcarry_u64(e2, hi, cc);
        let (lo, hi) = umull(a0, b3);
        let (e3, cc) = addcarry_u64(e3, lo, cc);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (lo, hi) = umull(a2, b3);
        let (e5, cc) = addcarry_u64(e5, lo, cc);
        let (e6, cc) = addcarry_u64(e6, hi, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo, hi) = umull(a1, b0);
        let (e1, cc) = addcarry_u64(e1, lo, 0);
        let (e2, cc) = addcarry_u64(e2, hi, cc);
        let (lo, hi) = umull(a3, b0);
        let (e3, cc) = addcarry_u64(e3, lo, cc);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (lo, hi) = umull(a3, b2);
        let (e5, cc) = addcarry_u64(e5, lo, cc);
        let (e6, cc) = addcarry_u64(e6, hi, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo, hi) = umull(a0, b2);
        let (e2, cc) = addcarry_u64(e2, lo, 0);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi) = umull(a1, b3);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (e6, cc) = addcarry_u64(e6, 0, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo, hi) = umull(a2, b0);
        let (e2, cc) = addcarry_u64(e2, lo, 0);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi) = umull(a3, b1);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (e6, cc) = addcarry_u64(e6, 0, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo, hi) = umull(a1, b2);
        let (lo2, hi2) = umull(a2, b1);
        let (lo, cc) = addcarry_u64(lo, lo2, 0);
        let (hi, tt) = addcarry_u64(hi, hi2, cc);
        let (e3, cc) = addcarry_u64(e3, lo, 0);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (e5, cc) = addcarry_u64(e5, tt as u64, cc);
        let (e6, cc) = addcarry_u64(e6, 0, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        // Montgomery reduction.
        //
        // The low part is lo(e) = e0..e3 (256 bits).
        // Let m = -1/q mod 2^256; the Montgomery reduction is adding
        // (lo(e)*m mod 2^256)*q to the high part g = e4..e7 (246 bits).
        //
        // We have m = 5*2^248 + 1; the product f = lo(e)*m mod 2^256 is
        // obtained by simply added e0*5 (mod 2^8) to the high byte of e3.
        let f0 = e0;
        let f1 = e1;
        let f2 = e2;
        let f3 = e3.wrapping_add(e0.wrapping_mul(5) << 56);

        // Compute g = f*q.
        let (g3, hi) = umull(f0, 5u64 << 56);
        let (g4, hi) = umull_add(f1, 5u64 << 56, hi);
        let (g5, hi) = umull_add(f2, 5u64 << 56, hi);
        let (g6, g7) = umull_add(f3, 5u64 << 56, hi);
        let (g0, cc) = subborrow_u64(0, f0, 0);
        let (g1, cc) = subborrow_u64(0, f1, cc);
        let (g2, cc) = subborrow_u64(0, f2, cc);
        let (g3, cc) = subborrow_u64(g3, f3, cc);
        let (g4, cc) = subborrow_u64(g4, 0, cc);
        let (g5, cc) = subborrow_u64(g5, 0, cc);
        let (g6, cc) = subborrow_u64(g6, 0, cc);
        let (g7, _) = subborrow_u64(g7, 0, cc);

        // We add g = f*q to e0..e7. Since e0..e7 < 2^502, and f < 2^256,
        // we know that the result is not less than
        // 2^502 + 2^256*5*2^248 < 6*2^504; it is also a multiple of
        // 2^256. After dividing by 2^256, we get a value which is
        // less than 6*2^248, i.e. already in our proper range.
        let (_, cc) = addcarry_u64(e0, g0, 0);
        let (_, cc) = addcarry_u64(e1, g1, cc);
        let (_, cc) = addcarry_u64(e2, g2, cc);
        let (_, cc) = addcarry_u64(e3, g3, cc);
        let (d0, cc) = addcarry_u64(e4, g4, cc);
        let (d1, cc) = addcarry_u64(e5, g5, cc);
        let (d2, cc) = addcarry_u64(e6, g6, cc);
        let (d3, _) = addcarry_u64(e7, g7, cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn set_square(&mut self) {
        let a = self.0;
        inline_asm::fp_sqr(&mut self.0, &a);
    }

    /// `fp_mul(x, x)` — kept for benchmarking and as a property-test oracle
    /// for `fp_sqr`.
    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn set_square_via_mul(&mut self) {
        let a = self.0;
        inline_asm::fp_mul(&mut self.0, &a, &a);
    }

    /// Extern-C `.S` `fp_mul_asm` — test oracle for the inline `fp_mul`.
    #[cfg(gf5_248_asm)]
    pub fn set_mul_via_extern(&mut self, rhs: &Self) {
        unsafe { fp_mul_asm(self, self, rhs) }
    }

    /// Extern-C `.S` `fp_sqr_asm` — test oracle for the inline `fp_sqr`.
    #[cfg(gf5_248_asm)]
    pub fn set_square_via_extern(&mut self) {
        unsafe { fp_sqr_asm(self, self) }
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    pub fn set_square(&mut self) {
        let (a0, a1, a2, a3) = (self.0[0], self.0[1], self.0[2], self.0[3]);

        // 1. Non-square products. Max intermediate value:
        //   a0*a1            * 2^64
        //   a0*a2            * 2^128
        //   (a0*a3 + a1*a2)  * 2^192
        //   a1*a3            * 2^256
        //   a2*a3            * 2^320
        // for a total which is slightly below 2^448, which means that
        // the value fits on e1..e6 (no possible carry into e7).
        let (e1, e2) = umull(a0, a1);
        let (e3, e4) = umull(a0, a3);
        let (e5, e6) = umull(a2, a3);
        let (lo, hi) = umull(a0, a2);
        let (e2, cc) = addcarry_u64(e2, lo, 0);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi) = umull(a1, a3);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (e6, _) = addcarry_u64(e6, 0, cc);
        let (lo, hi) = umull(a1, a2);
        let (e3, cc) = addcarry_u64(e3, lo, 0);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (e5, cc) = addcarry_u64(e5, 0, cc);
        let (e6, _) = addcarry_u64(e6, 0, cc);

        // 2. Double the intermediate value, then add the squares.
        let e7 = e6 >> 63;
        let e6 = (e6 << 1) | (e5 >> 63);
        let e5 = (e5 << 1) | (e4 >> 63);
        let e4 = (e4 << 1) | (e3 >> 63);
        let e3 = (e3 << 1) | (e2 >> 63);
        let e2 = (e2 << 1) | (e1 >> 63);
        let e1 = e1 << 1;

        let (e0, hi) = umull(a0, a0);
        let (e1, cc) = addcarry_u64(e1, hi, 0);
        let (lo, hi) = umull(a1, a1);
        let (e2, cc) = addcarry_u64(e2, lo, cc);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi) = umull(a2, a2);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (lo, hi) = umull(a3, a3);
        let (e6, cc) = addcarry_u64(e6, lo, cc);
        let (e7, _) = addcarry_u64(e7, hi, cc);

        // 3. Reduction (see set_mul() for details).
        let f0 = e0;
        let f1 = e1;
        let f2 = e2;
        let f3 = e3.wrapping_add(e0.wrapping_mul(5) << 56);

        let (g3, hi) = umull(f0, 5u64 << 56);
        let (g4, hi) = umull_add(f1, 5u64 << 56, hi);
        let (g5, hi) = umull_add(f2, 5u64 << 56, hi);
        let (g6, g7) = umull_add(f3, 5u64 << 56, hi);
        let (g0, cc) = subborrow_u64(0, f0, 0);
        let (g1, cc) = subborrow_u64(0, f1, cc);
        let (g2, cc) = subborrow_u64(0, f2, cc);
        let (g3, cc) = subborrow_u64(g3, f3, cc);
        let (g4, cc) = subborrow_u64(g4, 0, cc);
        let (g5, cc) = subborrow_u64(g5, 0, cc);
        let (g6, cc) = subborrow_u64(g6, 0, cc);
        let (g7, _) = subborrow_u64(g7, 0, cc);

        let (_, cc) = addcarry_u64(e0, g0, 0);
        let (_, cc) = addcarry_u64(e1, g1, cc);
        let (_, cc) = addcarry_u64(e2, g2, cc);
        let (_, cc) = addcarry_u64(e3, g3, cc);
        let (d0, cc) = addcarry_u64(e4, g4, cc);
        let (d1, cc) = addcarry_u64(e5, g5, cc);
        let (d2, cc) = addcarry_u64(e6, g6, cc);
        let (d3, _) = addcarry_u64(e7, g7, cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    // Square this value.
    #[inline(always)]
    pub fn square(self) -> Self {
        let mut r = self;
        r.set_square();
        r
    }

    // Square this value n times (in place).
    #[inline(always)]
    fn set_n_square(&mut self, n: u32) {
        for _ in 0..n {
            self.set_square();
        }
    }

    // Square this value n times.
    #[inline(always)]
    pub fn n_square(self, n: u32) -> Self {
        let mut r = self;
        r.set_n_square(n);
        r
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    fn sum_of_products_u64(a1: u64, b1: u64, a2: u64, b2: u64, c: u8) -> (u64, u64, u8) {
        let (lo1, hi1) = umull(a1, b1);
        let (lo2, hi2) = umull(a2, b2);
        let (lo, cc) = addcarry_u64(lo1, lo2, c);
        let (hi, tt) = addcarry_u64(hi1, hi2, cc);

        (lo, hi, tt)
    }

    /// `a + b` as raw limbs (no conditional reduction). For inputs < 2²⁵¹,
    /// result < 2²⁵²; safe as a `set_mul`/`fp_mul` operand.
    #[inline(always)]
    pub fn add_noreduce(a: &Self, b: &Self) -> Self {
        let (d0, cc) = addcarry_u64(a.0[0], b.0[0], 0);
        let (d1, cc) = addcarry_u64(a.0[1], b.0[1], cc);
        let (d2, cc) = addcarry_u64(a.0[2], b.0[2], cc);
        let (d3, _) = addcarry_u64(a.0[3], b.0[3], cc);
        Self([d0, d1, d2, d3])
    }

    /// `a − b + 2p` as raw limbs (no conditional reduction). For inputs
    /// < 2²⁵¹, result is in (0, 2²⁵²+2²⁵¹); safe as a `set_mul` operand.
    #[inline(always)]
    pub fn sub_2p_noreduce(a: &Self, b: &Self) -> Self {
        let (d0, cc) = subborrow_u64(a.0[0], b.0[0], 0);
        let (d1, cc) = subborrow_u64(a.0[1], b.0[1], cc);
        let (d2, cc) = subborrow_u64(a.0[2], b.0[2], cc);
        let (d3, _) = subborrow_u64(a.0[3], b.0[3], cc);
        let (d0, cc) = addcarry_u64(d0, Self::MOD_X2[0], 0);
        let (d1, cc) = addcarry_u64(d1, Self::MOD_X2[1], cc);
        let (d2, cc) = addcarry_u64(d2, Self::MOD_X2[2], cc);
        let (d3, _) = addcarry_u64(d3, Self::MOD_X2[3], cc);
        Self([d0, d1, d2, d3])
    }

    /// 2p − x as raw limbs (no reduction). For x < 2²⁵¹, result < 2²⁵².
    #[inline(always)]
    fn raw_2p_minus(x: &[u64; 4]) -> [u64; 4] {
        let (d0, cc) = subborrow_u64(Self::MOD_X2[0], x[0], 0);
        let (d1, cc) = subborrow_u64(Self::MOD_X2[1], x[1], cc);
        let (d2, cc) = subborrow_u64(Self::MOD_X2[2], x[2], cc);
        let (d3, _) = subborrow_u64(Self::MOD_X2[3], x[3], cc);
        [d0, d1, d2, d3]
    }

    /// One conditional `−p` round: subtract `p` iff bit 251 is set. Used by
    /// `set_add`/`set_sub` and the butterfly variants.
    #[inline(always)]
    fn cond_sub_p(d0: u64, d1: u64, d2: u64, d3: u64) -> (u64, u64, u64, u64) {
        let f = d3 >> 59;
        let m = f.wrapping_neg();
        let (d0, cc) = subborrow_u64(d0, m, 0);
        let (d1, cc) = subborrow_u64(d1, m, cc);
        let (d2, cc) = subborrow_u64(d2, m, cc);
        let (d3, _) = subborrow_u64(d3, Self::MODULUS[3] & m, cc);
        (d0, d1, d2, d3)
    }

    /// Butterfly: `(a+b, a−b)`, both fully reduced. One body so the loads of
    /// `a`/`b` are shared; LLVM emits ADC chains directly so no asm needed.
    #[inline(always)]
    pub fn butterfly(a: &Self, b: &Self) -> (Self, Self) {
        let (a0, a1, a2, a3) = (a.0[0], a.0[1], a.0[2], a.0[3]);
        let (b0, b1, b2, b3) = (b.0[0], b.0[1], b.0[2], b.0[3]);
        // sum = a+b
        let (s0, cc) = addcarry_u64(a0, b0, 0);
        let (s1, cc) = addcarry_u64(a1, b1, cc);
        let (s2, cc) = addcarry_u64(a2, b2, cc);
        let (s3, _) = addcarry_u64(a3, b3, cc);
        let (s0, s1, s2, s3) = Self::cond_sub_p(s0, s1, s2, s3);
        let (s0, s1, s2, s3) = Self::cond_sub_p(s0, s1, s2, s3);
        // diff = a−b; cond-add 2p on borrow, then one cond-sub-p
        let (d0, cc) = subborrow_u64(a0, b0, 0);
        let (d1, cc) = subborrow_u64(a1, b1, cc);
        let (d2, cc) = subborrow_u64(a2, b2, cc);
        let (d3, cc) = subborrow_u64(a3, b3, cc);
        let m = (cc as u64).wrapping_neg();
        let (d0, cc) = addcarry_u64(d0, m & Self::MOD_X2[0], 0);
        let (d1, cc) = addcarry_u64(d1, m, cc);
        let (d2, cc) = addcarry_u64(d2, m, cc);
        let (d3, _) = addcarry_u64(d3, m & Self::MOD_X2[3], cc);
        let (d0, d1, d2, d3) = Self::cond_sub_p(d0, d1, d2, d3);
        (Self([s0, s1, s2, s3]), Self([d0, d1, d2, d3]))
    }

    /// Lazy butterfly: `(a+b, a−b)` with both outputs < 2²⁵² (one reduction
    /// round each). Safe as `set_mul`/`set_square` operands; **not** as
    /// `sum_of_products` operands.
    #[inline(always)]
    pub fn butterfly_lazy(a: &Self, b: &Self) -> (Self, Self) {
        let (a0, a1, a2, a3) = (a.0[0], a.0[1], a.0[2], a.0[3]);
        let (b0, b1, b2, b3) = (b.0[0], b.0[1], b.0[2], b.0[3]);
        let (s0, cc) = addcarry_u64(a0, b0, 0);
        let (s1, cc) = addcarry_u64(a1, b1, cc);
        let (s2, cc) = addcarry_u64(a2, b2, cc);
        let (s3, _) = addcarry_u64(a3, b3, cc);
        let (s0, s1, s2, s3) = Self::cond_sub_p(s0, s1, s2, s3);
        let (d0, cc) = subborrow_u64(a0, b0, 0);
        let (d1, cc) = subborrow_u64(a1, b1, cc);
        let (d2, cc) = subborrow_u64(a2, b2, cc);
        let (d3, cc) = subborrow_u64(a3, b3, cc);
        let m = (cc as u64).wrapping_neg();
        let (d0, cc) = addcarry_u64(d0, m & Self::MOD_X2[0], 0);
        let (d1, cc) = addcarry_u64(d1, m, cc);
        let (d2, cc) = addcarry_u64(d2, m, cc);
        let (d3, _) = addcarry_u64(d3, m & Self::MOD_X2[3], cc);
        (Self([s0, s1, s2, s3]), Self([d0, d1, d2, d3]))
    }

    /// `[a₀,a₁,a₂,a₃, b₀,b₁,b₂,b₃]` for the contiguous-pair asm kernels.
    #[inline(always)]
    fn pair(a: &Self, b: &[u64; 4]) -> [u64; 8] {
        [a.0[0], a.0[1], a.0[2], a.0[3], b[0], b[1], b[2], b[3]]
    }

    /// `(re, im)` of `a × b` over GF(p²), where `a = [a.re, a.im]` and
    /// `b = [b.re, b.im]` are passed as the `Fp2` memory layout (8 contiguous
    /// limbs each). Zero-copy: the asm reads `a`/`b` in place; only the
    /// `[2p − b.im, b.re]` scratch is materialised.
    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn fp2_mul_kernel(a: &[u64; 8], b: &[u64; 8]) -> ([u64; 4], [u64; 4]) {
        // im = a.re·b.im + a.im·b.re — kernel(A,B) = A.lo·B.hi + A.hi·B.lo,
        // so kernel(a, b) directly.
        let mut im = [0u64; 4];
        inline_asm::fp_sumprod(&mut im, a, b);
        // re = a.re·b.re − a.im·b.im — kernel(a, [2p−b.im, b.re]).
        let bim: &[u64; 4] = (&b[4..8]).try_into().unwrap();
        let nb = Self::raw_2p_minus(bim);
        let bb = [nb[0], nb[1], nb[2], nb[3], b[0], b[1], b[2], b[3]];
        let mut re = [0u64; 4];
        inline_asm::fp_sumprod(&mut re, a, &bb);
        (re, im)
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn sum_of_products(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        // inline kernel computes a[0..4]·b[4..8] + a[4..8]·b[0..4]
        let aa = Self::pair(a1, &a2.0);
        let bb = Self::pair(b2, &b1.0);
        let mut r = [0u64; 4];
        inline_asm::fp_sumprod(&mut r, &aa, &bb);
        Self(r)
    }

    /// Extern-C `.S` `fp_sop_asm` — test oracle for the inline `fp_sumprod`.
    #[cfg(gf5_248_asm)]
    pub fn sum_of_products_via_extern(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        let mut u = Self::ZERO;
        let a = [*a2, *a1];
        let b = [*b1, *b2];
        unsafe { fp_sop_asm(&mut u, a.as_ptr(), b.as_ptr()) }
        u
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    pub fn sum_of_products(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        let (a10, a11, a12, a13) = (a1.0[0], a1.0[1], a1.0[2], a1.0[3]);
        let (b10, b11, b12, b13) = (b1.0[0], b1.0[1], b1.0[2], b1.0[3]);
        let (a20, a21, a22, a23) = (a2.0[0], a2.0[1], a2.0[2], a2.0[3]);
        let (b20, b21, b22, b23) = (b2.0[0], b2.0[1], b2.0[2], b2.0[3]);

        // Product: a1 * b1 + a2 * b2 -> 504 bits
        let (e0, e1, cc) = Self::sum_of_products_u64(a10, b10, a20, b20, 0);
        let (e2, e3, cc) = Self::sum_of_products_u64(a11, b11, a21, b21, cc);
        let (e4, e5, cc) = Self::sum_of_products_u64(a12, b12, a22, b22, cc);
        let (e6, e7, _) = Self::sum_of_products_u64(a13, b13, a23, b23, cc);

        let (lo, hi, cc1) = Self::sum_of_products_u64(a10, b11, a20, b21, 0);
        let (e1, cc) = addcarry_u64(e1, lo, 0);
        let (e2, cc) = addcarry_u64(e2, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a10, b13, a20, b23, cc1);
        let (e3, cc) = addcarry_u64(e3, lo, cc);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a12, b13, a22, b23, cc1);
        let (e5, cc) = addcarry_u64(e5, lo, cc);
        let (e6, cc) = addcarry_u64(e6, hi, cc);
        let (e7, _) = addcarry_u64(e7, cc1 as u64, cc);

        let (lo, hi, cc1) = Self::sum_of_products_u64(a11, b10, a21, b20, 0);
        let (e1, cc) = addcarry_u64(e1, lo, 0);
        let (e2, cc) = addcarry_u64(e2, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a13, b10, a23, b20, cc1);
        let (e3, cc) = addcarry_u64(e3, lo, cc);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a13, b12, a23, b22, cc1);
        let (e5, cc) = addcarry_u64(e5, lo, cc);
        let (e6, cc) = addcarry_u64(e6, hi, cc);
        let (e7, _) = addcarry_u64(e7, cc1 as u64, cc);

        let (lo, hi, cc1) = Self::sum_of_products_u64(a10, b12, a20, b22, 0);
        let (e2, cc) = addcarry_u64(e2, lo, 0);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a11, b13, a21, b23, cc1);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (e6, cc) = addcarry_u64(e6, cc1 as u64, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo, hi, cc1) = Self::sum_of_products_u64(a12, b10, a22, b20, 0);
        let (e2, cc) = addcarry_u64(e2, lo, 0);
        let (e3, cc) = addcarry_u64(e3, hi, cc);
        let (lo, hi, cc1) = Self::sum_of_products_u64(a13, b11, a23, b21, cc1);
        let (e4, cc) = addcarry_u64(e4, lo, cc);
        let (e5, cc) = addcarry_u64(e5, hi, cc);
        let (e6, cc) = addcarry_u64(e6, cc1 as u64, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        let (lo1, hi1, cc1) = Self::sum_of_products_u64(a11, b12, a21, b22, 0);
        let (lo2, hi2, cc2) = Self::sum_of_products_u64(a12, b11, a22, b21, 0);
        let (lo, cc) = addcarry_u64(lo1, lo2, 0);
        let (hi, tt) = addcarry_u64(hi1, hi2, cc);
        let (e3, cc) = addcarry_u64(e3, lo, 0);
        let (e4, cc) = addcarry_u64(e4, hi, cc);
        let (e5, cc) = addcarry_u64(e5, tt as u64, cc);
        let (e6, cc) = addcarry_u64(e6, 0, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        // Add in the additional carries from sum_of_products_u64
        let (e5, cc) = addcarry_u64(e5, cc1 as u64, cc2);
        let (e6, cc) = addcarry_u64(e6, 0, cc);
        let (e7, _) = addcarry_u64(e7, 0, cc);

        // Montgomery reduction, see set_mul
        let f0 = e0;
        let f1 = e1;
        let f2 = e2;
        let f3 = e3.wrapping_add(e0.wrapping_mul(5) << 56);

        // Compute g = f*q.
        let (g3, hi) = umull(f0, 5u64 << 56);
        let (g4, hi) = umull_add(f1, 5u64 << 56, hi);
        let (g5, hi) = umull_add(f2, 5u64 << 56, hi);
        let (g6, g7) = umull_add(f3, 5u64 << 56, hi);
        let (g0, cc) = subborrow_u64(0, f0, 0);
        let (g1, cc) = subborrow_u64(0, f1, cc);
        let (g2, cc) = subborrow_u64(0, f2, cc);
        let (g3, cc) = subborrow_u64(g3, f3, cc);
        let (g4, cc) = subborrow_u64(g4, 0, cc);
        let (g5, cc) = subborrow_u64(g5, 0, cc);
        let (g6, cc) = subborrow_u64(g6, 0, cc);
        let (g7, _) = subborrow_u64(g7, 0, cc);

        // We add g = f*q to e0..e7. Since e0..e7 < 2^504, and f < 2^256,
        // we know that the result is not less than
        // 2^502 + 2^256*5*2^248 < 6*2^504; it is also a multiple of
        // 2^256. After dividing by 2^256, we get a value which is
        // less than 6*2^248, i.e. already in our proper range.
        let (_, cc) = addcarry_u64(e0, g0, 0);
        let (_, cc) = addcarry_u64(e1, g1, cc);
        let (_, cc) = addcarry_u64(e2, g2, cc);
        let (_, cc) = addcarry_u64(e3, g3, cc);
        let (d0, cc) = addcarry_u64(e4, g4, cc);
        let (d1, cc) = addcarry_u64(e5, g5, cc);
        let (d2, cc) = addcarry_u64(e6, g6, cc);
        let (d3, _) = addcarry_u64(e7, g7, cc);

        let mut r = Self::ZERO;
        r.0[0] = d0;
        r.0[1] = d1;
        r.0[2] = d2;
        r.0[3] = d3;
        r
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn difference_of_products(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        // a1·b1 + a2·(2p − b2)
        let aa = Self::pair(a1, &a2.0);
        let bb = Self::pair(&Self(Self::raw_2p_minus(&b2.0)), &b1.0);
        let mut r = [0u64; 4];
        inline_asm::fp_sumprod(&mut r, &aa, &bb);
        Self(r)
    }

    /// Extern-C `.S` `fp_dop_asm` — test oracle.
    #[cfg(gf5_248_asm)]
    pub fn difference_of_products_via_extern(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        let mut u = Self::ZERO;
        let a = [*a1, *a2];
        let b = [*b1, *b2];
        unsafe { fp_dop_asm(&mut u, a.as_ptr(), b.as_ptr()) }
        u
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    pub fn difference_of_products(a1: &Self, b1: &Self, a2: &Self, b2: &Self) -> Self {
        let b2_minus = b2.neg();
        Self::sum_of_products(a1, b1, a2, &b2_minus)
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn fp2_sq_re(x0: &Self, x1: &Self) -> Self {
        // .S writes scratch at [rdi+32..63]; needs 8-limb dst.
        let mut u = [Self::ZERO; 2];
        let a = [*x0, *x1];
        unsafe { fp2_sqr_re_asm(u.as_mut_ptr(), a.as_ptr()) }
        u[0]
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    // TODO: optimise?
    pub fn fp2_sq_re(x0: &Self, x1: &Self) -> Self {
        (x0 - x1) * (x0 + x1)
    }

    #[cfg(gf5_248_asm)]
    #[inline(always)]
    pub fn fp2_sq_im(x0: &Self, x1: &Self) -> Self {
        let mut u = Self::ZERO;
        let a = [*x0, *x1];
        unsafe { fp2_sqr_im_asm(&mut u, a.as_ptr()) }
        u
    }

    #[cfg(not(gf5_248_asm))]
    #[inline(always)]
    // TODO: optimise?
    pub fn fp2_sq_im(x0: &Self, x1: &Self) -> Self {
        (x0 * x1).mul2()
    }

    // Ensure that the internal encoding of this value is in the 0..q-1
    // range.
    #[inline(always)]
    fn set_normalized(&mut self) {
        // Subtract q.
        let (d0, cc) = subborrow_u64(self.0[0], Self::MODULUS[0], 0);
        let (d1, cc) = subborrow_u64(self.0[1], Self::MODULUS[1], cc);
        let (d2, cc) = subborrow_u64(self.0[2], Self::MODULUS[2], cc);
        let (d3, cc) = subborrow_u64(self.0[3], Self::MODULUS[3], cc);

        // Add back q if there was a borrow.
        let m = (cc as u64).wrapping_neg();
        let (d0, cc) = addcarry_u64(d0, m & Self::MODULUS[0], 0);
        let (d1, cc) = addcarry_u64(d1, m & Self::MODULUS[1], cc);
        let (d2, cc) = addcarry_u64(d2, m & Self::MODULUS[2], cc);
        let (d3, _) = addcarry_u64(d3, m & Self::MODULUS[3], cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    // Set this value to u*f+v*g (with 'u' being self). Parameters f and g
    // are provided as u64, but they are signed integers in the -2^62..+2^62
    // range.
    #[inline]
    fn set_lin(&mut self, u: &Self, v: &Self, f: u64, g: u64) {
        // Make sure f is nonnegative, by negating it if necessary, and
        // also negating u in that case to keep u*f unchanged.
        let sf = sgnw(f);
        let f = (f ^ sf).wrapping_sub(sf);
        let tu = Self::select(u, &-u, sf as u32);

        // Same treatment for g and v.
        let sg = sgnw(g);
        let g = (g ^ sg).wrapping_sub(sg);
        let tv = Self::select(v, &-v, sg as u32);

        // Compute the linear combination on plain integers. Since f and
        // g are at most 2^62 each, intermediate 128-bit products cannot
        // overflow.
        let (d0, t) = umull_x2(tu.0[0], f, tv.0[0], g);
        let (d1, t) = umull_x2_add(tu.0[1], f, tv.0[1], g, t);
        let (d2, t) = umull_x2_add(tu.0[2], f, tv.0[2], g, t);
        let (d3, t) = umull_x2_add(tu.0[3], f, tv.0[3], g, t);

        // Reduction: we split the value into a low part (248 bits) and
        // a high part (71 bits, since t can be up to 63 bits). If the
        // high part is h, then:
        //   h*2^248 = (h mod 5)*2^248 + floor(h / 5) mod q
        // We can compute the division by 5 by multiplying by the right
        // constant, and right-shifting.
        let h0 = (d3 >> 56) | (t << 8);
        let h1 = t >> 56;
        let d3 = d3 & 0x00FFFFFFFFFFFFFF;
        let (_, z) = umull(h0, 0xCCCCCCCCCCCCCCCD);
        let quo0 = z >> 2;
        let rem0 = h0 - (5 * quo0);
        let quo1 = (h1 * 0xCD) >> 10;
        let rem1 = h1 - (5 * quo1);

        // h = rem0 + 5*quo0 + (rem1 + 5*quo1)*2^64
        //   = rem0 + rem1 + 5*(quo0 + quo1*2^64 + rem1*((2^64 - 1)/5))
        // We add rem0 and rem1 modulo 5, with an extra carry that
        // goes into the folded part (multiple of 5).
        let (e, cc) = addcarry_u64(rem0 + 0xFFFFFFFFFFFFFFFA, rem1, 0);
        let (f0, cc) = addcarry_u64(quo0, rem1 * 0x3333333333333333, cc);
        let (f1, _) = addcarry_u64(quo1, 0, cc);
        let e = e - 0xFFFFFFFFFFFFFFFA;

        // Now we conly have to add e*2^248 + f0:f1 to the 248-bit low part.
        let (d0, cc) = addcarry_u64(d0, f0, 0);
        let (d1, cc) = addcarry_u64(d1, f1, cc);
        let (d2, cc) = addcarry_u64(d2, 0, cc);
        let (d3, _) = addcarry_u64(d3, e << 56, cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    #[inline(always)]
    fn lin(a: &Self, b: &Self, f: u64, g: u64) -> Self {
        let mut r = Self::ZERO;
        r.set_lin(a, b, f, g);
        r
    }

    // Set this value to abs((a*f+b*g)/2^31). Values a and b are interpreted
    // as unsigned 256-bit integers. Coefficients f and g are provided as u64,
    // but they really are signed integers in the -2^31..+2^31 range
    // (inclusive). The low 31 bits are dropped (i.e. the division is assumed
    // to be exact). The result is assumed to fit in 256 bits (including the
    // sign bit) (otherwise, truncation occurs).
    //
    // Returned value is -1 (u64) if (a*f+b*g) was negative, 0 otherwise.
    #[inline]
    fn set_lindiv31abs(&mut self, a: &Self, b: &Self, f: u64, g: u64) -> u64 {
        // Replace f and g with abs(f) and abs(g), but remember the
        // original signs.
        let sf = sgnw(f);
        let f = (f ^ sf).wrapping_sub(sf);
        let sg = sgnw(g);
        let g = (g ^ sg).wrapping_sub(sg);

        // Apply the signs of f and g to the source operands.
        let (a0, cc) = subborrow_u64(a.0[0] ^ sf, sf, 0);
        let (a1, cc) = subborrow_u64(a.0[1] ^ sf, sf, cc);
        let (a2, cc) = subborrow_u64(a.0[2] ^ sf, sf, cc);
        let (a3, cc) = subborrow_u64(a.0[3] ^ sf, sf, cc);
        let (a4, _) = subborrow_u64(0, 0, cc);
        let (b0, cc) = subborrow_u64(b.0[0] ^ sg, sg, 0);
        let (b1, cc) = subborrow_u64(b.0[1] ^ sg, sg, cc);
        let (b2, cc) = subborrow_u64(b.0[2] ^ sg, sg, cc);
        let (b3, cc) = subborrow_u64(b.0[3] ^ sg, sg, cc);
        let (b4, _) = subborrow_u64(0, 0, cc);

        // Compute a*f+b*g into d0:d1:d2:d3:t. Since f and g are at
        // most 2^31, we can add two 128-bit products with no overflow.
        let (d0, t) = umull_x2(a0, f, b0, g);
        let (d1, t) = umull_x2_add(a1, f, b1, g, t);
        let (d2, t) = umull_x2_add(a2, f, b2, g, t);
        let (d3, t) = umull_x2_add(a3, f, b3, g, t);
        // d4 <- a4*f + b4*g + t; a4 and b4 can be only 0 or -1
        let d4 = t.wrapping_sub(a4 & f).wrapping_sub(b4 & g);

        // Shift-right the value by 31 bits.
        let d0 = (d0 >> 31) | (d1 << 33);
        let d1 = (d1 >> 31) | (d2 << 33);
        let d2 = (d2 >> 31) | (d3 << 33);
        let d3 = (d3 >> 31) | (d4 << 33);

        // If the result is negative, then negate it.
        let t = sgnw(d4);
        let (d0, cc) = subborrow_u64(d0 ^ t, t, 0);
        let (d1, cc) = subborrow_u64(d1 ^ t, t, cc);
        let (d2, cc) = subborrow_u64(d2 ^ t, t, cc);
        let (d3, _) = subborrow_u64(d3 ^ t, t, cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
        t
    }

    #[inline(always)]
    fn lindiv31abs(a: &Self, b: &Self, f: u64, g: u64) -> (Self, u64) {
        let mut r = Self::ZERO;
        let ng = r.set_lindiv31abs(a, b, f, g);
        (r, ng)
    }

    fn set_div(&mut self, y: &Self) {
        // Extended binary GCD:
        //
        //   a <- y
        //   b <- q (modulus)
        //   u <- x (self)
        //   v <- 0
        //
        // Value a is normalized (in the 0..q-1 range). Values a and b are
        // then considered as (signed) integers. Values u and v are field
        // elements.
        //
        // Invariants:
        //    a*x = y*u mod q
        //    b*x = y*v mod q
        //    b is always odd
        //
        // At each step:
        //    if a is even, then:
        //        a <- a/2, u <- u/2 mod q
        //    else:
        //        if a < b:
        //            (a, u, b, v) <- (b, v, a, u)
        //        a <- (a-b)/2, u <- (u-v)/2 mod q
        //
        // What we implement below is the optimized version of this
        // algorithm, as described in https://eprint.iacr.org/2020/972

        let mut a = *y;
        a.set_normalized();
        let mut b = Self(Self::MODULUS);
        let mut u = *self;
        let mut v = Self::ZERO;

        // Generic loop does 15*31 = 465 inner iterations.
        for _ in 0..15 {
            // Get approximations of a and b over 64 bits:
            //  - If len(a) <= 64 and len(b) <= 64, then we just use
            //    their values (low limbs).
            //  - Otherwise, with n = max(len(a), len(b)), we use:
            //       (a mod 2^31) + 2^31*floor(a / 2^(n - 33))
            //       (b mod 2^31) + 2^31*floor(b / 2^(n - 33))

            let m3 = a.0[3] | b.0[3];
            let m2 = a.0[2] | b.0[2];
            let m1 = a.0[1] | b.0[1];
            let tnz3 = sgnw(m3 | m3.wrapping_neg());
            let tnz2 = sgnw(m2 | m2.wrapping_neg()) & !tnz3;
            let tnz1 = sgnw(m1 | m1.wrapping_neg()) & !tnz3 & !tnz2;
            let tnzm = (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1);
            let tnza = (a.0[3] & tnz3) | (a.0[2] & tnz2) | (a.0[1] & tnz1);
            let tnzb = (b.0[3] & tnz3) | (b.0[2] & tnz2) | (b.0[1] & tnz1);
            let snza = (a.0[2] & tnz3) | (a.0[1] & tnz2) | (a.0[0] & tnz1);
            let snzb = (b.0[2] & tnz3) | (b.0[1] & tnz2) | (b.0[0] & tnz1);

            // If both len(a) <= 64 and len(b) <= 64, then:
            //    tnzm = 0
            //    tnza = 0, snza = 0, tnzb = 0, snzb = 0
            // Otherwise:
            //    tnzm != 0
            //    tnza contains the top non-zero limb of a
            //    snza contains the limb right below tnza
            //    tnzb contains the top non-zero limb of a
            //    snzb contains the limb right below tnzb
            //
            // We count the number of leading zero bits in tnzm:
            //  - If s <= 31, then the top 31 bits can be extracted from
            //    tnza and tnzb alone.
            //  - If 32 <= s <= 63, then we need some bits from snza and
            //    snzb as well.
            let s = lzcnt(tnzm);
            let sm = (31_i32.wrapping_sub(s as i32) >> 31) as u64;
            let tnza = tnza ^ (sm & (tnza ^ ((tnza << 32) | (snza >> 32))));
            let tnzb = tnzb ^ (sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32))));
            let s = s - (32 & (sm as u32));
            let tnza = tnza << s;
            let tnzb = tnzb << s;

            // At this point:
            //  - If len(a) <= 64 and len(b) <= 64, then:
            //       tnza = 0
            //       tnzb = 0
            //       tnz1 = tnz2 = tnz3 = 0
            //       we want to use the entire low words of a and b
            //  - Otherwise, we want to use the top 33 bits of tnza and
            //    tnzb, and the low 31 bits of the low words of a and b.
            let tzx = !(tnz1 | tnz2 | tnz3);
            let tnza = tnza | (a.0[0] & tzx);
            let tnzb = tnzb | (b.0[0] & tzx);
            let mut xa = (a.0[0] & 0x7FFFFFFF) | (tnza & 0xFFFFFFFF80000000);
            let mut xb = (b.0[0] & 0x7FFFFFFF) | (tnzb & 0xFFFFFFFF80000000);

            // Compute the 31 inner iterations on xa and xb.
            let mut fg0 = 1u64;
            let mut fg1 = 1u64 << 32;
            for _ in 0..31 {
                let a_odd = (xa & 1).wrapping_neg();
                let (_, cc) = subborrow_u64(xa, xb, 0);
                let swap = a_odd & (cc as u64).wrapping_neg();
                let t1 = swap & (xa ^ xb);
                xa ^= t1;
                xb ^= t1;
                let t2 = swap & (fg0 ^ fg1);
                fg0 ^= t2;
                fg1 ^= t2;
                xa = xa.wrapping_sub(a_odd & xb);
                fg0 = fg0.wrapping_sub(a_odd & fg1);
                xa >>= 1;
                fg1 <<= 1;
            }
            fg0 = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
            fg1 = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
            let f0 = (fg0 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g0 = (fg0 >> 32).wrapping_sub(0x7FFFFFFF);
            let f1 = (fg1 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g1 = (fg1 >> 32).wrapping_sub(0x7FFFFFFF);

            // Propagate updates to a, b, u and v.
            let (na, nega) = Self::lindiv31abs(&a, &b, f0, g0);
            let (nb, negb) = Self::lindiv31abs(&a, &b, f1, g1);
            let f0 = (f0 ^ nega).wrapping_sub(nega);
            let g0 = (g0 ^ nega).wrapping_sub(nega);
            let f1 = (f1 ^ negb).wrapping_sub(negb);
            let g1 = (g1 ^ negb).wrapping_sub(negb);
            let nu = Self::lin(&u, &v, f0, g0);
            let nv = Self::lin(&u, &v, f1, g1);
            a = na;
            b = nb;
            u = nu;
            v = nv;
        }

        // If y is invertible, then the final GCD is 1, and
        // len(a) + len(b) <= 37, so we can end the computation with
        // the low words directly. We only need 35 iterations to reach
        // the point where b = 1.
        //
        // If y is zero, then v is unchanged (hence zero) and none of
        // the subsequent iterations will change it either, so we get
        // 0 on output, which is what we want.
        let mut xa = a.0[0];
        let mut xb = b.0[0];
        let mut f0 = 1u64;
        let mut g0 = 0u64;
        let mut f1 = 0u64;
        let mut g1 = 1u64;
        for _ in 0..35 {
            let a_odd = (xa & 1).wrapping_neg();
            let (_, cc) = subborrow_u64(xa, xb, 0);
            let swap = a_odd & (cc as u64).wrapping_neg();
            let t1 = swap & (xa ^ xb);
            xa ^= t1;
            xb ^= t1;
            let t2 = swap & (f0 ^ f1);
            f0 ^= t2;
            f1 ^= t2;
            let t3 = swap & (g0 ^ g1);
            g0 ^= t3;
            g1 ^= t3;
            xa = xa.wrapping_sub(a_odd & xb);
            f0 = f0.wrapping_sub(a_odd & f1);
            g0 = g0.wrapping_sub(a_odd & g1);
            xa >>= 1;
            f1 <<= 1;
            g1 <<= 1;
        }

        self.set_lin(&u, &v, f1, g1);

        // At the point:
        //  - Numerator and denominator were both in Montgomery representation,
        //    but the two factors R canceled each other.
        //  - We have injected 31*15+35 = 500 extra factors of 2, hence we
        //    must divide the result by 2^500.
        //  - However, we also want to obtain the result in Montgomery
        //    representation, i.e. multiply by 2^256. We thus want to
        //    divide the current result by 2^(500 - 256) = 2^244.
        //  - We do this division by using a Montgomery multiplication with
        //    the Montgomery representation of 1/2^244, i.e. the integer
        //    2^256/2^244 = 4096.

        // At this point, we have injected extra factors of 2, one for
        // each of the 31*15+35 = 500 iterations, so we must divide by
        // 2^500 (mod q). This is done with a multiplication by the
        // appropriate constant.
        self.set_mul(&Self::INVT244);
    }

    pub fn set_invert(&mut self) {
        let r = *self;
        *self = Self::ONE;
        self.set_div(&r);
    }

    pub fn invert(self) -> Self {
        let mut r = Self::ONE;
        r.set_div(&self);
        r
    }

    // Perform a batch inversion of some elements. All elements of
    // the slice are replaced with their respective inverse (elements
    // of value zero are "inverted" into themselves).
    pub fn batch_invert(xx: &mut [Self]) {
        // We use Montgomery's trick:
        //   1/u = v*(1/(u*v))
        //   1/v = u*(1/(u*v))
        // Applied recursively on n elements, this computes an inversion
        // with a single inversion in the field, and 3*(n-1) multiplications.
        // We use batches of 200 elements; larger batches only yield
        // moderate improvements, while sticking to a fixed moderate batch
        // size allows stack-based allocation.
        let n = xx.len();
        let mut i = 0;
        while i < n {
            let blen = if (n - i) > 200 { 200 } else { n - i };
            let mut tt = [Self::ZERO; 200];
            tt[0] = xx[i];
            let zz0 = tt[0].is_zero();
            tt[0].set_cond(&Self::ONE, zz0);
            for j in 1..blen {
                tt[j] = xx[i + j];
                tt[j].set_cond(&Self::ONE, tt[j].is_zero());
                tt[j] *= tt[j - 1];
            }
            let mut k = Self::ONE / tt[blen - 1];
            for j in (1..blen).rev() {
                let mut x = xx[i + j];
                let zz = x.is_zero();
                x.set_cond(&Self::ONE, zz);
                xx[i + j].set_cond(&(k * tt[j - 1]), !zz);
                k *= x;
            }
            xx[i].set_cond(&k, !zz0);
            i += blen;
        }
    }

    // Compute the Legendre symbol on this value. Return value is:
    //   0   if this value is zero
    //  +1   if this value is a non-zero quadratic residue
    //  -1   if this value is not a quadratic residue
    pub fn legendre(self) -> i32 {
        // The algorithm is very similar to the optimized binary GCD that
        // is implemented in set_div(), with the following differences:
        //  - We do not keep track of the 'u' and 'v' values.
        //  - In each inner iteration, the running symbol value is
        //    adjusted, taking into account the low 2 or 3 bits of the
        //    involved values.
        //  - Since we need a couple of bits of look-ahead, we can only
        //    run 29 iterations in the inner loop, and we need an extra
        //    recomputation step for the next 2.
        // Otherwise, the 'a' and 'b' values are modified exactly as in
        // the binary GCD, so that we get the same guaranteed convergence
        // in a total of 510 iterations.

        let mut a = self;
        a.set_normalized();
        let mut b = Self(Self::MODULUS);
        let mut ls = 0u64; // running symbol information in the low bit

        // Outer loop
        for _ in 0..15 {
            // Get approximations of a and b over 64 bits.
            let m3 = a.0[3] | b.0[3];
            let m2 = a.0[2] | b.0[2];
            let m1 = a.0[1] | b.0[1];
            let tnz3 = sgnw(m3 | m3.wrapping_neg());
            let tnz2 = sgnw(m2 | m2.wrapping_neg()) & !tnz3;
            let tnz1 = sgnw(m1 | m1.wrapping_neg()) & !tnz3 & !tnz2;
            let tnzm = (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1);
            let tnza = (a.0[3] & tnz3) | (a.0[2] & tnz2) | (a.0[1] & tnz1);
            let tnzb = (b.0[3] & tnz3) | (b.0[2] & tnz2) | (b.0[1] & tnz1);
            let snza = (a.0[2] & tnz3) | (a.0[1] & tnz2) | (a.0[0] & tnz1);
            let snzb = (b.0[2] & tnz3) | (b.0[1] & tnz2) | (b.0[0] & tnz1);

            let s = lzcnt(tnzm);
            let sm = (31_i32.wrapping_sub(s as i32) >> 31) as u64;
            let tnza = tnza ^ (sm & (tnza ^ ((tnza << 32) | (snza >> 32))));
            let tnzb = tnzb ^ (sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32))));
            let s = s - (32 & (sm as u32));
            let tnza = tnza << s;
            let tnzb = tnzb << s;

            let tzx = !(tnz1 | tnz2 | tnz3);
            let tnza = tnza | (a.0[0] & tzx);
            let tnzb = tnzb | (b.0[0] & tzx);
            let mut xa = (a.0[0] & 0x7FFFFFFF) | (tnza & 0xFFFFFFFF80000000);
            let mut xb = (b.0[0] & 0x7FFFFFFF) | (tnzb & 0xFFFFFFFF80000000);

            // First 29 inner iterations.
            let mut fg0 = 1u64;
            let mut fg1 = 1u64 << 32;
            for _ in 0..29 {
                let a_odd = (xa & 1).wrapping_neg();
                let (_, cc) = subborrow_u64(xa, xb, 0);
                let swap = a_odd & (cc as u64).wrapping_neg();
                ls ^= swap & ((xa & xb) >> 1);
                let t1 = swap & (xa ^ xb);
                xa ^= t1;
                xb ^= t1;
                let t2 = swap & (fg0 ^ fg1);
                fg0 ^= t2;
                fg1 ^= t2;
                xa = xa.wrapping_sub(a_odd & xb);
                fg0 = fg0.wrapping_sub(a_odd & fg1);
                xa >>= 1;
                fg1 <<= 1;
                ls ^= xb.wrapping_add(2) >> 2;
            }

            // Compute the updated a and b (low words only) to get enough
            // bits for the next two iterations.
            let fg0z = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
            let fg1z = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
            let f0 = (fg0z & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g0 = (fg0z >> 32).wrapping_sub(0x7FFFFFFF);
            let f1 = (fg1z & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g1 = (fg1z >> 32).wrapping_sub(0x7FFFFFFF);
            let mut a0 = a.0[0]
                .wrapping_mul(f0)
                .wrapping_add(b.0[0].wrapping_mul(g0))
                >> 29;
            let mut b0 = a.0[0]
                .wrapping_mul(f1)
                .wrapping_add(b.0[0].wrapping_mul(g1))
                >> 29;
            for _ in 0..2 {
                let a_odd = (xa & 1).wrapping_neg();
                let (_, cc) = subborrow_u64(xa, xb, 0);
                let swap = a_odd & (cc as u64).wrapping_neg();
                ls ^= swap & ((a0 & b0) >> 1);
                let t1 = swap & (xa ^ xb);
                xa ^= t1;
                xb ^= t1;
                let t2 = swap & (fg0 ^ fg1);
                fg0 ^= t2;
                fg1 ^= t2;
                let t3 = swap & (a0 ^ b0);
                a0 ^= t3;
                b0 ^= t3;
                xa = xa.wrapping_sub(a_odd & xb);
                fg0 = fg0.wrapping_sub(a_odd & fg1);
                a0 = a0.wrapping_sub(a_odd & b0);
                xa >>= 1;
                fg1 <<= 1;
                a0 >>= 1;
                ls ^= b0.wrapping_add(2) >> 2;
            }

            // Propagate updates to a and b.
            fg0 = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
            fg1 = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
            let f0 = (fg0 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g0 = (fg0 >> 32).wrapping_sub(0x7FFFFFFF);
            let f1 = (fg1 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
            let g1 = (fg1 >> 32).wrapping_sub(0x7FFFFFFF);

            let (na, nega) = Self::lindiv31abs(&a, &b, f0, g0);
            let (nb, _) = Self::lindiv31abs(&a, &b, f1, g1);
            ls ^= nega & (nb.0[0] >> 1);
            a = na;
            b = nb;
        }

        // Final iterations: values are at most 37 bits now. We do not
        // need to keep track of update coefficients. Just like the GCD,
        // we need only 35 iterations, because after 35 iterations,
        // value a is 0 or 1, and b is 1, and no further modification to
        // the Legendre symbol may happen.
        let mut xa = a.0[0];
        let mut xb = b.0[0];
        for _ in 0..35 {
            let a_odd = (xa & 1).wrapping_neg();
            let (_, cc) = subborrow_u64(xa, xb, 0);
            let swap = a_odd & (cc as u64).wrapping_neg();
            ls ^= swap & ((xa & xb) >> 1);
            let t1 = swap & (xa ^ xb);
            xa ^= t1;
            xb ^= t1;
            xa = xa.wrapping_sub(a_odd & xb);
            xa >>= 1;
            ls ^= xb.wrapping_add(2) >> 2;
        }

        // At this point, if the source value was not zero, then the low
        // bit of ls contains the QR status (0 = square, 1 = non-square),
        // which we need to convert to the expected value (+1 or -1).
        // If y == 0, then we return 0, per the API.
        let r = 1u32.wrapping_sub(((ls as u32) & 1) << 1);
        (r & !self.is_zero()) as i32
    }

    /// Return `0xFFFFFFFF` when this value is a square in GF(p^2) and
    /// `0x00000000` otherwise.
    #[inline]
    pub fn is_square(self) -> u32 {
        !((self.legendre() >> 1) as u32)
    }

    // Set this value to its square root. Returned value is 0xFFFFFFFF
    // if the operation succeeded (value was indeed a quadratic
    // residue), 0 otherwise (value was not a quadratic residue). In the
    // latter case, this value is set to the square root of -self. In
    // all cases, the returned root is the one whose least significant
    // bit is 0 (when normalized in 0..q-1).
    fn set_sqrt_ext(&mut self) -> u32 {
        // Candidate root is self^((q+1)/4).
        // (q+1)/4 = 5*2^246
        let z = *self;
        let z5 = z * z.n_square(2);
        let mut y = z5.n_square(246);

        // Normalize y and negate it if necessary to set the low bit to 0.
        // We must take care to make the bit check on the value in normal
        // representation, not Montgomery representation.
        let mut yn = y;
        yn.set_montgomery_reduce();
        y.set_cond(&-y, ((yn.0[0] as u32) & 1).wrapping_neg());

        // Check that the candidate is indeed a square root.
        let r = y.square().equals(self);
        *self = y;
        r
    }

    // Set this value to its square root. Returned value is 0xFFFFFFFF
    // if the operation succeeded (value was indeed a quadratic
    // residue), 0 otherwise (value was not a quadratic residue). This
    // differs from set_sqrt_ext() in that this function sets the value
    // to zero if there is no square root.
    fn set_sqrt(&mut self) -> u32 {
        let r = self.set_sqrt_ext();
        self.set_cond(&Self::ZERO, !r);
        r
    }

    pub fn set_fourth_root(&mut self) -> u32 {
        // Candidate root is self^((q+1)/8).
        // (q+1)/4 = 5*2^246
        let z = *self;
        let z5 = z * z.n_square(2);
        let mut y = z5.n_square(245);

        // Normalize y and negate it if necessary to set the low bit to 0.
        // We must take care to make the bit check on the value in normal
        // representation, not Montgomery representation.
        let mut yn = y;
        yn.set_montgomery_reduce();
        y.set_cond(&-y, ((yn.0[0] as u32) & 1).wrapping_neg());

        // Check that the candidate is indeed a square root.
        let r = y.n_square(2).equals(self);
        *self = y;
        r
    }

    pub fn fourth_root(&self) -> (Self, u32) {
        let mut r = *self;
        let ok = r.set_fourth_root();
        r.set_cond(&Self::ZERO, !ok);
        (r, ok)
    }

    // Compute the square root of this value. Returned value are (y, r):
    //  - If this value is indeed a quadratic residue, then y is the
    //    square root whose least significant bit (when normalized in 0..q-1)
    //    is 0, and r is equal to 0xFFFFFFFF.
    //  - If this value is not a quadratic residue, then y is zero, and
    //    r is equal to 0.
    #[inline(always)]
    pub fn sqrt(self) -> (Self, u32) {
        let mut x = self;
        let r = x.set_sqrt();
        (x, r)
    }

    // Compute the square root of this value. Returned value are (y, r):
    //  - If this value is indeed a quadratic residue, then y is a
    //    square root of this value, and r is 0xFFFFFFFF.
    //  - If this value is not a quadratic residue, then y is set to
    //    a square root of -x, and r is 0x00000000.
    // In all cases, the returned root is normalized: the lest significant
    // bit of its integer representation (in the 0..q-1 range) is 0.
    #[inline(always)]
    pub fn sqrt_ext(self) -> (Self, u32) {
        let mut x = self;
        let r = x.set_sqrt_ext();
        (x, r)
    }

    /* TODO? Not very useful on this field.
    // Compute two signed integers (c0, c1) such that this self = c0/c1 in
    // the ring. The two returned values may range up to about 1.2*2^126
    // each (in absolute value).
    //
    // THIS FUNCTION IS NOT CONSTANT-TIME. It shall be used only for a
    // public source element.
    pub fn split_vartime(self) -> (i128, i128) {
        let mut k = self;
        k.set_normalized();
        lagrange253_vartime(&k.0, &Self::MODULUS)
    }
    */

    // Equality check between two field elements (constant-time);
    // returned value is 0xFFFFFFFF on equality, 0 otherwise.
    #[inline(always)]
    pub fn equals(self, rhs: &Self) -> u32 {
        (self - rhs).is_zero()
    }

    // Compare this value with zero (constant-time); returned value
    // is 0xFFFFFFFF if this element is zero, 0 otherwise.
    #[inline(always)]
    pub fn is_zero(self) -> u32 {
        // There are two possible representations for 0: 0 and q.
        let a0 = self.0[0];
        let a1 = self.0[1];
        let a2 = self.0[2];
        let a3 = self.0[3];
        let t0 = a0 | a1 | a2 | a3;
        let t1 = (a0 ^ Self::MODULUS[0])
            | (a1 ^ Self::MODULUS[1])
            | (a2 ^ Self::MODULUS[2])
            | (a3 ^ Self::MODULUS[3]);

        // Top bit of r is 0 if and only if one of t0 or t1 is zero.
        let r = (t0 | t0.wrapping_neg()) & (t1 | t1.wrapping_neg());
        ((r >> 63) as u32).wrapping_sub(1)
    }

    // This internal function decodes 32 bytes (exactly) into a 256-bit
    // integer.
    #[inline(always)]
    fn set_decode32_raw(&mut self, buf: &[u8]) {
        debug_assert!(buf.len() == 32);
        self.0[0] = u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[0..8]).unwrap());
        self.0[1] = u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[8..16]).unwrap());
        self.0[2] = u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[16..24]).unwrap());
        self.0[3] = u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[24..32]).unwrap());
    }

    // This internal function decodes 32 bytes (exactly) into a 256-bit
    // integer, and partially reduces it modulo q to ensure a 251-bit value.
    // WARNING: Montgomery representation is not applied.
    #[inline]
    fn set_decode32_reduce(&mut self, buf: &[u8]) {
        // Decode the 32 bytes as a raw integer.
        self.set_decode32_raw(buf);

        // Partial reduction.
        self.set_partial_reduce();
    }

    /// Decode an element from bytes, no check is made that the input
    /// value is reduced except that the buffer is of the excpected
    /// length of `Self::ENCODED_LENGTH`.
    pub const fn const_decode_no_check(buf: &[u8]) -> Self {
        let mut r = Self::ZERO;
        if buf.len() != Self::ENCODED_LENGTH {
            return r;
        }

        // TODO: clean up!
        // Fill the first N-1 elements
        let mut i = 0;
        while i < Self::N - 1 {
            r.0[i] = u64::from_le_bytes([
                buf[i * 8],
                buf[i * 8 + 1],
                buf[i * 8 + 2],
                buf[i * 8 + 3],
                buf[i * 8 + 4],
                buf[i * 8 + 5],
                buf[i * 8 + 6],
                buf[i * 8 + 7],
            ]);
            i += 1;
        }

        // Fill the last element
        let mut w = 0u64;
        let mut j = 0;
        while j < Self::ENCODED_LENGTH - (Self::N - 1) * 8 {
            w |= (buf[(Self::N - 1) * 8 + j] as u64) << (8 * j);
            j += 1;
        }
        r.0[Self::N - 1] = w;

        Self::const_mmul(r, Self::R2)
    }

    // Apply partial reduction modulo q so that the internal state fits
    // on 251 bits.
    // WARNING: Montgomery representation is not applied.
    #[inline(always)]
    fn set_partial_reduce(&mut self) {
        // We need to reduce the value into an acceptable range (251 bits).
        // Extract low 248-bit part, and the high part (8 bits).
        let h = self.0[3] >> 56;
        let d3 = self.0[3] & 0x00FFFFFFFFFFFFFF;

        // Fold value h by dividing it by 5. Since 5*2^248 = 1 mod q,
        // We add floor(h/5) + (h mod 5)*2^248 to the value.
        let quo = (h * 0xCD) >> 10;
        let rem = h - (5 * quo);
        let (d0, cc) = addcarry_u64(self.0[0], quo, 0);
        let (d1, cc) = addcarry_u64(self.0[1], 0, cc);
        let (d2, cc) = addcarry_u64(self.0[2], 0, cc);
        let (d3, _) = addcarry_u64(d3, rem << 56, cc);

        self.0[0] = d0;
        self.0[1] = d1;
        self.0[2] = d2;
        self.0[3] = d3;
    }

    // Encode this value over exactly 32 bytes. Encoding is always canonical
    // (little-endian encoding of the value in the 0..q-1 range, high four
    // bits of the last byte are always 0).
    #[inline(always)]
    pub fn encode32(self) -> [u8; 32] {
        let mut r = self;
        r.set_montgomery_reduce();
        let mut d = [0u8; 32];
        d[0..8].copy_from_slice(&r.0[0].to_le_bytes());
        d[8..16].copy_from_slice(&r.0[1].to_le_bytes());
        d[16..24].copy_from_slice(&r.0[2].to_le_bytes());
        d[24..32].copy_from_slice(&r.0[3].to_le_bytes());
        d
    }

    // Encode this value over exactly 32 bytes. Encoding is always canonical
    // (little-endian encoding of the value in the 0..q-1 range, high four
    // bits of the last byte are always 0).
    #[inline(always)]
    pub fn encode(self) -> [u8; 32] {
        self.encode32()
    }

    // Decode the field element from the provided bytes. If the source
    // slice does not have length exactly 32 bytes, or if the encoding
    // is non-canonical (i.e. does not represent an integer in the 0
    // to q-1 range), then this element is set to zero, and 0 is returned.
    // Otherwise, this element is set to the decoded value, and 0xFFFFFFFF
    // is returned.
    #[inline]
    pub fn set_decode_ct(&mut self, buf: &[u8]) -> u32 {
        if buf.len() != 32 {
            *self = Self::ZERO;
            return 0;
        }

        self.set_decode32_raw(buf);

        // Try to subtract q from the value; if that does not yield a
        // borrow, then the encoding was not canonical.
        let (_, cc) = subborrow_u64(self.0[0], Self::MODULUS[0], 0);
        let (_, cc) = subborrow_u64(self.0[1], Self::MODULUS[1], cc);
        let (_, cc) = subborrow_u64(self.0[2], Self::MODULUS[2], cc);
        let (_, cc) = subborrow_u64(self.0[3], Self::MODULUS[3], cc);

        // Clear the value if not canonical.
        let cc = (cc as u64).wrapping_neg();
        self.0[0] &= cc;
        self.0[1] &= cc;
        self.0[2] &= cc;
        self.0[3] &= cc;

        // Convert to Montgomery representation.
        self.set_mul(&Self::R2);

        cc as u32
    }

    // Decode a field element from 32 bytes. On success, this returns
    // (r, cc), where cc has value 0xFFFFFFFF. If the source encoding is not
    // canonical (i.e. the unsigned little-endian interpretation of the
    // 32 bytes yields an integer with is not lower than q), then this
    // returns (0, 0).
    #[inline(always)]
    pub fn decode_ct(buf: &[u8]) -> (Self, u32) {
        let mut r = Self::ZERO;
        let cc = r.set_decode_ct(buf);
        (r, cc)
    }

    // Decode a field element from 32 bytes. On success, this returns
    // (r, cc), where cc has value 0xFFFFFFFF. If the source encoding is not
    // canonical (i.e. the unsigned little-endian interpretation of the
    // 32 bytes yields an integer with is not lower than q), then this
    // returns (0, 0).
    #[inline]
    pub fn decode(buf: &[u8]) -> (Self, u32) {
        Self::decode_ct(buf)
    }

    // Decode a field element from some bytes. The bytes are interpreted
    // in unsigned little-endian convention, and the resulting integer
    // is reduced modulo q. This process never fails.
    pub fn set_decode_reduce(&mut self, buf: &[u8]) {
        *self = Self::ZERO;
        let mut n = buf.len();
        if n == 0 {
            return;
        }
        if (n & 31) != 0 {
            // If input size is not a multiple of 32, then we decode a
            // partial block and the value is already less than 2^248.
            let k = n & !(31_usize);
            let mut tmp = [0u8; 32];
            tmp[..(n - k)].copy_from_slice(&buf[k..]);
            n = k;
            self.set_decode32_raw(&tmp);
        } else {
            // If input size is a multiple of 32, then we decode a full
            // 32-byte block, and we must reduce it to get a 251-bit value.
            n -= 32;
            self.set_decode32_reduce(&buf[n..]);
        }

        // Process all remaining blocks, in descending address order.
        while n > 0 {
            // For each block, we multiply the current value by 2^256
            // (i.e. we Montgomery-multiply by 2^512) then we add the
            // next block.
            let k = n - 32;
            let mut v = Self::ZERO;
            v.set_decode32_reduce(&buf[k..k + 32]);
            self.set_mul(&Self::R2);
            self.set_add(&v);
            n = k;
        }

        // We decoded the value as a plain integer, we must convert it
        // to Montgomery representation.
        self.set_mul(&Self::R2);
    }

    // Decode a field element from some bytes. The bytes are interpreted
    // in unsigned little-endian convention, and the resulting integer
    // is reduced modulo q. This process never fails.
    #[inline(always)]
    pub fn decode_reduce(buf: &[u8]) -> Self {
        let mut r = Self::ZERO;
        r.set_decode_reduce(buf);
        r
    }

    /// Get the "hash" of the value (low 64 bits of the Montgomery
    /// representation).
    pub fn hashcode(self) -> u64 {
        self.0[0]
    }

    // Partial reduction modulo q, as a function usable in constant
    // contexts. It is safe (constant-time) and thus also usable at
    // runtime, but less efficient than set_mul() since it does not
    // leverage intrinsics.
    const fn const_mred(a0: u64, a1: u64, a2: u64, a3: u64) -> [u64; 4] {
        // Custom add-with-carry.
        const fn adc(x: u64, y: u64, cc: u64) -> (u64, u64) {
            let z = (x as u128).wrapping_add(y as u128).wrapping_add(cc as u128);
            (z as u64, (z >> 64) as u64)
        }

        // Extract low 248-bit part, and the high part (8 bits).
        let h = a3 >> 56;
        let a3 = a3 & 0x00FFFFFFFFFFFFFF;

        // Fold value h by dividing it by 5. Since 5*2^248 = 1 mod q,
        // We add floor(h/5) + (h mod 5)*2^248 to the value.
        let quo = (h * 0xCD) >> 10;
        let rem = h - (5 * quo);
        let (d0, cc) = adc(a0, quo, 0);
        let (d1, cc) = adc(a1, 0, cc);
        let (d2, cc) = adc(a2, 0, cc);
        let (d3, _) = adc(a3, rem << 56, cc);

        [d0, d1, d2, d3]
    }

    // Montgomery multiplication, as a function usable in constant
    // contexts. It is safe (constant-time) and thus also usable at
    // runtime, but less efficient than set_mul() since it does not
    // leverage intrinsics.
    //
    // Special: this function also tolerates that input a[] ranges
    // up to 2^256-1 (full 256-bit word), as long as b[] is less than q.
    const fn const_mmul(a: Self, b: Self) -> Self {
        // Custom add-with-carry.
        const fn adc(x: u64, y: u64, cc: u64) -> (u64, u64) {
            let z = (x as u128).wrapping_add(y as u128).wrapping_add(cc as u128);
            (z as u64, (z >> 64) as u64)
        }

        // Compute x*y + a + b, returned over two words (lo, hi).
        const fn umaal(x: u64, y: u64, a: u64, b: u64) -> (u64, u64) {
            let z = (x as u128) * (y as u128) + (a as u128) + (b as u128);
            (z as u64, (z >> 64) as u64)
        }

        // Given d0..d3 (with d <= 2*q-1), operand b[] (b <= q-1) and
        // multiplier aj, return ((d + aj*b) / 2^64) mod m, partially
        // reduced (output is at most 2*q-1).
        const fn mmul1(
            aj: u64,
            b: [u64; 4],
            d0: u64,
            d1: u64,
            d2: u64,
            d3: u64,
        ) -> (u64, u64, u64, u64) {
            // d <- d + a*bj (may range up to (2^64+1)*q, needs 5 words)
            let (d0, hi) = umaal(aj, b[0], d0, 0);
            let (d1, hi) = umaal(aj, b[1], d1, hi);
            let (d2, hi) = umaal(aj, b[2], d2, hi);
            let (d3, d4) = umaal(aj, b[3], d3, hi);
            let f = d0;
            let (_, hi) = umaal(f, 0xFFFFFFFFFFFFFFFF, d0, 0);
            let (d0, hi) = umaal(f, 0xFFFFFFFFFFFFFFFF, d1, hi);
            let (d1, hi) = umaal(f, 0xFFFFFFFFFFFFFFFF, d2, hi);
            let (d2, hi) = umaal(f, 0x04FFFFFFFFFFFFFF, d3, hi);
            let (d3, _) = adc(d4, hi, 0);
            (d0, d1, d2, d3)
        }

        let (d0, d1, d2, d3) = (0u64, 0u64, 0u64, 0u64);
        let (d0, d1, d2, d3) = mmul1(a.0[0], b.0, d0, d1, d2, d3);
        let (d0, d1, d2, d3) = mmul1(a.0[1], b.0, d0, d1, d2, d3);
        let (d0, d1, d2, d3) = mmul1(a.0[2], b.0, d0, d1, d2, d3);
        let (d0, d1, d2, d3) = mmul1(a.0[3], b.0, d0, d1, d2, d3);

        Self(Self::const_mred(d0, d1, d2, d3))
    }
}

// ========================================================================
// Implementations of all the traits needed to use the simple operators
// (+, *, /...) on field element instances, with or without references.

impl Add<GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn add(self, other: GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_add(&other);
        r
    }
}

impl Add<&GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn add(self, other: &GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_add(other);
        r
    }
}

impl Add<GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn add(self, other: GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_add(&other);
        r
    }
}

impl Add<&GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn add(self, other: &GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_add(other);
        r
    }
}

impl AddAssign<GF5_248> for GF5_248 {
    #[inline(always)]
    fn add_assign(&mut self, other: GF5_248) {
        self.set_add(&other);
    }
}

impl AddAssign<&GF5_248> for GF5_248 {
    #[inline(always)]
    fn add_assign(&mut self, other: &GF5_248) {
        self.set_add(other);
    }
}

impl Div<GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn div(self, other: GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_div(&other);
        r
    }
}

impl Div<&GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn div(self, other: &GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_div(other);
        r
    }
}

impl Div<GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn div(self, other: GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_div(&other);
        r
    }
}

impl Div<&GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn div(self, other: &GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_div(other);
        r
    }
}

impl DivAssign<GF5_248> for GF5_248 {
    #[inline(always)]
    fn div_assign(&mut self, other: GF5_248) {
        self.set_div(&other);
    }
}

impl DivAssign<&GF5_248> for GF5_248 {
    #[inline(always)]
    fn div_assign(&mut self, other: &GF5_248) {
        self.set_div(other);
    }
}

impl Mul<GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn mul(self, other: GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_mul(&other);
        r
    }
}

impl Mul<&GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn mul(self, other: &GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_mul(other);
        r
    }
}

impl Mul<GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn mul(self, other: GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_mul(&other);
        r
    }
}

impl Mul<&GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn mul(self, other: &GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_mul(other);
        r
    }
}

impl MulAssign<GF5_248> for GF5_248 {
    #[inline(always)]
    fn mul_assign(&mut self, other: GF5_248) {
        self.set_mul(&other);
    }
}

impl MulAssign<&GF5_248> for GF5_248 {
    #[inline(always)]
    fn mul_assign(&mut self, other: &GF5_248) {
        self.set_mul(other);
    }
}

impl Neg for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn neg(self) -> GF5_248 {
        let mut r = self;
        r.set_neg();
        r
    }
}

impl Neg for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn neg(self) -> GF5_248 {
        let mut r = *self;
        r.set_neg();
        r
    }
}

impl Sub<GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn sub(self, other: GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_sub(&other);
        r
    }
}

impl Sub<&GF5_248> for GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn sub(self, other: &GF5_248) -> GF5_248 {
        let mut r = self;
        r.set_sub(other);
        r
    }
}

impl Sub<GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn sub(self, other: GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_sub(&other);
        r
    }
}

impl Sub<&GF5_248> for &GF5_248 {
    type Output = GF5_248;

    #[inline(always)]
    fn sub(self, other: &GF5_248) -> GF5_248 {
        let mut r = *self;
        r.set_sub(other);
        r
    }
}

impl SubAssign<GF5_248> for GF5_248 {
    #[inline(always)]
    fn sub_assign(&mut self, other: GF5_248) {
        self.set_sub(&other);
    }
}

impl SubAssign<&GF5_248> for GF5_248 {
    #[inline(always)]
    fn sub_assign(&mut self, other: &GF5_248) {
        self.set_sub(other);
    }
}

impl ::std::fmt::Display for GF5_248 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        // Encode the value to get a canonical value
        let v_bytes = self.encode();

        // Convert the bytes into a base-16 string
        let mut v_hex_string: String = v_bytes
            .iter()
            .rev()
            .map(|byte| format!("{:02x}", byte)) // Format each byte as a two-digit hex
            .collect();

        // Remove leading zeros
        v_hex_string = v_hex_string.trim_start_matches("0").to_string();

        // If the value was zero, we need to add a zero back
        if v_hex_string.is_empty() {
            v_hex_string = String::from("0");
        }

        write!(f, "0x{}", v_hex_string)
    }
}

impl From<u64> for GF5_248 {
    fn from(x: u64) -> GF5_248 {
        let mut r = Self::ZERO;
        r.0[0] = x;
        r.set_mul(&Self::R2);
        r
    }
}

impl From<i64> for GF5_248 {
    fn from(x: i64) -> GF5_248 {
        let sx = (x >> 63) as u64;
        let ax = ((x as u64) ^ sx).wrapping_sub(sx);
        let mut r = Self::from(ax);
        r.set_cond_neg(sx as u32);
        r
    }
}

impl From<u32> for GF5_248 {
    fn from(x: u32) -> GF5_248 {
        Self::from(x as u64)
    }
}

impl From<i32> for GF5_248 {
    fn from(x: i32) -> GF5_248 {
        let mut r = Self::ONE;
        r.set_mul_small(x);
        r
    }
}
