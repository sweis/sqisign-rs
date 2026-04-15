//! GF(p) backend via the `fp2` crate's `define_fp_core!` macro
//! (saturated 4-limb Montgomery, p = 5·2²⁴⁸ − 1).
//!
//! This file replaces `fp.rs` when feature `gf-fp2crate` is enabled.
//! `fp2.rs` is unchanged and builds GF(p²) on top of this `Fp`.

#![allow(clippy::needless_range_loop)]

use core::fmt;

#[cfg(any(feature = "lvl3", feature = "lvl5"))]
compile_error!("gf-fp2crate currently supports lvl1 only");

#[cfg(feature = "sign")]
compile_error!(
    "gf-fp2crate is verify-only for now: sign_data precomp constants are stored in the \
     modarith internal representation and need re-encoding for this backend"
);

const MODULUS: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0x04FF_FFFF_FFFF_FFFF,
];

::fp2::define_fp_core!(typename = FpInner, modulus = MODULUS,);

// (p − 3)/4, the progenitor exponent shared by inv and sqrt.
const EXP_3DIV4: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0x013F_FFFF_FFFF_FFFF,
];

/// Word type used by the C reference (`digit_t` for RADIX_64).
pub type Digit = u64;
pub const RADIX: u32 = 64;

pub const NWORDS_FIELD: usize = FpInner::N;
pub const NWORDS_ORDER: usize = 4;
pub const BITS: u32 = 256;
pub const LOG2P: u32 = 8;
pub const FP_ENCODED_BYTES: usize = FpInner::ENCODED_LENGTH;

/// An element of GF(p), Montgomery form, saturated 4×64-bit limbs (`fp2` crate).
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Fp(pub FpInner);

pub const ZERO: Fp = Fp(FpInner::ZERO);
pub const ONE: Fp = Fp(FpInner::ONE);
pub const MINUS_ONE: Fp = Fp(FpInner::MINUS_ONE);

impl PartialEq for Fp {
    fn eq(&self, other: &Self) -> bool {
        self.0.equals(&other.0) == u32::MAX
    }
}
impl Eq for Fp {}

impl Default for Fp {
    fn default() -> Self {
        ZERO
    }
}

// ---- operator wiring (fp2 crate already provides ops on FpInner) ----

macro_rules! fp_binop {
    ($t:ty, $tr:ident, $f:ident, $atr:ident, $af:ident) => {
        impl core::ops::$tr for $t {
            type Output = $t;
            #[inline]
            fn $f(mut self, rhs: $t) -> $t {
                core::ops::$atr::$af(&mut self, &rhs);
                self
            }
        }
        impl core::ops::$tr<&$t> for $t {
            type Output = $t;
            #[inline]
            fn $f(mut self, rhs: &$t) -> $t {
                core::ops::$atr::$af(&mut self, rhs);
                self
            }
        }
        impl core::ops::$tr<$t> for &$t {
            type Output = $t;
            #[inline]
            fn $f(self, rhs: $t) -> $t {
                let mut r = *self;
                core::ops::$atr::$af(&mut r, &rhs);
                r
            }
        }
        impl core::ops::$tr<&$t> for &$t {
            type Output = $t;
            #[inline]
            fn $f(self, rhs: &$t) -> $t {
                let mut r = *self;
                core::ops::$atr::$af(&mut r, rhs);
                r
            }
        }
        impl core::ops::$tr<$t> for &mut $t {
            type Output = $t;
            #[inline]
            fn $f(self, rhs: $t) -> $t {
                let mut r = *self;
                core::ops::$atr::$af(&mut r, &rhs);
                r
            }
        }
        impl core::ops::$atr<$t> for $t {
            #[inline]
            fn $af(&mut self, rhs: $t) {
                core::ops::$atr::$af(self, &rhs);
            }
        }
    };
}

macro_rules! fp_ops {
    ($t:ty) => {
        fp_binop!($t, Add, add, AddAssign, add_assign);
        fp_binop!($t, Sub, sub, SubAssign, sub_assign);
        fp_binop!($t, Mul, mul, MulAssign, mul_assign);
        impl core::ops::Neg for $t {
            type Output = $t;
            #[inline]
            fn neg(mut self) -> $t {
                self.neg_ip();
                self
            }
        }
        impl core::ops::Neg for &$t {
            type Output = $t;
            #[inline]
            fn neg(self) -> $t {
                let mut r = *self;
                r.neg_ip();
                r
            }
        }
    };
}
pub(crate) use fp_binop;
pub(crate) use fp_ops;

impl core::ops::AddAssign<&Fp> for Fp {
    #[inline]
    fn add_assign(&mut self, rhs: &Fp) {
        self.0 += rhs.0;
    }
}
impl core::ops::SubAssign<&Fp> for Fp {
    #[inline]
    fn sub_assign(&mut self, rhs: &Fp) {
        self.0 -= rhs.0;
    }
}
impl core::ops::MulAssign<&Fp> for Fp {
    #[inline]
    fn mul_assign(&mut self, rhs: &Fp) {
        self.0 *= rhs.0;
    }
}
fp_ops!(Fp);

impl Fp {
    pub const ZERO: Self = ZERO;
    pub const ONE: Self = ONE;
    pub const MINUS_ONE: Self = MINUS_ONE;

    /// Construct from canonical little-endian bytes (no range check).
    pub const fn from_le_bytes_unchecked(buf: &[u8; FP_ENCODED_BYTES]) -> Self {
        Fp(FpInner::const_decode_no_check(buf))
    }

    #[inline]
    pub fn from_small(val: Digit) -> Self {
        Fp(FpInner::from(val))
    }

    // ---- arithmetic ----

    #[inline]
    #[must_use]
    pub fn square(self) -> Self {
        Fp(self.0.square())
    }
    #[inline]
    pub fn square_ip(&mut self) {
        self.0.set_square();
    }
    #[inline]
    #[must_use]
    pub fn dbl(self) -> Self {
        Fp(self.0.mul2())
    }
    #[inline]
    pub fn dbl_ip(&mut self) {
        self.0.set_mul2();
    }
    #[inline]
    pub fn neg_ip(&mut self) {
        self.0.set_neg();
    }
    #[inline]
    #[must_use]
    pub fn half(self) -> Self {
        Fp(self.0.half())
    }
    #[inline]
    #[must_use]
    pub fn div3(self) -> Self {
        static INV3: std::sync::OnceLock<FpInner> = std::sync::OnceLock::new();
        Fp(self.0 * *INV3.get_or_init(|| FpInner::THREE.invert()))
    }
    #[inline]
    #[must_use]
    pub fn mul_small(self, val: u32) -> Self {
        Fp(self.0.mul_small(val as i32))
    }
    #[inline]
    #[must_use]
    pub fn inv(self) -> Self {
        Fp(self.0.invert())
    }
    /// `self^((p-3)/4)`: the progenitor used for both inversion and sqrt.
    #[inline]
    #[must_use]
    pub fn exp_3div4(self) -> Self {
        Fp(self.0.pow_pubexp(&EXP_3DIV4))
    }
    /// √self, assuming self is a square; result is undefined otherwise.
    /// Returns `self^((p+1)/4)` (matches the modarith backend's `modsqrt`).
    #[inline]
    #[must_use]
    pub fn sqrt(self) -> Self {
        // fp2 crate's sqrt zeroes the result if non-square; we want the
        // raw exponentiation to match the modarith backend on undefined input.
        let h = self.exp_3div4();
        h * self
    }

    // ---- predicates ----

    #[inline]
    pub fn is_zero_ct(&self) -> u32 {
        self.0.is_zero()
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero() == u32::MAX
    }
    #[inline]
    pub fn is_equal_ct(&self, other: &Self) -> u32 {
        self.0.equals(&other.0)
    }
    #[inline]
    pub fn is_square_ct(&self) -> u32 {
        // legendre() returns 0, 1, or -1; map {0,1}→all-ones, -1→0 (matches modarith modqr).
        !((self.0.legendre() >> 1) as u32)
    }
    #[inline]
    pub fn is_square(&self) -> bool {
        self.is_square_ct() == u32::MAX
    }

    // ---- conditional / encode ----

    #[inline]
    pub fn select(a0: &Self, a1: &Self, ctl: u32) -> Self {
        Fp(FpInner::select(&a0.0, &a1.0, ctl))
    }
    #[inline]
    pub fn cswap(a: &mut Self, b: &mut Self, ctl: u32) {
        FpInner::cond_swap(&mut a.0, &mut b.0, ctl);
    }
    pub fn encode(&self) -> [u8; FP_ENCODED_BYTES] {
        self.0.encode()
    }
    pub fn try_decode(src: &[u8]) -> Option<Self> {
        let (v, ok) = FpInner::decode(src);
        (ok == u32::MAX).then_some(Fp(v))
    }
    pub fn decode_reduce(src: &[u8]) -> Self {
        Fp(FpInner::decode_reduce(src))
    }
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp(0x")?;
        for b in self.encode().iter().rev() {
            write!(f, "{b:02x}")?;
        }
        write!(f, ")")
    }
}
