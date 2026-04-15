//! GF(p) backend via the vendored `GF5_248` type (4-limb Montgomery,
//! p = 5·2²⁴⁸ − 1, with hand-written x86_64 BMI2/ADX asm including a
//! dedicated `fp_sqr_asm`).
//!
//! Default GF(p) backend at lvl1. `fp2.rs` is unchanged.

use core::fmt;

#[path = "gf5_248_vendored/mod.rs"]
mod gf5_248;
pub use gf5_248::GF5_248 as FpInner;

/// Word type used by the C reference (`digit_t` for RADIX_64).
pub type Digit = u64;
pub const RADIX: u32 = 64;

pub const NWORDS_FIELD: usize = FpInner::N;
pub const NWORDS_ORDER: usize = 4;
pub const BITS: u32 = 256;
pub const LOG2P: u32 = 8;
pub const FP_ENCODED_BYTES: usize = FpInner::ENCODED_LENGTH;

/// An element of GF(p), Montgomery form, saturated 4×64-bit limbs (vendored asm).
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

// ---- operator wiring ----

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
    /// This backend has a fused `a·b ± c·d` (one Montgomery reduction).
    pub const HAS_FUSED_SUMPROD: bool = true;

    /// `a·b + c·d` with one Montgomery reduction.
    #[inline]
    pub fn mul_add(a: &Self, b: &Self, c: &Self, d: &Self) -> Self {
        Fp(FpInner::sum_of_products(&a.0, &b.0, &c.0, &d.0))
    }
    /// `a·b − c·d` with one Montgomery reduction.
    #[inline]
    pub fn mul_sub(a: &Self, b: &Self, c: &Self, d: &Self) -> Self {
        Fp(FpInner::difference_of_products(&a.0, &b.0, &c.0, &d.0))
    }
    /// `(re, im)` of the GF(p²) product `a × b`, reading both operands as
    /// the `Fp2` memory layout (`[re.limbs, im.limbs]`). Zero-copy on asm.
    #[cfg(gf5_248_asm)]
    #[inline]
    pub fn fp2_mul_kernel(a: &[u64; 8], b: &[u64; 8]) -> (Self, Self) {
        let (re, im) = FpInner::fp2_mul_kernel(a, b);
        (Fp(re), Fp(im))
    }
    /// `a + b` without conditional reduction; result may exceed 2²⁵¹ but is a
    /// valid `Mul`/`mul_add`/`mul_sub` operand. Used by `Fp2::square`.
    #[inline]
    pub fn add_noreduce(a: &Self, b: &Self) -> Self {
        Fp(FpInner::add_noreduce(&a.0, &b.0))
    }
    /// `a − b + 2p` without conditional reduction. Used by `Fp2::square`.
    #[inline]
    pub fn sub_2p_noreduce(a: &Self, b: &Self) -> Self {
        Fp(FpInner::sub_2p_noreduce(&a.0, &b.0))
    }

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
    /// (p-3)/4 = 5·2²⁴⁶ − 1.
    #[must_use]
    pub fn exp_3div4(self) -> Self {
        let z = self.0;
        let z2 = z.square();
        let z4 = z2.square();
        let e2 = z2 * z; // 2²−1
        let e4 = e2.n_square(2) * e2; // 2⁴−1
        let e8 = e4.n_square(4) * e4;
        let e16 = e8.n_square(8) * e8;
        let e32 = e16.n_square(16) * e16;
        let e64 = e32.n_square(32) * e32;
        let e128 = e64.n_square(64) * e64;
        let e192 = e128.n_square(64) * e64;
        let e224 = e192.n_square(32) * e32;
        let e240 = e224.n_square(16) * e16;
        let e244 = e240.n_square(4) * e4;
        let w = e244.n_square(2) * e2; // 2²⁴⁶−1
        Fp(w.n_square(2) * w * z4) // 5·(2²⁴⁶−1) + 4 = (p−3)/4
    }
    /// √self, assuming self is a square; result is undefined otherwise.
    /// Returns `self^((p+1)/4)` (matches the modarith backend's `modsqrt`).
    #[inline]
    #[must_use]
    pub fn sqrt(self) -> Self {
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
        // legendre() returns 0, 1, or -1; map {0,1}→all-ones, -1→0.
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

#[cfg(all(test, gf5_248_asm))]
mod tests {
    use super::*;
    use crate::test_util::Prng;

    /// `fp_sqr` (inline asm) must be limb-identical to `fp_mul(x, x)` (inline
    /// asm) and to the extern-C `.S` `fp_sqr_asm`/`fp_mul_asm` for all inputs
    /// in the working range (< 2²⁵²).
    #[test]
    fn sqr_matches_mul() {
        let mut prng = Prng(0x5248_5248);
        for _ in 0..100_000 {
            let x = prng.fp();
            let mut s = x.0;
            s.set_square();
            let mut m = x.0;
            m.set_square_via_mul();
            let mut se = x.0;
            se.set_square_via_extern();
            let mut me = x.0;
            me.set_mul_via_extern(&x.0);
            assert_eq!(s.encode(), m.encode());
            assert_eq!(s.encode(), se.encode());
            assert_eq!(s.encode(), me.encode());
        }
        // Edge cases.
        for x in [Fp::ZERO, Fp::ONE, Fp::MINUS_ONE, Fp::from_small(u64::MAX)] {
            let mut s = x.0;
            s.set_square();
            let mut m = x.0;
            m.set_square_via_mul();
            assert_eq!(s.encode(), m.encode());
        }
    }

    /// `mul_add`/`mul_sub` (inline asm) must equal `a*b ± c*d` via separate
    /// muls, and equal the extern-C `.S` `fp_sop_asm`/`fp_dop_asm`.
    #[test]
    fn sumprod_matches_separate() {
        let mut prng = Prng(0x50_50_50_50);
        for _ in 0..100_000 {
            let (a, b, c, d) = (prng.fp(), prng.fp(), prng.fp(), prng.fp());
            let add_ref = (a * b + c * d).encode();
            let sub_ref = (a * b - c * d).encode();
            assert_eq!(Fp::mul_add(&a, &b, &c, &d).encode(), add_ref);
            assert_eq!(Fp::mul_sub(&a, &b, &c, &d).encode(), sub_ref);
            let se = FpInner::sum_of_products_via_extern(&a.0, &b.0, &c.0, &d.0);
            let de = FpInner::difference_of_products_via_extern(&a.0, &b.0, &c.0, &d.0);
            assert_eq!(se.encode(), add_ref);
            assert_eq!(de.encode(), sub_ref);
        }
        // Edge cases including max-magnitude inputs.
        let m = Fp::MINUS_ONE; // largest canonical value
        for &(a, b, c, d) in &[
            (Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO),
            (Fp::ONE, Fp::ONE, Fp::ONE, Fp::ONE),
            (m, m, m, m),
            (m, m, Fp::ZERO, Fp::ZERO),
            (Fp::ZERO, Fp::ZERO, m, m),
            (m, Fp::ONE, Fp::ONE, m),
        ] {
            assert_eq!(
                Fp::mul_add(&a, &b, &c, &d).encode(),
                (a * b + c * d).encode()
            );
            assert_eq!(
                Fp::mul_sub(&a, &b, &c, &d).encode(),
                (a * b - c * d).encode()
            );
        }
    }

    #[test]
    fn sqr_algebraic() {
        let mut prng = Prng(0x1111);
        for _ in 0..1000 {
            let x = prng.fp();
            assert_eq!((x * x).encode(), x.square().encode());
        }
    }

    #[test]
    fn exp_3div4_is_progenitor() {
        // x^((p-3)/4) · x² · x = x^((p+1)/4)² · x⁻¹ · ... — simplest check:
        // (x^((p-3)/4))² · x³ = x^((p-3)/2 + 3) = x^((p+3)/2) = x · x^((p+1)/2) = x · legendre(x).
        // For squares: exp_3div4(x)·x = sqrt(x), so (exp_3div4(x)·x)² = x.
        let mut prng = Prng(0x2222);
        for _ in 0..200 {
            let x = prng.fp().square(); // ensure square
            let r = x.exp_3div4() * x;
            assert_eq!(r.square().encode(), x.encode());
        }
    }
}
