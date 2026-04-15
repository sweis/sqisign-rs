//! GF(p²) arithmetic, modulo X² + 1. Port of `lvlx/fp2.c`.

use super::fp::{decode_masked, fp_binop, fp_ops, Digit, Fp, FP_ENCODED_BYTES, RADIX};
use core::fmt;
use core::ops::{AddAssign, MulAssign, SubAssign};

/// `FP2_ENCODED_BYTES`.
pub const FP2_ENCODED_BYTES: usize = 2 * FP_ENCODED_BYTES;

/// An element of GF(p²) ≅ GF(p)[i]/(i²+1).
#[derive(Clone, Copy, Default, Eq, PartialEq)]
#[repr(C)]
pub struct Fp2 {
    pub re: Fp,
    pub im: Fp,
}

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} + i*{:?}", self.re, self.im)
    }
}

impl Fp2 {
    pub const ZERO: Self = Self {
        re: Fp::ZERO,
        im: Fp::ZERO,
    };
    pub const ONE: Self = Self {
        re: Fp::ONE,
        im: Fp::ZERO,
    };

    #[inline]
    pub fn from_small(val: Digit) -> Self {
        Self {
            re: Fp::from_small(val),
            im: Fp::ZERO,
        }
    }

    // ---- arithmetic ----

    /// Self². Karatsuba: re = (a+b)(a-b), im = 2ab.
    #[inline]
    #[must_use]
    pub fn square(self) -> Self {
        let sum = self.re + self.im;
        let diff = self.re - self.im;
        Self {
            re: sum * diff,
            im: (self.re * self.im).dbl(),
        }
    }
    #[inline]
    pub fn square_ip(&mut self) {
        *self = self.square();
    }

    #[inline]
    #[must_use]
    pub fn dbl(self) -> Self {
        Self {
            re: self.re.dbl(),
            im: self.im.dbl(),
        }
    }
    #[inline]
    pub fn dbl_ip(&mut self) {
        self.re.dbl_ip();
        self.im.dbl_ip();
    }

    #[inline]
    pub fn neg_ip(&mut self) {
        self.re.neg_ip();
        self.im.neg_ip();
    }

    #[inline]
    #[must_use]
    pub fn half(self) -> Self {
        Self {
            re: self.re.half(),
            im: self.im.half(),
        }
    }

    #[inline]
    #[must_use]
    pub fn mul_small(self, n: u32) -> Self {
        Self {
            re: self.re.mul_small(n),
            im: self.im.mul_small(n),
        }
    }

    /// 1/self; returns 0 for input 0.
    #[must_use]
    pub fn inv(self) -> Self {
        let n = (self.re.square() + self.im.square()).inv();
        Self {
            re: self.re * n,
            im: -(self.im * n),
        }
    }

    /// Frobenius / complex conjugation: self^p.
    #[inline]
    #[must_use]
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Variable-time square-and-multiply: self^exp, exp little-endian limbs.
    #[must_use]
    pub fn pow_vartime(self, exp: &[Digit]) -> Self {
        let mut acc = self;
        let mut out = Fp2::ONE;
        for &word in exp {
            for i in 0..RADIX {
                if (word >> i) & 1 == 1 {
                    out *= acc;
                }
                acc.square_ip();
            }
        }
        out
    }

    /// Square root in GF(p²) following Aardal et al. (ePrint 2024/1563),
    /// canonicalized so the result has even real part (or, if re==0, even imag).
    /// Result is undefined when self is not a square.
    #[must_use]
    pub fn sqrt(self) -> Self {
        let s = (self.re.square() + self.im.square()).sqrt();
        let mut x0 = Fp::select(&s, &self.re, self.im.is_zero_ct());
        x0 += self.re;
        let t0 = x0.dbl();

        let mut x1 = t0.exp_3div4();

        x0 *= x1;
        x1 *= self.im;
        let t1 = x0.dbl().square();
        let f = (t0 - t1).is_zero_ct();
        let t0 = Fp::select(&x1, &x0, f);
        let t1 = Fp::select(&(-x0), &x1, f);

        let t0_is_zero = t0.is_zero_ct();
        let t0_is_odd = ((t0.encode()[0] & 1) as u32).wrapping_neg();
        let t1_is_odd = ((t1.encode()[0] & 1) as u32).wrapping_neg();

        let negate = t0_is_odd | (t0_is_zero & t1_is_odd);
        Self {
            re: Fp::select(&t0, &(-t0), negate),
            im: Fp::select(&t1, &(-t1), negate),
        }
    }

    /// √self if self is a square, else `None`. Constant-time in the value of self.
    #[must_use]
    pub fn sqrt_verify(self) -> Option<Self> {
        let r = self.sqrt();
        (r.square().is_equal_ct(&self) != 0).then_some(r)
    }

    /// `self + 1`.
    #[inline]
    #[must_use]
    pub fn add_one(self) -> Self {
        Self {
            re: self.re + Fp::ONE,
            im: self.im,
        }
    }

    // ---- predicates ----
    // _ct variants return all-ones / all-zeros u32 masks for constant-time select.

    #[inline]
    pub fn is_zero_ct(&self) -> u32 {
        self.re.is_zero_ct() & self.im.is_zero_ct()
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.is_zero_ct() != 0
    }
    #[inline]
    pub fn is_one_ct(&self) -> u32 {
        self.re.is_equal_ct(&Fp::ONE) & self.im.is_zero_ct()
    }
    #[inline]
    pub fn is_one(&self) -> bool {
        self.is_one_ct() != 0
    }
    #[inline]
    pub fn is_equal_ct(&self, other: &Self) -> u32 {
        self.re.is_equal_ct(&other.re) & self.im.is_equal_ct(&other.im)
    }
    #[inline]
    pub fn is_square_ct(&self) -> u32 {
        (self.re.square() + self.im.square()).is_square_ct()
    }
    #[inline]
    pub fn is_square(&self) -> bool {
        self.is_square_ct() != 0
    }

    // ---- conditional / encode ----

    #[inline]
    pub fn select(a0: &Self, a1: &Self, ctl: u32) -> Self {
        Self {
            re: Fp::select(&a0.re, &a1.re, ctl),
            im: Fp::select(&a0.im, &a1.im, ctl),
        }
    }
    #[inline]
    pub fn cswap(a: &mut Self, b: &mut Self, ctl: u32) {
        Fp::cswap(&mut a.re, &mut b.re, ctl);
        Fp::cswap(&mut a.im, &mut b.im, ctl);
    }

    #[inline]
    pub fn encode(&self) -> [u8; FP2_ENCODED_BYTES] {
        let mut out = [0u8; FP2_ENCODED_BYTES];
        out[..FP_ENCODED_BYTES].copy_from_slice(&self.re.encode());
        out[FP_ENCODED_BYTES..].copy_from_slice(&self.im.encode());
        out
    }

    /// Decode canonical little-endian bytes; returns `None` if either half ≥ p.
    #[inline]
    pub fn try_decode(src: &[u8]) -> Option<Self> {
        let mut x = Self::ZERO;
        let re = decode_masked(&mut x.re, &src[..FP_ENCODED_BYTES]);
        let im = decode_masked(&mut x.im, &src[FP_ENCODED_BYTES..2 * FP_ENCODED_BYTES]);
        ((re & im) == 0xFFFF_FFFF).then_some(x)
    }
}

impl AddAssign<&Fp2> for Fp2 {
    #[inline]
    fn add_assign(&mut self, rhs: &Fp2) {
        self.re += &rhs.re;
        self.im += &rhs.im;
    }
}
impl SubAssign<&Fp2> for Fp2 {
    #[inline]
    fn sub_assign(&mut self, rhs: &Fp2) {
        self.re -= &rhs.re;
        self.im -= &rhs.im;
    }
}
/// (a + bi)(c + di), via Karatsuba.
impl MulAssign<&Fp2> for Fp2 {
    #[inline]
    fn mul_assign(&mut self, rhs: &Fp2) {
        let t0 = (self.re + self.im) * (rhs.re + rhs.im);
        let t1 = self.im * rhs.im;
        self.re *= &rhs.re;
        self.im = t0 - t1 - self.re;
        self.re -= t1;
    }
}
fp_ops!(Fp2);

/// In-place batched inversion (Montgomery's trick).
pub fn fp2_batched_inv(x: &mut [Fp2]) {
    let len = x.len();
    if len == 0 {
        return;
    }
    let mut t1 = vec![Fp2::ZERO; len];
    t1[0] = x[0];
    for i in 1..len {
        t1[i] = t1[i - 1] * x[i];
    }
    let mut acc = t1[len - 1].inv();
    for i in (1..len).rev() {
        let r = acc * t1[i - 1];
        acc *= x[i];
        x[i] = r;
    }
    x[0] = acc;
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    use crate::test_util::assert_hex;
    use crate::test_util::Prng;

    fn fp2_random(prng: &mut Prng) -> Fp2 {
        prng.fp2()
    }

    const ITERS: usize = 300;

    #[test]
    fn addition_laws() {
        let mut prng = Prng(0x10);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            assert_eq!((a + b) + c, a + (b + c));
            assert_eq!(a + b, b + a);
            assert!((a + (-a)).is_zero());
            assert_eq!(a + Fp2::ONE, a.add_one());
        }
    }

    #[test]
    fn subtraction_laws() {
        let mut prng = Prng(0x11);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            assert_eq!((a - b) - c, a - (b + c));
            assert_eq!(a - b, -(b - a));
            assert!((a - a).is_zero());
        }
    }

    #[test]
    fn multiplication_laws() {
        let mut prng = Prng(0x12);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a * b, b * a);
            assert_eq!(a * Fp2::ONE, a);
            assert!((a * Fp2::ZERO).is_zero());
        }
    }

    #[test]
    fn mul_matches_schoolbook() {
        // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        let mut prng = Prng(0x13);
        for _ in 0..ITERS {
            let y = fp2_random(&mut prng);
            let z = fp2_random(&mut prng);
            let x = y * z;
            assert_eq!(x.re, y.re * z.re - y.im * z.im);
            assert_eq!(x.im, y.re * z.im + y.im * z.re);
        }
    }

    #[test]
    fn mul_small() {
        let mut prng = Prng(0x14);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let val = (prng.next() as u32) & 0x7FFF_FFFF;
            assert_eq!(a.mul_small(val), a * Fp2::from_small(val as Digit));
        }
    }

    #[test]
    fn assign_ops_match_value_ops() {
        let mut prng = Prng(0x21);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let mut r = a;
            r *= b;
            assert_eq!(r, a * b);
            let mut r = a;
            r.square_ip();
            assert_eq!(r, a.square());
            let mut r = a;
            r += b;
            assert_eq!(r, a + b);
            let mut r = a;
            r -= b;
            assert_eq!(r, a - b);
            let mut r = a;
            r.neg_ip();
            assert_eq!(r, -a);
            let mut r = a;
            r.dbl_ip();
            assert_eq!(r, a + a);
            assert_eq!(r, a.dbl());
        }
    }

    #[test]
    fn squaring() {
        let mut prng = Prng(0x15);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn inversion() {
        let mut prng = Prng(0x16);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            assert_eq!(a * a.inv(), Fp2::ONE);
        }
        assert!(Fp2::ZERO.inv().is_zero());
    }

    #[test]
    fn sqrt_and_is_square() {
        let mut prng = Prng(0x17);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let c = a.square();
            assert!(c.is_square());
            assert_eq!(c.sqrt_verify().map(|r| r.square()), Some(c));
            let r = c.sqrt();
            assert!(r == a || r == -a);
        }
    }

    #[test]
    fn batched_inv() {
        let mut prng = Prng(0x18);
        let mut xs: Vec<Fp2> = (0..8).map(|_| fp2_random(&mut prng)).collect();
        let orig = xs.clone();
        fp2_batched_inv(&mut xs);
        for (x, o) in xs.iter().zip(orig.iter()) {
            assert_eq!(*x * *o, Fp2::ONE);
        }
    }

    /// Golden vectors extracted from the C reference implementation
    /// (libsqisign_gf_lvl1.a, RADIX_64, ref build).
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    #[test]
    fn golden_c_vectors() {
        let x = Fp2 {
            re: Fp::from_small(3),
            im: Fp::from_small(4),
        };
        let y = x.square();
        assert_hex(
            &y.encode(),
            "f8ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04\
             1800000000000000000000000000000000000000000000000000000000000000",
        );
        let s = y.sqrt();
        assert_hex(
            &s.encode(),
            "fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04\
             fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
        );
        assert_hex(
            &s.inv().encode(),
            "6666666666666666666666666666666666666666666666666666666666666601\
             cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc04",
        );
    }

    #[test]
    fn encode_decode_roundtrip() {
        let mut prng = Prng(0x19);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            assert_eq!(Fp2::try_decode(&a.encode()), Some(a));
        }
    }
}
