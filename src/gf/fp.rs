//! GF(p) arithmetic for the active SQIsign prime.
//!
//! This file is the level-agnostic layer. The unrolled multiply/square,
//! the progenitor addition chain, and all per-prime constants live in the
//! `super::backend` module (one of `fp_p5248`/`fp_p65376`/`fp_p27500`),
//! selected via cargo feature.
//!
//! All values are kept in Montgomery form internally; `encode`/`decode`
//! convert between canonical little-endian bytes and the internal form.

#![allow(clippy::needless_range_loop)]

use core::fmt;

use super::backend::{
    modmul, modpro, modsqr, partial_reduce, LIMB_BITS, MASK, NBYTES, NLIMBS, NRES_C, PR_WORDS,
    P_TOP, R2, THREE_INV, TWO_INV, TWO_P_TOP,
};
pub use super::backend::{
    BITS, FP_ENCODED_BYTES, LOG2P, MINUS_ONE, NWORDS_FIELD, NWORDS_ORDER, ONE, ZERO,
};

/// Word type used by the C reference (`digit_t` for RADIX_64).
pub type Digit = u64;
/// `RADIX` from `tutil.h`.
pub const RADIX: u32 = 64;

/// An element of GF(p), in Montgomery form, unsaturated `NLIMBS` × `LIMB_BITS`-bit limbs.
#[derive(Clone, Copy)]
#[repr(C)]
pub struct Fp(pub [u64; NWORDS_FIELD]);

// The internal representation is not unique; compare canonical forms.
impl PartialEq for Fp {
    fn eq(&self, other: &Self) -> bool {
        modcmp(&self.0, &other.0) != 0
    }
}
impl Eq for Fp {}

impl Default for Fp {
    fn default() -> Self {
        Fp([0; NWORDS_FIELD])
    }
}

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

impl Fp {
    pub const ZERO: Self = ZERO;
    pub const ONE: Self = ONE;
    pub const MINUS_ONE: Self = MINUS_ONE;
    /// No fused `a·b ± c·d` here; `Fp2::mul` uses 3-mul Karatsuba instead.
    pub const HAS_FUSED_SUMPROD: bool = false;

    /// `a·b + c·d`. (No fused path on this backend.)
    #[inline]
    pub fn mul_add(a: &Self, b: &Self, c: &Self, d: &Self) -> Self {
        *a * b + *c * d
    }
    /// `a·b − c·d`. (No fused path on this backend.)
    #[inline]
    pub fn mul_sub(a: &Self, b: &Self, c: &Self, d: &Self) -> Self {
        *a * b - *c * d
    }
    #[inline]
    pub fn add_noreduce(a: &Self, b: &Self) -> Self {
        *a + b
    }
    #[inline]
    pub fn sub_2p_noreduce(a: &Self, b: &Self) -> Self {
        *a - b
    }

    #[inline]
    pub fn from_small(val: Digit) -> Self {
        let mut x = Self::ZERO;
        modint(val as i32, &mut x.0);
        x
    }

    // ---- arithmetic ----

    #[inline]
    #[must_use]
    pub fn square(self) -> Self {
        let mut r = Self::ZERO;
        modsqr(&self.0, &mut r.0);
        r
    }
    #[inline]
    pub fn square_ip(&mut self) {
        let a = self.0;
        modsqr(&a, &mut self.0);
    }
    #[inline]
    #[must_use]
    pub fn dbl(self) -> Self {
        let mut r = Self::ZERO;
        modadd(&self.0, &self.0, &mut r.0);
        r
    }
    #[inline]
    pub fn dbl_ip(&mut self) {
        let a = self.0;
        modadd(&a, &a, &mut self.0);
    }
    #[inline]
    pub fn neg_ip(&mut self) {
        let a = self.0;
        modneg(&a, &mut self.0);
    }
    #[inline]
    #[must_use]
    pub fn half(self) -> Self {
        let mut r = Self::ZERO;
        modmul(&TWO_INV.0, &self.0, &mut r.0);
        r
    }
    #[inline]
    #[must_use]
    pub fn div3(self) -> Self {
        let mut r = Self::ZERO;
        modmul(&THREE_INV.0, &self.0, &mut r.0);
        r
    }
    #[inline]
    #[must_use]
    pub fn mul_small(self, val: u32) -> Self {
        let mut r = Self::ZERO;
        modmli(&self.0, val as i32, &mut r.0);
        r
    }
    #[inline]
    #[must_use]
    pub fn inv(self) -> Self {
        let mut r = Self::ZERO;
        modinv(&self.0, None, &mut r.0);
        r
    }
    /// `self^((p-3)/4)`: the progenitor used for both inversion and sqrt.
    #[inline]
    #[must_use]
    pub fn exp_3div4(self) -> Self {
        let mut r = Self::ZERO;
        modpro(&self.0, &mut r.0);
        r
    }
    /// √self, assuming self is a square; result is undefined otherwise.
    #[inline]
    #[must_use]
    pub fn sqrt(self) -> Self {
        let mut r = Self::ZERO;
        modsqrt(&self.0, None, &mut r.0);
        r
    }

    // ---- predicates ----
    // _ct variants return all-ones / all-zeros u32 masks for constant-time select.

    #[inline]
    pub fn is_zero_ct(&self) -> u32 {
        (modis0(&self.0) as u32).wrapping_neg()
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        modis0(&self.0) != 0
    }
    #[inline]
    pub fn is_equal_ct(&self, other: &Self) -> u32 {
        (modcmp(&self.0, &other.0) as u32).wrapping_neg()
    }
    #[inline]
    pub fn is_square_ct(&self) -> u32 {
        (modqr(None, &self.0) as u32).wrapping_neg()
    }
    #[inline]
    pub fn is_square(&self) -> bool {
        modqr(None, &self.0) != 0
    }

    // ---- conditional / encode ----

    /// Constant-time select: returns `a0` if `ctl == 0`, `a1` if `ctl == 0xFFFFFFFF`.
    #[inline]
    pub fn select(a0: &Self, a1: &Self, ctl: u32) -> Self {
        let cw = (ctl as i32) as i64 as u64;
        let mut d = Self::ZERO;
        for i in 0..NWORDS_FIELD {
            d.0[i] = a0.0[i] ^ (cw & (a0.0[i] ^ a1.0[i]));
        }
        d
    }
    /// Constant-time conditional swap.
    #[inline]
    pub fn cswap(a: &mut Self, b: &mut Self, ctl: u32) {
        modcsw((ctl & 1) as i32, &mut a.0, &mut b.0);
    }

    pub fn encode(&self) -> [u8; FP_ENCODED_BYTES] {
        let mut dst = [0u8; FP_ENCODED_BYTES];
        let mut c = [0u64; NLIMBS];
        redc(&self.0, &mut c);
        for b in dst.iter_mut().take(NBYTES) {
            *b = (c[0] & 0xff) as u8;
            let _ = modshr(8, &mut c);
        }
        dst
    }

    /// Decode canonical little-endian bytes; returns `None` if value ≥ p.
    pub fn try_decode(src: &[u8]) -> Option<Self> {
        let mut d = Self::ZERO;
        (decode_masked(&mut d, src) == 0xFFFF_FFFF).then_some(d)
    }

    /// Decode an arbitrary-length little-endian byte string and reduce mod p.
    pub fn decode_reduce(src: &[u8]) -> Self {
        let mut d = Fp::ZERO;
        let mut len = src.len();
        if len == 0 {
            return d;
        }
        let mut tmp = [0u8; NBYTES];

        let rem = len % NBYTES;
        if rem != 0 {
            let k = len - rem;
            tmp[..rem].copy_from_slice(&src[k..len]);
            tmp[rem..].fill(0);
            decode_masked(&mut d, &tmp);
            len = k;
        }

        while len > 0 {
            d *= R2;
            len -= NBYTES;
            let mut t = [0u64; PR_WORDS];
            for j in 0..PR_WORDS {
                t[j] = u64::from_le_bytes(src[len + 8 * j..len + 8 * j + 8].try_into().unwrap());
            }
            let tin = t;
            partial_reduce(&mut t, &tin);
            for j in 0..PR_WORDS {
                tmp[8 * j..8 * j + 8].copy_from_slice(&t[j].to_le_bytes());
            }
            let mut a = Fp::ZERO;
            decode_masked(&mut a, &tmp);
            d += a;
        }
        d
    }
}

impl core::ops::AddAssign<&Fp> for Fp {
    #[inline]
    fn add_assign(&mut self, rhs: &Fp) {
        let a = self.0;
        modadd(&a, &rhs.0, &mut self.0);
    }
}
impl core::ops::SubAssign<&Fp> for Fp {
    #[inline]
    fn sub_assign(&mut self, rhs: &Fp) {
        let a = self.0;
        modsub(&a, &rhs.0, &mut self.0);
    }
}
impl core::ops::MulAssign<&Fp> for Fp {
    #[inline]
    fn mul_assign(&mut self, rhs: &Fp) {
        let a = self.0;
        modmul(&a, &rhs.0, &mut self.0);
    }
}
fp_ops!(Fp);

#[inline]
fn decode_masked(d: &mut Fp, src: &[u8]) -> u32 {
    debug_assert!(src.len() >= FP_ENCODED_BYTES);
    d.0 = [0u64; NLIMBS];
    for i in (0..NBYTES).rev() {
        modshl(8, &mut d.0);
        d.0[0] = d.0[0].wrapping_add(src[i] as u64);
    }
    let res = (modfsb(&mut d.0) as u64).wrapping_neg();
    let t = d.0;
    nres(&t, &mut d.0);
    for i in 0..NLIMBS {
        d.0[i] &= res;
    }
    res as u32
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

// ===========================================================================
// Internal modarith primitives shared across primes (loop-based; the unrolled
// modmul/modsqr/modpro live in the backend).
// ===========================================================================

#[inline]
fn prop(n: &mut [u64; NLIMBS]) -> u64 {
    let mut carry: i64 = (n[0] as i64) >> LIMB_BITS;
    n[0] &= MASK;
    for i in 1..(NLIMBS - 1) {
        carry = carry.wrapping_add(n[i] as i64);
        n[i] = (carry as u64) & MASK;
        carry >>= LIMB_BITS;
    }
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(carry as u64);
    ((n[NLIMBS - 1] >> 1) >> 62).wrapping_neg()
}

#[inline]
fn flatten(n: &mut [u64; NLIMBS]) -> i32 {
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(1u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(P_TOP & carry);
    let _ = prop(n);
    (carry & 1) as i32
}

#[inline]
fn modfsb(n: &mut [u64; NLIMBS]) -> i32 {
    n[0] = n[0].wrapping_add(1);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_sub(P_TOP);
    flatten(n)
}

#[inline]
fn modadd(a: &[u64; NLIMBS], b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = a[i].wrapping_add(b[i]);
    }
    n[0] = n[0].wrapping_add(2);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_sub(TWO_P_TOP);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(TWO_P_TOP & carry);
    let _ = prop(n);
}

#[inline]
fn modsub(a: &[u64; NLIMBS], b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = a[i].wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(TWO_P_TOP & carry);
    let _ = prop(n);
}

#[inline]
fn modneg(b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = 0u64.wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(TWO_P_TOP & carry);
    let _ = prop(n);
}

#[inline]
fn modcpy(a: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    *c = *a;
}

#[inline]
fn modnsqr(a: &mut [u64; NLIMBS], n: u32) {
    for _ in 0..n {
        let t = *a;
        modsqr(&t, a);
    }
}

fn modinv(x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>, z: &mut [u64; NLIMBS]) {
    let mut s = [0u64; NLIMBS];
    let mut t = [0u64; NLIMBS];
    match h {
        None => modpro(x, &mut t),
        Some(hh) => modcpy(hh, &mut t),
    }
    modcpy(x, &mut s);
    modnsqr(&mut t, 2);
    modmul(&s, &t, z);
}

#[inline]
fn nres(m: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    modmul(m, &NRES_C, n);
}

#[inline]
fn redc(n: &[u64; NLIMBS], m: &mut [u64; NLIMBS]) {
    let mut c = [0u64; NLIMBS];
    c[0] = 1;
    modmul(n, &c, m);
    let _ = modfsb(m);
}

fn modis1(a: &[u64; NLIMBS]) -> i32 {
    let mut c = [0u64; NLIMBS];
    redc(a, &mut c);
    let mut d: u64 = 0;
    for i in 1..NLIMBS {
        d |= c[i];
    }
    let c0 = c[0];
    (1u64 & (d.wrapping_sub(1) >> LIMB_BITS) & ((c0 ^ 1).wrapping_sub(1) >> LIMB_BITS)) as i32
}

fn modis0(a: &[u64; NLIMBS]) -> i32 {
    let mut c = [0u64; NLIMBS];
    redc(a, &mut c);
    let mut d: u64 = 0;
    for i in 0..NLIMBS {
        d |= c[i];
    }
    (1u64 & (d.wrapping_sub(1) >> LIMB_BITS)) as i32
}

#[inline]
fn modint(x: i32, a: &mut [u64; NLIMBS]) {
    *a = [0u64; NLIMBS];
    a[0] = x as u64;
    let t = *a;
    nres(&t, a);
}

#[inline]
fn modmli(a: &[u64; NLIMBS], b: i32, c: &mut [u64; NLIMBS]) {
    let mut t = [0u64; NLIMBS];
    modint(b, &mut t);
    modmul(a, &t, c);
}

fn modqr(h: Option<&[u64; NLIMBS]>, x: &[u64; NLIMBS]) -> i32 {
    let mut r = [0u64; NLIMBS];
    match h {
        None => {
            modpro(x, &mut r);
            let t = r;
            modsqr(&t, &mut r);
        }
        Some(hh) => modsqr(hh, &mut r),
    }
    let t = r;
    modmul(&t, x, &mut r);
    modis1(&r) | modis0(x)
}

#[inline(never)]
fn modcsw(b: i32, g: &mut [u64; NLIMBS], f: &mut [u64; NLIMBS]) {
    let r: u64 = 0x3cc3_c33c_5aa5_a55a;
    let c0 = (1u64.wrapping_sub(b as u64)).wrapping_add(r);
    let c1 = (b as u64).wrapping_add(r);
    for i in 0..NLIMBS {
        let s = g[i];
        let t = f[i];
        let w = r.wrapping_mul(t.wrapping_add(s));
        f[i] = c0
            .wrapping_mul(t)
            .wrapping_add(c1.wrapping_mul(s))
            .wrapping_sub(w);
        g[i] = c0
            .wrapping_mul(s)
            .wrapping_add(c1.wrapping_mul(t))
            .wrapping_sub(w);
    }
}

fn modsqrt(x: &[u64; NLIMBS], h: Option<&[u64; NLIMBS]>, r: &mut [u64; NLIMBS]) {
    let mut s = [0u64; NLIMBS];
    let mut y = [0u64; NLIMBS];
    match h {
        None => modpro(x, &mut y),
        Some(hh) => modcpy(hh, &mut y),
    }
    modmul(&y, x, &mut s);
    modcpy(&s, r);
}

#[inline]
fn modshl(n: u32, a: &mut [u64; NLIMBS]) {
    a[NLIMBS - 1] = (a[NLIMBS - 1] << n) | (a[NLIMBS - 2] >> (LIMB_BITS - n));
    for i in (1..(NLIMBS - 1)).rev() {
        a[i] = ((a[i] << n) & MASK) | (a[i - 1] >> (LIMB_BITS - n));
    }
    a[0] = (a[0] << n) & MASK;
}

#[inline]
fn modshr(n: u32, a: &mut [u64; NLIMBS]) -> u64 {
    let r = a[0] & ((1u64 << n) - 1);
    for i in 0..(NLIMBS - 1) {
        a[i] = (a[i] >> n) | ((a[i + 1] << (LIMB_BITS - n)) & MASK);
    }
    a[NLIMBS - 1] >>= n;
    r
}

fn modcmp(a: &[u64; NLIMBS], b: &[u64; NLIMBS]) -> i32 {
    let mut c = [0u64; NLIMBS];
    let mut d = [0u64; NLIMBS];
    redc(a, &mut c);
    redc(b, &mut d);
    let mut eq: u64 = 1;
    for i in 0..NLIMBS {
        eq &= ((c[i] ^ d[i]).wrapping_sub(1) >> LIMB_BITS) & 1;
    }
    eq as i32
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

    fn fp_random(prng: &mut Prng) -> Fp {
        prng.fp()
    }

    const ITERS: usize = 500;

    #[test]
    fn equality_and_zero() {
        let mut prng = Prng(0xDEAD_BEEF_1234_5678);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = a + Fp::ONE;
            assert_eq!(a.is_equal_ct(&a), 0xFFFF_FFFF);
            assert_eq!(a.is_equal_ct(&b), 0);
            assert!(Fp::ZERO.is_zero());
            assert!(!Fp::ONE.is_zero());
        }
    }

    #[test]
    fn mul_small() {
        let mut prng = Prng(0x1111_2222_3333_4444);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let val = (prng.next() as u32) & 0x7FFF_FFFF;
            assert_eq!(a.mul_small(val), a * Fp::from_small(val as u64));
        }
    }

    #[test]
    fn half_and_div3() {
        let mut prng = Prng(0xAAAA_BBBB_CCCC_DDDD);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            assert_eq!((a + a).half(), a);
            assert_eq!((a + a + a).div3(), a);
        }
    }

    #[test]
    fn set_small() {
        assert_eq!(Fp::ONE + Fp::ONE, Fp::from_small(2));
        assert_eq!(Fp::ONE.dbl().dbl(), Fp::from_small(4));
    }

    #[test]
    fn addition_laws() {
        let mut prng = Prng(0x1);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            assert_eq!((a + b) + c, a + (b + c));
            assert_eq!(a + b, b + a);
            assert_eq!(a + Fp::ZERO, a);
            assert!((a + (-a)).is_zero());
        }
    }

    #[test]
    fn subtraction_laws() {
        let mut prng = Prng(0x2);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            assert_eq!((a - b) - c, a - (b + c));
            assert_eq!(a - b, -(b - a));
            assert!((a - a).is_zero());
        }
    }

    #[test]
    fn multiplication_laws() {
        let mut prng = Prng(0x3);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a * b, b * a);
            assert_eq!(a * Fp::ONE, a);
            assert!((a * Fp::ZERO).is_zero());
        }
    }

    #[test]
    fn squaring() {
        let mut prng = Prng(0x4);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            assert_eq!(a.square(), a * a);
        }
        assert!(Fp::ZERO.square().is_zero());
    }

    #[test]
    fn inversion() {
        let mut prng = Prng(0x5);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            assert_eq!(a * a.inv(), Fp::ONE);
        }
        assert!(Fp::ZERO.inv().is_zero());
    }

    #[test]
    fn sqrt_and_is_square() {
        let mut prng = Prng(0x6);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let c = a.square();
            assert!(c.is_square());
            let r = c.sqrt();
            assert!(r == a || r == -a);
        }
    }

    #[test]
    fn encode_decode_roundtrip() {
        let mut prng = Prng(0x7);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            assert_eq!(Fp::try_decode(&a.encode()), Some(a));
        }
    }

    #[test]
    fn decode_rejects_noncanonical() {
        let buf = [0xFFu8; FP_ENCODED_BYTES];
        assert_eq!(Fp::try_decode(&buf), None);
    }

    /// Golden vectors extracted from the C reference (libsqisign_gf_lvl1.a).
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    #[test]
    fn golden_c_vectors() {
        assert_hex(
            &Fp::ONE.encode(),
            "0100000000000000000000000000000000000000000000000000000000000000",
        );
        let c = Fp::from_small(5) * Fp::from_small(7);
        assert_hex(
            &c.encode(),
            "2300000000000000000000000000000000000000000000000000000000000000",
        );
        assert_hex(
            &c.inv().encode(),
            "9224499224499224499224499224499224499224499224499224499224499201",
        );
        assert_hex(
            &Fp::from_small(12345).square().sqrt().encode(),
            "c6cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
        );

        let mut src = [0u8; 64];
        for i in 0..64 {
            src[i] = (i as u8).wrapping_mul(7).wrapping_add(1);
        }
        assert_hex(
            &Fp::decode_reduce(&src).encode(),
            "a43cd712e8565f34a3ab4d895ecdd577b388f7ffa1ddb2212acc07dd4b54f602",
        );
        assert_hex(
            &Fp::decode_reduce(&src[..47]).encode(),
            "2c35d712e8565f34a3ab4d895ecdd57771787f868d949ba2a9b0b7bec5ccd303",
        );
    }

    #[test]
    fn precomp_fp_constants_match() {
        assert_eq!(Fp::ONE, crate::precomp::FP2_CONSTANTS[1].re);
        assert_eq!(-Fp::ONE, crate::precomp::FP2_CONSTANTS[3].re);
    }

    #[test]
    fn select_and_cswap() {
        let mut prng = Prng(0x8);
        let a = fp_random(&mut prng);
        let b = fp_random(&mut prng);
        assert_eq!(Fp::select(&a, &b, 0), a);
        assert_eq!(Fp::select(&a, &b, 0xFFFF_FFFF), b);

        let mut x = a;
        let mut y = b;
        Fp::cswap(&mut x, &mut y, 0);
        assert_eq!((x, y), (a, b));
        Fp::cswap(&mut x, &mut y, 0xFFFF_FFFF);
        assert_eq!((x, y), (b, a));
    }
}
