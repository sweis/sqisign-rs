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
#[derive(Clone, Copy, Eq, PartialEq)]
#[repr(C)]
pub struct Fp(pub [u64; NWORDS_FIELD]);

impl Default for Fp {
    fn default() -> Self {
        Fp([0; NWORDS_FIELD])
    }
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buf = [0u8; FP_ENCODED_BYTES];
        fp_encode(&mut buf, self);
        write!(f, "Fp(0x")?;
        for b in buf.iter().rev() {
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
fn modzer(a: &mut [u64; NLIMBS]) {
    *a = [0u64; NLIMBS];
}

#[inline]
fn modone(a: &mut [u64; NLIMBS]) {
    *a = [0u64; NLIMBS];
    a[0] = 1;
    let t = *a;
    nres(&t, a);
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
// Public fp_* API (port of the wrapper section + lvlx/fp.c)
// ===========================================================================

#[inline]
pub fn fp_set_small(x: &mut Fp, val: Digit) {
    modint(val as i32, &mut x.0);
}

#[inline]
pub fn fp_mul_small(x: &mut Fp, a: &Fp, val: u32) {
    modmli(&a.0, val as i32, &mut x.0);
}

#[inline]
pub fn fp_set_zero(x: &mut Fp) {
    modzer(&mut x.0);
}

#[inline]
pub fn fp_set_one(x: &mut Fp) {
    modone(&mut x.0);
}

#[inline]
pub fn fp_is_equal(a: &Fp, b: &Fp) -> u32 {
    (modcmp(&a.0, &b.0) as u32).wrapping_neg()
}

#[inline]
pub fn fp_is_zero(a: &Fp) -> u32 {
    (modis0(&a.0) as u32).wrapping_neg()
}

#[inline]
pub fn fp_copy(out: &mut Fp, a: &Fp) {
    out.0 = a.0;
}

#[inline]
pub fn fp_cswap(a: &mut Fp, b: &mut Fp, ctl: u32) {
    modcsw((ctl & 1) as i32, &mut a.0, &mut b.0);
}

#[inline]
pub fn fp_add(out: &mut Fp, a: &Fp, b: &Fp) {
    modadd(&a.0, &b.0, &mut out.0);
}

#[inline]
pub fn fp_sub(out: &mut Fp, a: &Fp, b: &Fp) {
    modsub(&a.0, &b.0, &mut out.0);
}

#[inline]
pub fn fp_neg(out: &mut Fp, a: &Fp) {
    modneg(&a.0, &mut out.0);
}

#[inline]
pub fn fp_sqr(out: &mut Fp, a: &Fp) {
    modsqr(&a.0, &mut out.0);
}

#[inline]
pub fn fp_mul(out: &mut Fp, a: &Fp, b: &Fp) {
    modmul(&a.0, &b.0, &mut out.0);
}

// ---------------------------------------------------------------------------
// In-place variants. `modmul`/`modsqr` write `c[k]` only after they've finished
// reading the input limbs needed for column k, so they are alias-safe; the
// stack copy here exists only to satisfy the borrow checker and SROA's away.
// ---------------------------------------------------------------------------

#[inline]
pub fn fp_mul_ip(r: &mut Fp, b: &Fp) {
    let a = r.0;
    modmul(&a, &b.0, &mut r.0);
}

#[inline]
pub fn fp_sqr_ip(r: &mut Fp) {
    let a = r.0;
    modsqr(&a, &mut r.0);
}

#[inline]
pub fn fp_add_ip(r: &mut Fp, b: &Fp) {
    let a = r.0;
    modadd(&a, &b.0, &mut r.0);
}

#[inline]
pub fn fp_sub_ip(r: &mut Fp, b: &Fp) {
    let a = r.0;
    modsub(&a, &b.0, &mut r.0);
}

#[inline]
pub fn fp_neg_ip(r: &mut Fp) {
    let a = r.0;
    modneg(&a, &mut r.0);
}

#[inline]
pub fn fp_dbl_ip(r: &mut Fp) {
    let a = r.0;
    modadd(&a, &a, &mut r.0);
}

#[inline]
pub fn fp_inv(x: &mut Fp) {
    let a = x.0;
    modinv(&a, None, &mut x.0);
}

#[inline]
pub fn fp_is_square(a: &Fp) -> u32 {
    (modqr(None, &a.0) as u32).wrapping_neg()
}

#[inline]
pub fn fp_sqrt(a: &mut Fp) {
    let x = a.0;
    modsqrt(&x, None, &mut a.0);
}

#[inline]
pub fn fp_half(out: &mut Fp, a: &Fp) {
    modmul(&TWO_INV.0, &a.0, &mut out.0);
}

#[inline]
pub fn fp_exp3div4(out: &mut Fp, a: &Fp) {
    modpro(&a.0, &mut out.0);
}

#[inline]
pub fn fp_div3(out: &mut Fp, a: &Fp) {
    modmul(&THREE_INV.0, &a.0, &mut out.0);
}

#[inline]
pub fn fp_select(d: &mut Fp, a0: &Fp, a1: &Fp, ctl: u32) {
    let cw = (ctl as i32) as i64 as u64;
    for i in 0..NWORDS_FIELD {
        d.0[i] = a0.0[i] ^ (cw & (a0.0[i] ^ a1.0[i]));
    }
}

pub fn fp_encode(dst: &mut [u8], a: &Fp) {
    debug_assert!(dst.len() >= FP_ENCODED_BYTES);
    let mut c = [0u64; NLIMBS];
    redc(&a.0, &mut c);
    for i in 0..NBYTES {
        dst[i] = (c[0] & 0xff) as u8;
        let _ = modshr(8, &mut c);
    }
}

pub fn fp_decode(d: &mut Fp, src: &[u8]) -> u32 {
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

/// Decode an arbitrary-length little-endian byte string and reduce mod p.
pub fn fp_decode_reduce(d: &mut Fp, src: &[u8]) {
    fp_set_zero(d);
    let mut len = src.len();
    if len == 0 {
        return;
    }
    let mut tmp = [0u8; NBYTES];

    let rem = len % NBYTES;
    if rem != 0 {
        let k = len - rem;
        tmp[..rem].copy_from_slice(&src[k..len]);
        tmp[rem..].fill(0);
        fp_decode(d, &tmp);
        len = k;
    }

    while len > 0 {
        let prev = *d;
        fp_mul(d, &prev, &R2);
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
        let mut a = Fp::default();
        fp_decode(&mut a, &tmp);
        let prev = *d;
        fp_add(d, &prev, &a);
    }
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
    fn one_constant_matches_set_one() {
        let mut a = Fp::default();
        fp_set_one(&mut a);
        assert_eq!(fp_is_equal(&a, &ONE), 0xFFFF_FFFF);
    }

    #[test]
    fn equality_and_zero() {
        let mut prng = Prng(0xDEAD_BEEF_1234_5678);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut b = Fp::default();
            fp_add(&mut b, &a, &ONE);
            let mut c = Fp::default();
            fp_set_zero(&mut c);

            assert_ne!(fp_is_equal(&a, &a), 0);
            assert_eq!(fp_is_equal(&a, &b), 0);
            assert_ne!(fp_is_equal(&c, &ZERO), 0);
            assert_ne!(fp_is_zero(&ZERO), 0);
            assert_eq!(fp_is_zero(&ONE), 0);
        }
    }

    #[test]
    fn mul_small() {
        let mut prng = Prng(0x1111_2222_3333_4444);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let val = (prng.next() as u32) & 0x7FFF_FFFF;
            let mut b = Fp::default();
            fp_mul_small(&mut b, &a, val);
            let mut c = Fp::default();
            fp_set_small(&mut c, val as u64);
            let mut d = Fp::default();
            fp_mul(&mut d, &a, &c);
            assert_ne!(fp_is_equal(&b, &d), 0);
        }
    }

    #[test]
    fn half_and_div3() {
        let mut prng = Prng(0xAAAA_BBBB_CCCC_DDDD);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut b = Fp::default();
            fp_add(&mut b, &a, &a);
            let mut c = Fp::default();
            fp_half(&mut c, &b);
            assert_ne!(fp_is_equal(&a, &c), 0);

            let mut b = Fp::default();
            fp_add(&mut b, &a, &a);
            let b2 = b;
            fp_add(&mut b, &b2, &a);
            let mut c = Fp::default();
            fp_div3(&mut c, &b);
            assert_ne!(fp_is_equal(&a, &c), 0);
        }
    }

    #[test]
    fn set_small() {
        let mut a = Fp::default();
        fp_set_one(&mut a);
        let mut b = Fp::default();
        fp_add(&mut b, &a, &a);
        let mut c = Fp::default();
        fp_set_small(&mut c, 2);
        assert_ne!(fp_is_equal(&b, &c), 0);

        let b2 = b;
        fp_add(&mut b, &b2, &b2);
        fp_set_small(&mut c, 4);
        assert_ne!(fp_is_equal(&b, &c), 0);
    }

    #[test]
    fn addition_laws() {
        let mut prng = Prng(0x1);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            let mut d = Fp::default();
            let mut e = Fp::default();
            let mut f = Fp::default();

            fp_add(&mut d, &a, &b);
            fp_add(&mut e, &d, &c);
            fp_add(&mut d, &b, &c);
            fp_add(&mut f, &d, &a);
            assert_ne!(fp_is_equal(&e, &f), 0);

            fp_add(&mut d, &a, &b);
            fp_add(&mut e, &b, &a);
            assert_ne!(fp_is_equal(&d, &e), 0);

            let z = ZERO;
            fp_add(&mut d, &a, &z);
            assert_ne!(fp_is_equal(&a, &d), 0);

            fp_neg(&mut d, &a);
            fp_add(&mut e, &a, &d);
            assert_ne!(fp_is_zero(&e), 0);
        }
    }

    #[test]
    fn subtraction_laws() {
        let mut prng = Prng(0x2);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            let mut d = Fp::default();
            let mut e = Fp::default();
            let mut f = Fp::default();

            fp_sub(&mut d, &a, &b);
            fp_sub(&mut e, &d, &c);
            fp_add(&mut d, &b, &c);
            fp_sub(&mut f, &a, &d);
            assert_ne!(fp_is_equal(&e, &f), 0);

            fp_sub(&mut d, &a, &b);
            fp_sub(&mut e, &b, &a);
            let e2 = e;
            fp_neg(&mut e, &e2);
            assert_ne!(fp_is_equal(&d, &e), 0);

            fp_sub(&mut e, &a, &a);
            assert_ne!(fp_is_zero(&e), 0);
        }
    }

    #[test]
    fn multiplication_laws() {
        let mut prng = Prng(0x3);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let b = fp_random(&mut prng);
            let c = fp_random(&mut prng);
            let mut d = Fp::default();
            let mut e = Fp::default();
            let mut f = Fp::default();

            fp_mul(&mut d, &a, &b);
            fp_mul(&mut e, &d, &c);
            fp_mul(&mut d, &b, &c);
            fp_mul(&mut f, &d, &a);
            assert_ne!(fp_is_equal(&e, &f), 0);

            fp_add(&mut d, &b, &c);
            fp_mul(&mut e, &a, &d);
            fp_mul(&mut d, &a, &b);
            fp_mul(&mut f, &a, &c);
            let f2 = f;
            fp_add(&mut f, &d, &f2);
            assert_ne!(fp_is_equal(&e, &f), 0);

            fp_mul(&mut d, &a, &b);
            fp_mul(&mut e, &b, &a);
            assert_ne!(fp_is_equal(&d, &e), 0);

            fp_mul(&mut d, &a, &ONE);
            assert_ne!(fp_is_equal(&a, &d), 0);

            fp_mul(&mut d, &a, &ZERO);
            assert_ne!(fp_is_zero(&d), 0);
        }
    }

    #[test]
    fn squaring() {
        let mut prng = Prng(0x4);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut b = Fp::default();
            let mut c = Fp::default();
            fp_sqr(&mut b, &a);
            fp_mul(&mut c, &a, &a);
            assert_ne!(fp_is_equal(&b, &c), 0);
        }
        let mut d = Fp::default();
        fp_sqr(&mut d, &ZERO);
        assert_ne!(fp_is_zero(&d), 0);
    }

    #[test]
    fn inversion() {
        let mut prng = Prng(0x5);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut b = a;
            fp_inv(&mut b);
            let mut c = Fp::default();
            fp_mul(&mut c, &a, &b);
            assert_ne!(fp_is_equal(&c, &ONE), 0);
        }
        let mut z = ZERO;
        fp_inv(&mut z);
        assert_ne!(fp_is_zero(&z), 0);
    }

    #[test]
    fn sqrt_and_is_square() {
        let mut prng = Prng(0x6);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut c = Fp::default();
            fp_sqr(&mut c, &a);
            assert_ne!(fp_is_square(&c), 0);
            fp_sqrt(&mut c);
            let mut d = Fp::default();
            fp_neg(&mut d, &c);
            assert!(fp_is_equal(&a, &c) != 0 || fp_is_equal(&a, &d) != 0);
        }
    }

    #[test]
    fn encode_decode_roundtrip() {
        let mut prng = Prng(0x7);
        for _ in 0..ITERS {
            let a = fp_random(&mut prng);
            let mut buf = [0u8; FP_ENCODED_BYTES];
            fp_encode(&mut buf, &a);
            let mut b = Fp::default();
            let ok = fp_decode(&mut b, &buf);
            assert_eq!(ok, 0xFFFF_FFFF);
            assert_ne!(fp_is_equal(&a, &b), 0);
        }
    }

    #[test]
    fn decode_rejects_noncanonical() {
        let buf = [0xFFu8; FP_ENCODED_BYTES];
        let mut d = Fp::default();
        let ok = fp_decode(&mut d, &buf);
        assert_eq!(ok, 0);
        assert_ne!(fp_is_zero(&d), 0);
    }

    /// Golden vectors extracted from the C reference (libsqisign_gf_lvl1.a).
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    #[test]
    fn golden_c_vectors() {
        let mut enc = [0u8; FP_ENCODED_BYTES];

        let mut a = Fp::default();
        fp_set_one(&mut a);
        fp_encode(&mut enc, &a);
        assert_hex(
            &enc,
            "0100000000000000000000000000000000000000000000000000000000000000",
        );

        let mut b = Fp::default();
        let mut c = Fp::default();
        fp_set_small(&mut a, 5);
        fp_set_small(&mut b, 7);
        fp_mul(&mut c, &a, &b);
        fp_encode(&mut enc, &c);
        assert_hex(
            &enc,
            "2300000000000000000000000000000000000000000000000000000000000000",
        );

        fp_inv(&mut c);
        fp_encode(&mut enc, &c);
        assert_hex(
            &enc,
            "9224499224499224499224499224499224499224499224499224499224499201",
        );

        fp_set_small(&mut a, 12345);
        fp_sqr(&mut b, &a);
        fp_sqrt(&mut b);
        fp_encode(&mut enc, &b);
        assert_hex(
            &enc,
            "c6cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
        );

        let mut src = [0u8; 64];
        for i in 0..64 {
            src[i] = (i as u8).wrapping_mul(7).wrapping_add(1);
        }
        fp_decode_reduce(&mut a, &src);
        fp_encode(&mut enc, &a);
        assert_hex(
            &enc,
            "a43cd712e8565f34a3ab4d895ecdd577b388f7ffa1ddb2212acc07dd4b54f602",
        );

        fp_decode_reduce(&mut a, &src[..47]);
        fp_encode(&mut enc, &a);
        assert_hex(
            &enc,
            "2c35d712e8565f34a3ab4d895ecdd57771787f868d949ba2a9b0b7bec5ccd303",
        );
    }

    #[test]
    fn precomp_fp_constants_match() {
        let mut neg1 = Fp::default();
        fp_set_one(&mut neg1);
        let one = neg1;
        fp_neg(&mut neg1, &one);
        assert_eq!(
            fp_is_equal(&one, &crate::precomp::FP2_CONSTANTS[1].re),
            0xFFFF_FFFF
        );
        assert_eq!(
            fp_is_equal(&neg1, &crate::precomp::FP2_CONSTANTS[3].re),
            0xFFFF_FFFF
        );
    }

    #[test]
    fn select_and_cswap() {
        let mut prng = Prng(0x8);
        let a = fp_random(&mut prng);
        let b = fp_random(&mut prng);
        let mut d = Fp::default();
        fp_select(&mut d, &a, &b, 0);
        assert_ne!(fp_is_equal(&d, &a), 0);
        fp_select(&mut d, &a, &b, 0xFFFF_FFFF);
        assert_ne!(fp_is_equal(&d, &b), 0);

        let mut x = a;
        let mut y = b;
        fp_cswap(&mut x, &mut y, 0);
        assert_ne!(fp_is_equal(&x, &a), 0);
        assert_ne!(fp_is_equal(&y, &b), 0);
        fp_cswap(&mut x, &mut y, 0xFFFF_FFFF);
        assert_ne!(fp_is_equal(&x, &b), 0);
        assert_ne!(fp_is_equal(&y, &a), 0);
    }
}
