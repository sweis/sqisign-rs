// SPDX-License-Identifier: Apache-2.0
//! GF(p²) arithmetic, modulo X² + 1. Port of `lvlx/fp2.c`.

use super::fp::*;
use core::fmt;

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

#[inline]
pub fn fp2_set_small(x: &mut Fp2, val: Digit) {
    fp_set_small(&mut x.re, val);
    fp_set_zero(&mut x.im);
}

#[inline]
pub fn fp2_mul_small(x: &mut Fp2, y: &Fp2, n: u32) {
    fp_mul_small(&mut x.re, &y.re, n);
    fp_mul_small(&mut x.im, &y.im, n);
}

#[inline]
pub fn fp2_set_one(x: &mut Fp2) {
    fp_set_one(&mut x.re);
    fp_set_zero(&mut x.im);
}

#[inline]
pub fn fp2_set_zero(x: &mut Fp2) {
    fp_set_zero(&mut x.re);
    fp_set_zero(&mut x.im);
}

#[inline]
pub fn fp2_is_zero(a: &Fp2) -> u32 {
    fp_is_zero(&a.re) & fp_is_zero(&a.im)
}

#[inline]
pub fn fp2_is_equal(a: &Fp2, b: &Fp2) -> u32 {
    fp_is_equal(&a.re, &b.re) & fp_is_equal(&a.im, &b.im)
}

#[inline]
pub fn fp2_is_one(a: &Fp2) -> u32 {
    fp_is_equal(&a.re, &ONE) & fp_is_zero(&a.im)
}

#[inline]
pub fn fp2_copy(x: &mut Fp2, y: &Fp2) {
    *x = *y;
}

#[inline]
pub fn fp2_add(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    fp_add(&mut x.re, &y.re, &z.re);
    fp_add(&mut x.im, &y.im, &z.im);
}

#[inline]
pub fn fp2_add_one(x: &mut Fp2, y: &Fp2) {
    fp_add(&mut x.re, &y.re, &ONE);
    fp_copy(&mut x.im, &y.im);
}

#[inline]
pub fn fp2_sub(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    fp_sub(&mut x.re, &y.re, &z.re);
    fp_sub(&mut x.im, &y.im, &z.im);
}

#[inline]
pub fn fp2_neg(x: &mut Fp2, y: &Fp2) {
    fp_neg(&mut x.re, &y.re);
    fp_neg(&mut x.im, &y.im);
}

/// (y.re + y.im·i)(z.re + z.im·i), via Karatsuba.
pub fn fp2_mul(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    let mut t0 = Fp::default();
    let mut t1 = Fp::default();

    fp_add(&mut t0, &y.re, &y.im);
    fp_add(&mut t1, &z.re, &z.im);
    let s0 = t0;
    fp_mul(&mut t0, &s0, &t1);
    fp_mul(&mut t1, &y.im, &z.im);
    fp_mul(&mut x.re, &y.re, &z.re);
    fp_sub(&mut x.im, &t0, &t1);
    let im = x.im;
    fp_sub(&mut x.im, &im, &x.re);
    let re = x.re;
    fp_sub(&mut x.re, &re, &t1);
}

pub fn fp2_sqr(x: &mut Fp2, y: &Fp2) {
    let mut sum = Fp::default();
    let mut diff = Fp::default();

    fp_add(&mut sum, &y.re, &y.im);
    fp_sub(&mut diff, &y.re, &y.im);
    fp_mul(&mut x.im, &y.re, &y.im);
    let im = x.im;
    fp_add(&mut x.im, &im, &im);
    fp_mul(&mut x.re, &sum, &diff);
}

pub fn fp2_inv(x: &mut Fp2) {
    let mut t0 = Fp::default();
    let mut t1 = Fp::default();

    fp_sqr(&mut t0, &x.re);
    fp_sqr(&mut t1, &x.im);
    let s = t0;
    fp_add(&mut t0, &s, &t1);
    fp_inv(&mut t0);
    let re = x.re;
    fp_mul(&mut x.re, &re, &t0);
    let im = x.im;
    fp_mul(&mut x.im, &im, &t0);
    let im = x.im;
    fp_neg(&mut x.im, &im);
}

pub fn fp2_is_square(x: &Fp2) -> u32 {
    let mut t0 = Fp::default();
    let mut t1 = Fp::default();
    fp_sqr(&mut t0, &x.re);
    fp_sqr(&mut t1, &x.im);
    let s = t0;
    fp_add(&mut t0, &s, &t1);
    fp_is_square(&t0)
}

/// Square root in GF(p²) following Aardal et al. (ePrint 2024/1563),
/// canonicalized so the result has even real part (or, if re==0, even imag).
pub fn fp2_sqrt(a: &mut Fp2) {
    let mut x0 = Fp::default();
    let mut x1 = Fp::default();
    let mut t0 = Fp::default();
    let mut t1 = Fp::default();

    fp_sqr(&mut x0, &a.re);
    fp_sqr(&mut x1, &a.im);
    let s = x0;
    fp_add(&mut x0, &s, &x1);
    fp_sqrt(&mut x0);
    let s = x0;
    fp_select(&mut x0, &s, &a.re, fp_is_zero(&a.im));
    let s = x0;
    fp_add(&mut x0, &s, &a.re);
    fp_add(&mut t0, &x0, &x0);

    fp_exp3div4(&mut x1, &t0);

    let s = x0;
    fp_mul(&mut x0, &s, &x1);
    let s = x1;
    fp_mul(&mut x1, &s, &a.im);
    fp_add(&mut t1, &x0, &x0);
    let s = t1;
    fp_sqr(&mut t1, &s);
    let s = t0;
    fp_sub(&mut t0, &s, &t1);
    let f = fp_is_zero(&t0);
    fp_neg(&mut t1, &x0);
    fp_copy(&mut t0, &x1);
    let s0 = t0;
    fp_select(&mut t0, &s0, &x0, f);
    let s1 = t1;
    fp_select(&mut t1, &s1, &x1, f);

    let t0_is_zero = fp_is_zero(&t0);

    let mut tmp = [0u8; FP_ENCODED_BYTES];
    fp_encode(&mut tmp, &t0);
    let t0_is_odd = ((tmp[0] & 1) as u32).wrapping_neg();
    fp_encode(&mut tmp, &t1);
    let t1_is_odd = ((tmp[0] & 1) as u32).wrapping_neg();

    let negate = t0_is_odd | (t0_is_zero & t1_is_odd);
    fp_neg(&mut x0, &t0);
    fp_select(&mut a.re, &t0, &x0, negate);
    fp_neg(&mut x0, &t1);
    fp_select(&mut a.im, &t1, &x0, negate);
}

/// Compute sqrt(a) into a; returns 0xFFFFFFFF iff the result squares back to the input.
pub fn fp2_sqrt_verify(a: &mut Fp2) -> u32 {
    let t0 = *a;
    fp2_sqrt(a);
    let mut t1 = Fp2::default();
    fp2_sqr(&mut t1, a);
    fp2_is_equal(&t0, &t1)
}

#[inline]
pub fn fp2_half(x: &mut Fp2, y: &Fp2) {
    fp_half(&mut x.re, &y.re);
    fp_half(&mut x.im, &y.im);
}

/// In-place batched inversion (Montgomery's trick).
pub fn fp2_batched_inv(x: &mut [Fp2]) {
    let len = x.len();
    if len == 0 {
        return;
    }
    let mut t1 = vec![Fp2::default(); len];
    let mut t2 = vec![Fp2::default(); len];

    t1[0] = x[0];
    for i in 1..len {
        let prev = t1[i - 1];
        fp2_mul(&mut t1[i], &prev, &x[i]);
    }
    let mut inverse = t1[len - 1];
    fp2_inv(&mut inverse);

    t2[0] = inverse;
    for i in 1..len {
        let prev = t2[i - 1];
        fp2_mul(&mut t2[i], &prev, &x[len - i]);
    }
    x[0] = t2[len - 1];
    for i in 1..len {
        let a = t1[i - 1];
        let b = t2[len - i - 1];
        fp2_mul(&mut x[i], &a, &b);
    }
}

/// Variable-time square-and-multiply: out = x^exp, where exp is little-endian limbs.
pub fn fp2_pow_vartime(out: &mut Fp2, x: &Fp2, exp: &[Digit]) {
    let mut acc = *x;
    fp2_set_one(out);
    for &word in exp.iter() {
        for i in 0..RADIX {
            if (word >> i) & 1 == 1 {
                let o = *out;
                fp2_mul(out, &o, &acc);
            }
            let a = acc;
            fp2_sqr(&mut acc, &a);
        }
    }
}

pub fn fp2_encode(dst: &mut [u8], a: &Fp2) {
    fp_encode(&mut dst[..FP_ENCODED_BYTES], &a.re);
    fp_encode(&mut dst[FP_ENCODED_BYTES..2 * FP_ENCODED_BYTES], &a.im);
}

pub fn fp2_decode(d: &mut Fp2, src: &[u8]) -> u32 {
    let re = fp_decode(&mut d.re, &src[..FP_ENCODED_BYTES]);
    let im = fp_decode(&mut d.im, &src[FP_ENCODED_BYTES..2 * FP_ENCODED_BYTES]);
    re & im
}

#[inline]
pub fn fp2_select(d: &mut Fp2, a0: &Fp2, a1: &Fp2, ctl: u32) {
    fp_select(&mut d.re, &a0.re, &a1.re, ctl);
    fp_select(&mut d.im, &a0.im, &a1.im, ctl);
}

#[inline]
pub fn fp2_cswap(a: &mut Fp2, b: &mut Fp2, ctl: u32) {
    fp_cswap(&mut a.re, &mut b.re, ctl);
    fp_cswap(&mut a.im, &mut b.im, ctl);
}

/// Frobenius endomorphism: out = in^p = conj(in) in GF(p²).
#[inline]
pub fn fp2_frob(out: &mut Fp2, input: &Fp2) {
    fp_copy(&mut out.re, &input.re);
    fp_neg(&mut out.im, &input.im);
}

pub fn fp2_print(name: &str, a: &Fp2) {
    let mut buf = [0u8; FP_ENCODED_BYTES];
    fp_encode(&mut buf, &a.re);
    print!("{}0x", name);
    for b in buf.iter().rev() {
        print!("{:02x}", b);
    }
    print!(" + i*0x");
    fp_encode(&mut buf, &a.im);
    for b in buf.iter().rev() {
        print!("{:02x}", b);
    }
    println!();
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    struct Prng(u64);
    impl Prng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn fill(&mut self, buf: &mut [u8]) {
            for chunk in buf.chunks_mut(8) {
                let x = self.next().to_le_bytes();
                chunk.copy_from_slice(&x[..chunk.len()]);
            }
        }
    }

    fn fp2_random(prng: &mut Prng) -> Fp2 {
        let mut buf = [0u8; FP_ENCODED_BYTES];
        let mut a = Fp2::default();
        prng.fill(&mut buf);
        fp_decode_reduce(&mut a.re, &buf);
        prng.fill(&mut buf);
        fp_decode_reduce(&mut a.im, &buf);
        a
    }

    const ITERS: usize = 300;

    #[test]
    fn addition_laws() {
        let mut prng = Prng(0x10);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            let mut d = Fp2::default();
            let mut e = Fp2::default();
            let mut f = Fp2::default();

            fp2_add(&mut d, &a, &b);
            fp2_add(&mut e, &d, &c);
            fp2_add(&mut d, &b, &c);
            fp2_add(&mut f, &d, &a);
            assert_ne!(fp2_is_equal(&e, &f), 0);

            fp2_add(&mut d, &a, &b);
            fp2_add(&mut e, &b, &a);
            assert_ne!(fp2_is_equal(&d, &e), 0);

            fp2_neg(&mut d, &a);
            fp2_add(&mut e, &a, &d);
            assert_ne!(fp2_is_zero(&e), 0);

            let mut one = Fp2::default();
            fp2_set_one(&mut one);
            fp2_add(&mut e, &a, &one);
            fp2_add_one(&mut f, &a);
            assert_ne!(fp2_is_equal(&e, &f), 0);
        }
    }

    #[test]
    fn subtraction_laws() {
        let mut prng = Prng(0x11);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            let mut d = Fp2::default();
            let mut e = Fp2::default();
            let mut f = Fp2::default();

            fp2_sub(&mut d, &a, &b);
            fp2_sub(&mut e, &d, &c);
            fp2_add(&mut d, &b, &c);
            fp2_sub(&mut f, &a, &d);
            assert_ne!(fp2_is_equal(&e, &f), 0);

            fp2_sub(&mut d, &a, &b);
            fp2_sub(&mut e, &b, &a);
            let e2 = e;
            fp2_neg(&mut e, &e2);
            assert_ne!(fp2_is_equal(&d, &e), 0);

            fp2_sub(&mut e, &a, &a);
            assert_ne!(fp2_is_zero(&e), 0);
        }
    }

    #[test]
    fn multiplication_laws() {
        let mut prng = Prng(0x12);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let b = fp2_random(&mut prng);
            let c = fp2_random(&mut prng);
            let mut d = Fp2::default();
            let mut e = Fp2::default();
            let mut f = Fp2::default();

            fp2_mul(&mut d, &a, &b);
            fp2_mul(&mut e, &d, &c);
            fp2_mul(&mut d, &b, &c);
            fp2_mul(&mut f, &d, &a);
            assert_ne!(fp2_is_equal(&e, &f), 0);

            fp2_add(&mut d, &b, &c);
            fp2_mul(&mut e, &a, &d);
            fp2_mul(&mut d, &a, &b);
            fp2_mul(&mut f, &a, &c);
            let f2 = f;
            fp2_add(&mut f, &d, &f2);
            assert_ne!(fp2_is_equal(&e, &f), 0);

            fp2_mul(&mut d, &a, &b);
            fp2_mul(&mut e, &b, &a);
            assert_ne!(fp2_is_equal(&d, &e), 0);

            let mut one = Fp2::default();
            fp2_set_one(&mut one);
            fp2_mul(&mut d, &a, &one);
            assert_ne!(fp2_is_equal(&a, &d), 0);

            let mut zero = Fp2::default();
            fp2_set_zero(&mut zero);
            fp2_mul(&mut d, &a, &zero);
            assert_ne!(fp2_is_zero(&d), 0);
        }
    }

    #[test]
    fn mul_matches_schoolbook() {
        // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        let mut prng = Prng(0x13);
        for _ in 0..ITERS {
            let y = fp2_random(&mut prng);
            let z = fp2_random(&mut prng);
            let mut x = Fp2::default();
            fp2_mul(&mut x, &y, &z);

            let mut ac = Fp::default();
            let mut bd = Fp::default();
            let mut ad = Fp::default();
            let mut bc = Fp::default();
            fp_mul(&mut ac, &y.re, &z.re);
            fp_mul(&mut bd, &y.im, &z.im);
            fp_mul(&mut ad, &y.re, &z.im);
            fp_mul(&mut bc, &y.im, &z.re);
            let mut re = Fp::default();
            let mut im = Fp::default();
            fp_sub(&mut re, &ac, &bd);
            fp_add(&mut im, &ad, &bc);
            assert_ne!(fp_is_equal(&x.re, &re), 0);
            assert_ne!(fp_is_equal(&x.im, &im), 0);
        }
    }

    #[test]
    fn mul_small() {
        let mut prng = Prng(0x14);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let val = (prng.next() as u32) & 0x7FFF_FFFF;
            let mut b = Fp2::default();
            fp2_mul_small(&mut b, &a, val);
            let mut c = Fp2::default();
            fp2_set_small(&mut c, val as Digit);
            let mut d = Fp2::default();
            fp2_mul(&mut d, &a, &c);
            assert_ne!(fp2_is_equal(&b, &d), 0);
        }
    }

    #[test]
    fn squaring() {
        let mut prng = Prng(0x15);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let mut b = Fp2::default();
            let mut c = Fp2::default();
            fp2_sqr(&mut b, &a);
            fp2_mul(&mut c, &a, &a);
            assert_ne!(fp2_is_equal(&b, &c), 0);
        }
    }

    #[test]
    fn inversion() {
        let mut prng = Prng(0x16);
        let mut one = Fp2::default();
        fp2_set_one(&mut one);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let mut b = a;
            fp2_inv(&mut b);
            let mut c = Fp2::default();
            fp2_mul(&mut c, &a, &b);
            assert_ne!(fp2_is_equal(&c, &one), 0);
        }
        let mut z = Fp2::default();
        fp2_inv(&mut z);
        assert_ne!(fp2_is_zero(&z), 0);
    }

    #[test]
    fn sqrt_and_is_square() {
        let mut prng = Prng(0x17);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let mut c = Fp2::default();
            fp2_sqr(&mut c, &a);
            assert_ne!(fp2_is_square(&c), 0);

            let mut b = c;
            assert_ne!(fp2_sqrt_verify(&mut b), 0);

            fp2_sqrt(&mut c);
            let mut d = Fp2::default();
            fp2_neg(&mut d, &c);
            assert!(fp2_is_equal(&a, &c) != 0 || fp2_is_equal(&a, &d) != 0);
        }
    }

    #[test]
    fn batched_inv() {
        let mut prng = Prng(0x18);
        let mut one = Fp2::default();
        fp2_set_one(&mut one);
        let mut xs: Vec<Fp2> = (0..8).map(|_| fp2_random(&mut prng)).collect();
        let orig = xs.clone();
        fp2_batched_inv(&mut xs);
        for (x, o) in xs.iter().zip(orig.iter()) {
            let mut p = Fp2::default();
            fp2_mul(&mut p, x, o);
            assert_ne!(fp2_is_equal(&p, &one), 0);
        }
    }

    fn assert_hex(buf: &[u8], expected: &str) {
        let mut got = String::new();
        for b in buf {
            got.push_str(&format!("{:02x}", b));
        }
        assert_eq!(got, expected);
    }

    /// Golden vectors extracted from the C reference implementation
    /// (libsqisign_gf_lvl1.a, RADIX_64, ref build).
    #[test]
    fn golden_c_vectors() {
        let mut enc = [0u8; FP2_ENCODED_BYTES];
        let mut x = Fp2::default();
        fp_set_small(&mut x.re, 3);
        fp_set_small(&mut x.im, 4);
        let mut y = Fp2::default();
        fp2_sqr(&mut y, &x);
        fp2_encode(&mut enc, &y);
        assert_hex(
            &enc,
            "f8ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04\
             1800000000000000000000000000000000000000000000000000000000000000",
        );

        fp2_sqrt(&mut y);
        fp2_encode(&mut enc, &y);
        assert_hex(
            &enc,
            "fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04\
             fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04",
        );

        fp2_inv(&mut y);
        fp2_encode(&mut enc, &y);
        assert_hex(
            &enc,
            "6666666666666666666666666666666666666666666666666666666666666601\
             cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc04",
        );
    }

    #[test]
    fn encode_decode_roundtrip() {
        let mut prng = Prng(0x19);
        for _ in 0..ITERS {
            let a = fp2_random(&mut prng);
            let mut buf = [0u8; FP2_ENCODED_BYTES];
            fp2_encode(&mut buf, &a);
            let mut b = Fp2::default();
            assert_eq!(fp2_decode(&mut b, &buf), 0xFFFF_FFFF);
            assert_ne!(fp2_is_equal(&a, &b), 0);
        }
    }
}
