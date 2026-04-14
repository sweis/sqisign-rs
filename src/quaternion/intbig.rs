// SPDX-License-Identifier: Apache-2.0
//! Big-integer wrapper over `rug::Integer`, ported from `intbig.c`.
//!
//! All `ibz_*` functions mirror the C signatures (out-params first) so the
//! rest of the quaternion module can be ported mechanically. This is
//! intentionally un-idiomatic; it will be refactored once KAT passes.

use crate::common::ctrdrbg::randombytes;
use crate::mp::Digit;
use rug::integer::{IsPrime, Order};
use rug::ops::{NegAssign, Pow, RemRounding};
use rug::{Assign, Complete, Integer};
use std::cmp::Ordering;
use std::sync::OnceLock;

/// Arbitrary-precision signed integer.
pub type Ibz = Integer;

// ---------------------------------------------------------------------------
// Constants

fn const_int(v: u32) -> &'static Ibz {
    // Each constant gets its own OnceLock; the array index picks it.
    static CELLS: [OnceLock<Ibz>; 4] = [
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
    ];
    CELLS[v as usize].get_or_init(|| Integer::from(v))
}

pub fn ibz_const_zero() -> &'static Ibz {
    const_int(0)
}
pub fn ibz_const_one() -> &'static Ibz {
    const_int(1)
}
pub fn ibz_const_two() -> &'static Ibz {
    const_int(2)
}
pub fn ibz_const_three() -> &'static Ibz {
    const_int(3)
}

// ---------------------------------------------------------------------------
// Init / finalize (no-ops in Rust beyond Default; kept for C parity)

#[inline]
pub fn ibz_init(x: &mut Ibz) {
    x.assign(0);
}
#[inline]
pub fn ibz_finalize(_x: &mut Ibz) {}

// ---------------------------------------------------------------------------
// Basic arithmetic

#[inline]
pub fn ibz_add(sum: &mut Ibz, a: &Ibz, b: &Ibz) {
    sum.assign(a + b);
}
#[inline]
pub fn ibz_sub(diff: &mut Ibz, a: &Ibz, b: &Ibz) {
    diff.assign(a - b);
}
#[inline]
pub fn ibz_mul(prod: &mut Ibz, a: &Ibz, b: &Ibz) {
    prod.assign(a * b);
}
#[inline]
pub fn ibz_neg(neg: &mut Ibz, a: &Ibz) {
    neg.assign(a);
    neg.neg_assign();
}
#[inline]
pub fn ibz_abs(abs: &mut Ibz, a: &Ibz) {
    abs.assign(a.abs_ref());
}

/// Truncating division: `a = quotient*b + remainder`, quotient rounded toward zero
/// (matches `mpz_tdiv_qr`).
pub fn ibz_div(quotient: &mut Ibz, remainder: &mut Ibz, a: &Ibz, b: &Ibz) {
    let (q, r) = a.div_rem_ref(b).complete();
    quotient.assign(q);
    remainder.assign(r);
}

/// Truncating right-shift of |a| with sign preserved (matches `mpz_tdiv_q_2exp`).
#[inline]
pub fn ibz_div_2exp(quotient: &mut Ibz, a: &Ibz, exp: u32) {
    // mpz_tdiv_q_2exp rounds toward zero; rug's >> on Integer is arithmetic
    // (floor for negatives), so emulate via abs/shift/restore-sign.
    let neg = a.is_negative();
    quotient.assign(a.abs_ref());
    *quotient >>= exp;
    if neg {
        quotient.neg_assign();
    }
}

/// Floor division (matches `mpz_fdiv_qr`).
pub fn ibz_div_floor(q: &mut Ibz, r: &mut Ibz, n: &Ibz, d: &Ibz) {
    let (qq, rr) = n.div_rem_floor_ref(d).complete();
    q.assign(qq);
    r.assign(rr);
}

/// `r = a mod |b|`, always non-negative (matches `mpz_mod`).
pub fn ibz_mod(r: &mut Ibz, a: &Ibz, b: &Ibz) {
    let abs_b = b.clone().abs();
    r.assign(a.rem_floor(&abs_b));
}

/// Floor remainder by an unsigned long (matches `mpz_fdiv_ui`).
#[inline]
pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    debug_assert!(d != 0);
    // rug's mod_u takes u32; use a temporary for larger moduli.
    if d <= u32::MAX as u64 {
        n.mod_u(d as u32) as u64
    } else {
        let dm = Integer::from(d);
        let r: Integer = n.rem_floor(&dm).into();
        r.to_u64().unwrap()
    }
}

#[inline]
pub fn ibz_divides(a: &Ibz, b: &Ibz) -> i32 {
    a.is_divisible(b) as i32
}

#[inline]
pub fn ibz_pow(out: &mut Ibz, x: &Ibz, e: u32) {
    out.assign(x.pow(e));
}

pub fn ibz_pow_mod(out: &mut Ibz, x: &Ibz, e: &Ibz, m: &Ibz) {
    out.assign(x);
    out.pow_mod_mut(e, m)
        .expect("ibz_pow_mod: x not invertible for negative exponent");
}

/// Position of the lowest set bit (`mpz_scan1(x, 0)`).
#[inline]
pub fn ibz_two_adic(x: &Ibz) -> i32 {
    x.find_one(0).map(|b| b as i32).unwrap_or(-1)
}

#[inline]
pub fn ibz_cmp(a: &Ibz, b: &Ibz) -> i32 {
    match a.cmp(b) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}
#[inline]
pub fn ibz_is_zero(x: &Ibz) -> i32 {
    x.is_zero() as i32
}
#[inline]
pub fn ibz_is_one(x: &Ibz) -> i32 {
    (*x == 1) as i32
}
#[inline]
pub fn ibz_cmp_int32(x: &Ibz, y: i32) -> i32 {
    match x.partial_cmp(&y).unwrap() {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}
#[inline]
pub fn ibz_is_even(x: &Ibz) -> i32 {
    (!x.get_bit(0)) as i32
}
#[inline]
pub fn ibz_is_odd(x: &Ibz) -> i32 {
    x.get_bit(0) as i32
}

#[inline]
pub fn ibz_set(i: &mut Ibz, x: i32) {
    i.assign(x);
}

pub fn ibz_convert_to_str(i: &Ibz, base: i32) -> Option<String> {
    if base != 10 && base != 16 {
        return None;
    }
    Some(i.to_string_radix(base))
}

#[inline]
pub fn ibz_print(num: &Ibz, base: i32) {
    debug_assert!(base == 10 || base == 16);
    print!("{}", num.to_string_radix(base));
}

pub fn ibz_set_from_str(i: &mut Ibz, s: &str, base: i32) -> i32 {
    match Integer::from_str_radix(s, base) {
        Ok(v) => {
            i.assign(v);
            1
        }
        Err(_) => 0,
    }
}

#[inline]
pub fn ibz_copy(target: &mut Ibz, value: &Ibz) {
    target.assign(value);
}
#[inline]
pub fn ibz_swap(a: &mut Ibz, b: &mut Ibz) {
    core::mem::swap(a, b);
}

/// Low 31 bits + sign bit, matching the C `LONG_MAX > INT32_MAX` branch.
/// Calls `mpz_get_si` directly to exactly match C's truncation semantics
/// (sign+magnitude, not two's-complement) for KAT determinism.
pub fn ibz_get(i: &Ibz) -> i32 {
    // SAFETY: rug's Integer is a transparent wrapper over mpz_t.
    let t: i64 = unsafe { gmp_mpfr_sys::gmp::mpz_get_si(i.as_raw()) };
    let sign_bit = ((t >> 32) as i32) & i32::MIN;
    let low = (t as i32) & i32::MAX;
    sign_bit | low
}

// ---------------------------------------------------------------------------
// Digit / byte conversion

/// Read `dig_len` little-endian u64 limbs into `target` (matches `mpz_import` call in C).
pub fn ibz_copy_digits(target: &mut Ibz, dig: &[Digit]) {
    target.assign_digits(dig, Order::Lsf);
}

/// Write `ibz` to `target` as little-endian u64 limbs, zero-padding to len.
/// `ibz` must be non-negative.
pub fn ibz_to_digits(target: &mut [Digit], ibz: &Ibz) {
    debug_assert!(!ibz.is_negative());
    for d in target.iter_mut() {
        *d = 0;
    }
    let n = ibz.significant_digits::<Digit>();
    debug_assert!(n <= target.len(), "target too small for ibz");
    if n > 0 {
        ibz.write_digits(&mut target[..n], Order::Lsf);
    }
}

// ---------------------------------------------------------------------------
// Random sampling (deterministic via the KAT DRBG)

/// Uniform sample in `[a, b]` by rejection. Matches the C byte-sampling loop
/// (length and mask computed from `b - a`, big-endian masking on the top limb).
pub fn ibz_rand_interval(rand: &mut Ibz, a: &Ibz, b: &Ibz) -> i32 {
    let mut bmina = (b - a).complete();
    if bmina.is_zero() {
        rand.assign(a);
        return 1;
    }
    debug_assert!(bmina.is_positive());

    let len_bits = bmina.significant_bits() as usize;
    let len_bytes = (len_bits + 7) / 8;
    const LIMB: usize = core::mem::size_of::<Digit>();
    const LIMB_BITS: usize = LIMB * 8;
    let len_limbs = (len_bytes + LIMB - 1) / LIMB;

    let shift = (LIMB_BITS - len_bits % LIMB_BITS) % LIMB_BITS;
    let mask: Digit = (!0u64) >> shift;

    let mut buf = vec![0u8; len_limbs * LIMB];
    let mut tmp = Integer::new();
    loop {
        // C calls randombytes for len_bytes only; the trailing bytes of the
        // top limb stay whatever was there before. We zero the buffer once
        // (C's stack array is uninitialized in release, but the mask renders
        // those bits irrelevant — and in C's debug path the bytes are pre-set
        // to 0xFF then overwritten; the mask still discards them). To match
        // KAT bytes consumed, we randombytes exactly len_bytes.
        randombytes(&mut buf[..len_bytes]);
        // Reinterpret as little-endian limbs.
        let mut limbs = vec![0u64; len_limbs];
        for (i, chunk) in buf.chunks(LIMB).enumerate() {
            let mut arr = [0u8; LIMB];
            arr[..chunk.len()].copy_from_slice(chunk);
            limbs[i] = u64::from_le_bytes(arr);
        }
        limbs[len_limbs - 1] &= mask;
        tmp.assign_digits(&limbs, Order::Lsf);
        if tmp <= bmina {
            break;
        }
    }
    bmina.assign(&tmp + a);
    rand.assign(bmina);
    1
}

/// Uniform sample in `[a, b]` for small non-negative `a < b`.
pub fn ibz_rand_interval_i(rand: &mut Ibz, a: i32, b: i32) -> i32 {
    assert!(a >= 0 && b >= 0 && b > a, "a={a} b={b}");
    let diff = (b - a) as u32;
    let lz = diff.leading_zeros();
    let mask: u32 = if lz == 0 {
        u32::MAX
    } else {
        (1u32 << (32 - lz)) - 1
    };
    debug_assert!(mask >= diff && (mask as u64) < 2 * diff as u64);

    let mut buf = [0u8; 4];
    let r = loop {
        randombytes(&mut buf);
        let r = u32::from_le_bytes(buf) & mask;
        if r <= diff {
            break r as i32;
        }
    };
    ibz_set(rand, r + a);
    1
}

/// Uniform sample in `[-m, m]`.
pub fn ibz_rand_interval_minm_m(rand: &mut Ibz, m: i32) -> i32 {
    let two_m = Integer::from(m) * 2;
    let ret = ibz_rand_interval(rand, ibz_const_zero(), &two_m);
    *rand -= m;
    ret
}

/// Uniform sample in `[-2^m, 2^m]`.
pub fn ibz_rand_interval_bits(rand: &mut Ibz, m: u32) -> i32 {
    let hi = Integer::from(1) << m;
    let lo = (-&hi).complete();
    let ret = ibz_rand_interval(rand, &lo, &hi);
    if ret != 1 {
        return ret;
    }
    // Match the C code's trailing `mpz_sub_ui(*rand, *rand, m)` exactly.
    *rand -= m;
    ret
}

#[inline]
pub fn ibz_bitsize(a: &Ibz) -> i32 {
    a.significant_bits() as i32
}
#[inline]
pub fn ibz_size_in_base(a: &Ibz, base: i32) -> i32 {
    // GMP's mpz_sizeinbase returns 1 for 0; rug's to_string_radix gives "0".
    // significant_bits handles base 2; for others, fall back to string length.
    if base == 2 {
        a.significant_bits().max(1) as i32
    } else {
        a.to_string_radix(base).trim_start_matches('-').len().max(1) as i32
    }
}

// ---------------------------------------------------------------------------
// Number theory

#[inline]
pub fn ibz_gcd(gcd: &mut Ibz, a: &Ibz, b: &Ibz) {
    gcd.assign(a.gcd_ref(b));
}

pub fn ibz_invmod(inv: &mut Ibz, a: &Ibz, m: &Ibz) -> i32 {
    inv.assign(a);
    match inv.invert_mut(m) {
        Ok(()) => 1,
        Err(()) => 0,
    }
}

pub fn ibz_probab_prime(n: &Ibz, reps: i32) -> i32 {
    match n.is_probably_prime(reps.max(0) as u32) {
        IsPrime::No => 0,
        IsPrime::Probably => 1,
        IsPrime::Yes => 2,
    }
}

#[inline]
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    a.legendre(p)
}

/// Perfect-square root: returns 1 and sets `sqrt` if `a` is a perfect square.
pub fn ibz_sqrt(sqrt: &mut Ibz, a: &Ibz) -> i32 {
    if a.is_perfect_square() {
        sqrt.assign(a.sqrt_ref());
        1
    } else {
        0
    }
}

#[inline]
pub fn ibz_sqrt_floor(sqrt: &mut Ibz, a: &Ibz) {
    sqrt.assign(a.sqrt_ref());
}

/// Tonelli–Shanks square root modulo an odd prime `p`.
/// Returns 1 and sets `sqrt` to a root of `a` mod `p` if one exists.
pub fn ibz_sqrt_mod_p(sqrt: &mut Ibz, a: &Ibz, p: &Ibz) -> i32 {
    debug_assert!(ibz_probab_prime(p, 30) != 0);

    let amod = {
        let mut t = Integer::new();
        ibz_mod(&mut t, a, p);
        t
    };
    if amod.is_zero() {
        sqrt.assign(0);
        // C falls through but legendre(0,p)=0 so it returns 0; we mirror that.
    }
    if amod.legendre(p) != 1 {
        return 0;
    }

    let pm1: Integer = (p - 1u32).complete();
    let mut tmp = Integer::new();

    if p.mod_u(4) == 3 {
        // sqrt = a^{(p+1)/4}
        tmp.assign((p + 1u32).complete() >> 2);
        sqrt.assign(amod.pow_mod_ref(&tmp, p).unwrap());
    } else if p.mod_u(8) == 5 {
        // Atkin's trick.
        tmp.assign(&pm1 >> 2);
        let t = Integer::from(amod.pow_mod_ref(&tmp, p).unwrap());
        if t == 1 {
            tmp.assign((p + 3u32).complete() >> 3);
            sqrt.assign(amod.pow_mod_ref(&tmp, p).unwrap());
        } else {
            tmp.assign((p - 5u32).complete() >> 3);
            let a4: Integer = (&amod << 2u32).complete();
            let mut r = Integer::from(a4.pow_mod_ref(&tmp, p).unwrap());
            let a2: Integer = (&amod << 1u32).complete();
            r *= &a2;
            ibz_mod(sqrt, &r, p);
        }
    } else {
        // p ≡ 1 (mod 8): full Tonelli–Shanks.
        let mut q: Integer = pm1.clone();
        let mut e: u32 = 0;
        while !q.get_bit(e) {
            e += 1;
        }
        q >>= e;

        // Find a non-residue.
        let mut qnr = Integer::from(2);
        while qnr.legendre(p) != -1 {
            qnr += 1;
        }
        let mut z = Integer::from(qnr.pow_mod_ref(&q, p).unwrap());

        let mut y = Integer::from(amod.pow_mod_ref(&q, p).unwrap());
        let half_q1: Integer = (&q + 1u32).complete() >> 1;
        let mut x = Integer::from(amod.pow_mod_ref(&half_q1, p).unwrap());

        let mut exp = Integer::from(1) << (e - 2);
        let two = Integer::from(2);
        for _ in 0..e {
            let b = Integer::from(y.pow_mod_ref(&exp, p).unwrap());
            if b == pm1 {
                x *= &z;
                ibz_mod(&mut tmp, &x, p);
                x.assign(&tmp);
                y *= &z;
                y *= &z;
                ibz_mod(&mut tmp, &y, p);
                y.assign(&tmp);
            }
            // z = z^2 mod p
            z.pow_mod_mut(&two, p).unwrap();
            exp >>= 1;
        }
        sqrt.assign(x);
    }
    1
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }

    #[test]
    fn arithmetic_basics() {
        let mut r = Ibz::new();
        ibz_add(&mut r, &z(7), &z(-3));
        assert_eq!(r, 4);
        ibz_sub(&mut r, &z(7), &z(-3));
        assert_eq!(r, 10);
        ibz_mul(&mut r, &z(7), &z(-3));
        assert_eq!(r, -21);
        ibz_neg(&mut r, &z(-5));
        assert_eq!(r, 5);
        ibz_abs(&mut r, &z(-9));
        assert_eq!(r, 9);
    }

    #[test]
    fn div_tdiv_semantics() {
        // mpz_tdiv_qr rounds toward zero.
        let (mut q, mut r) = (Ibz::new(), Ibz::new());
        ibz_div(&mut q, &mut r, &z(-7), &z(3));
        assert_eq!((q.to_i32().unwrap(), r.to_i32().unwrap()), (-2, -1));
        ibz_div(&mut q, &mut r, &z(7), &z(-3));
        assert_eq!((q.to_i32().unwrap(), r.to_i32().unwrap()), (-2, 1));
    }

    #[test]
    fn div_floor_semantics() {
        let (mut q, mut r) = (Ibz::new(), Ibz::new());
        ibz_div_floor(&mut q, &mut r, &z(-7), &z(3));
        assert_eq!((q.to_i32().unwrap(), r.to_i32().unwrap()), (-3, 2));
    }

    #[test]
    fn div_2exp_trunc() {
        let mut q = Ibz::new();
        ibz_div_2exp(&mut q, &z(-7), 1);
        assert_eq!(q, -3); // tdiv: -7/2 → -3, not -4
        ibz_div_2exp(&mut q, &z(7), 1);
        assert_eq!(q, 3);
    }

    #[test]
    fn mod_nonneg() {
        let mut r = Ibz::new();
        ibz_mod(&mut r, &z(-7), &z(3));
        assert_eq!(r, 2);
        ibz_mod(&mut r, &z(-7), &z(-3));
        assert_eq!(r, 2); // sign of divisor ignored
        ibz_mod(&mut r, &z(7), &z(3));
        assert_eq!(r, 1);
    }

    #[test]
    fn predicates() {
        assert_eq!(ibz_is_zero(&z(0)), 1);
        assert_eq!(ibz_is_one(&z(1)), 1);
        assert_eq!(ibz_is_even(&z(6)), 1);
        assert_eq!(ibz_is_odd(&z(7)), 1);
        assert_eq!(ibz_cmp(&z(5), &z(3)), 1);
        assert_eq!(ibz_cmp_int32(&z(-2), -2), 0);
        assert_eq!(ibz_divides(&z(12), &z(4)), 1);
        assert_eq!(ibz_divides(&z(12), &z(5)), 0);
    }

    #[test]
    fn two_adic_and_bitsize() {
        assert_eq!(ibz_two_adic(&z(12)), 2); // 12 = 0b1100
        assert_eq!(ibz_two_adic(&z(1)), 0);
        assert_eq!(ibz_bitsize(&z(255)), 8);
        assert_eq!(ibz_bitsize(&z(256)), 9);
    }

    #[test]
    fn ibz_get_low_bits() {
        // Preserves parity and low bits even for large values.
        let big = Integer::from_str_radix("1234567890ABCDEF1234567890ABCDEF", 16).unwrap();
        let g = ibz_get(&big);
        assert_eq!(g & 1, (big.get_bit(0) as i32));
        assert_eq!(g & 3, (big.mod_u(4) as i32));
        // For values fitting in i32, exact.
        assert_eq!(ibz_get(&z(42)), 42);
        assert_eq!(ibz_get(&z(-42)), -42);
    }

    #[test]
    fn digits_roundtrip() {
        let limbs: [u64; 3] = [0xDEADBEEF, 0xCAFEBABE, 0x1];
        let mut x = Ibz::new();
        ibz_copy_digits(&mut x, &limbs);
        let mut out = [0u64; 4];
        ibz_to_digits(&mut out, &x);
        assert_eq!(&out[..3], &limbs);
        assert_eq!(out[3], 0);
        // zero
        let zero = Ibz::new();
        let mut out2 = [99u64; 2];
        ibz_to_digits(&mut out2, &zero);
        assert_eq!(out2, [0, 0]);
    }

    #[test]
    fn gcd_invmod() {
        let mut g = Ibz::new();
        ibz_gcd(&mut g, &z(462), &z(1071));
        assert_eq!(g, 21);
        let mut inv = Ibz::new();
        assert_eq!(ibz_invmod(&mut inv, &z(3), &z(7)), 1);
        assert_eq!(inv, 5);
        assert_eq!(ibz_invmod(&mut inv, &z(4), &z(8)), 0);
    }

    #[test]
    fn sqrt_perfect() {
        let mut s = Ibz::new();
        assert_eq!(ibz_sqrt(&mut s, &z(144)), 1);
        assert_eq!(s, 12);
        assert_eq!(ibz_sqrt(&mut s, &z(145)), 0);
        ibz_sqrt_floor(&mut s, &z(145));
        assert_eq!(s, 12);
    }

    #[test]
    fn sqrt_mod_p_all_branches() {
        let mut s = Ibz::new();
        // p ≡ 3 mod 4: p=7, a=2 → 4²=16≡2
        assert_eq!(ibz_sqrt_mod_p(&mut s, &z(2), &z(7)), 1);
        let s2 = (&s * &s).complete();
        assert_eq!(s2.mod_u(7), 2);
        // p ≡ 5 mod 8: p=13, a=3 → 4²=16≡3
        assert_eq!(ibz_sqrt_mod_p(&mut s, &z(3), &z(13)), 1);
        let s2 = (&s * &s).complete();
        assert_eq!(s2.mod_u(13), 3);
        // p=13, a=10 (the t≠1 sub-branch): 6²=36≡10
        assert_eq!(ibz_sqrt_mod_p(&mut s, &z(10), &z(13)), 1);
        let s2 = (&s * &s).complete();
        assert_eq!(s2.mod_u(13), 10);
        // p ≡ 1 mod 8: p=17, a=2 → 6²=36≡2
        assert_eq!(ibz_sqrt_mod_p(&mut s, &z(2), &z(17)), 1);
        let s2 = (&s * &s).complete();
        assert_eq!(s2.mod_u(17), 2);
        // Non-residue: a=3 mod 7
        assert_eq!(ibz_sqrt_mod_p(&mut s, &z(3), &z(7)), 0);
    }

    #[test]
    fn legendre_prime() {
        assert_eq!(ibz_legendre(&z(2), &z(7)), 1);
        assert_eq!(ibz_legendre(&z(3), &z(7)), -1);
        assert!(ibz_probab_prime(&z(97), 20) > 0);
        assert_eq!(ibz_probab_prime(&z(100), 20), 0);
    }

    #[test]
    fn pow_and_pow_mod() {
        let mut r = Ibz::new();
        ibz_pow(&mut r, &z(3), 5);
        assert_eq!(r, 243);
        ibz_pow_mod(&mut r, &z(3), &z(5), &z(7));
        assert_eq!(r, 5); // 243 mod 7 = 5
    }

    #[test]
    fn str_roundtrip() {
        let mut x = Ibz::new();
        assert_eq!(ibz_set_from_str(&mut x, "-12345678901234567890", 10), 1);
        assert_eq!(ibz_convert_to_str(&x, 10).unwrap(), "-12345678901234567890");
    }
}
