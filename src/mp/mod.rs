// SPDX-License-Identifier: Apache-2.0
//! Saturated-limb fixed-precision multiprecision arithmetic.
//!
//! Direct port of `src/mp/ref/generic/{mp.c, include/mp.h}` from the C
//! reference. Operates on little-endian arrays of 64-bit digits.
//!
//! All functions preserve the exact semantics (including quirks) of the C
//! implementation so that downstream modules produce bit-identical results.

#![allow(clippy::needless_range_loop)]

use core::fmt::Write as _;

/// One saturated limb. The C reference is built with `RADIX_64`.
pub type Digit = u64;
/// Signed companion of [`Digit`].
pub type SDigit = i64;
/// Bits per limb.
pub const RADIX: u32 = 64;
/// log2(RADIX).
pub const LOG2RADIX: u32 = 6;

// ---------------------------------------------------------------------------
// Constant-time digit predicates (from mp.h)
// ---------------------------------------------------------------------------

/// Returns 1 if `x != 0`, else 0. Constant-time.
#[inline(always)]
pub const fn is_digit_nonzero_ct(x: Digit) -> u32 {
    ((x | x.wrapping_neg()) >> (RADIX - 1)) as u32
}

/// Returns 1 if `x == 0`, else 0. Constant-time.
#[inline(always)]
pub const fn is_digit_zero_ct(x: Digit) -> u32 {
    1 ^ is_digit_nonzero_ct(x)
}

/// Returns 1 if `x < y`, else 0. Constant-time.
#[inline(always)]
pub const fn is_digit_lessthan_ct(x: Digit, y: Digit) -> u32 {
    ((x ^ ((x ^ y) | (x.wrapping_sub(y) ^ y))) >> (RADIX - 1)) as u32
}

// ---------------------------------------------------------------------------
// Digit-level primitives (ADDC / SUBC / SHIFTR / SHIFTL / MUL macros in C)
// ---------------------------------------------------------------------------

/// Full-adder: returns `(sum, carry_out)` where
/// `sum = a + b + carry_in (mod 2^64)` and `carry_out ∈ {0,1}`.
#[inline(always)]
pub const fn addc(a: Digit, b: Digit, carry_in: u32) -> (Digit, u32) {
    let temp = a.wrapping_add(carry_in as Digit);
    let sum = b.wrapping_add(temp);
    let carry_out = is_digit_lessthan_ct(temp, carry_in as Digit) | is_digit_lessthan_ct(sum, temp);
    (sum, carry_out)
}

/// Full-subtractor: returns `(diff, borrow_out)` where
/// `diff = a - b - borrow_in (mod 2^64)` and `borrow_out ∈ {0,1}`.
#[inline(always)]
pub const fn subc(a: Digit, b: Digit, borrow_in: u32) -> (Digit, u32) {
    let temp = a.wrapping_sub(b);
    let borrow_out = is_digit_lessthan_ct(a, b) | ((borrow_in) & is_digit_zero_ct(temp));
    let diff = temp.wrapping_sub(borrow_in as Digit);
    (diff, borrow_out)
}

/// Funnel shift right: `(high:low) >> shift`, low word. `shift` must be in `1..RADIX`.
#[inline(always)]
pub const fn shiftr(high: Digit, low: Digit, shift: u32) -> Digit {
    (low >> shift) ^ (high << (RADIX - shift))
}

/// Funnel shift left: `(high:low) << shift`, high word. `shift` must be in `1..RADIX`.
#[inline(always)]
pub const fn shiftl(high: Digit, low: Digit, shift: u32) -> Digit {
    (high << shift) ^ (low >> (RADIX - shift))
}

/// 64×64 → 128-bit multiply. `out[0]` = low word, `out[1]` = high word.
#[inline(always)]
pub fn mul(out: &mut [Digit], a: Digit, b: Digit) {
    let r = (a as u128) * (b as u128);
    out[0] = r as Digit;
    out[1] = (r >> RADIX) as Digit;
}

/// Convenience wrapper around [`mul`] returning a `[low, high]` pair.
#[inline(always)]
pub fn mul_pair(a: Digit, b: Digit) -> [Digit; 2] {
    let r = (a as u128) * (b as u128);
    [r as Digit, (r >> RADIX) as Digit]
}

// ---------------------------------------------------------------------------
// Multiprecision routines (from mp.c)
// ---------------------------------------------------------------------------

/// `c = a + b` over `nwords` limbs (wrapping). `c` must not alias `a` or `b`.
#[inline]
pub fn mp_add(c: &mut [Digit], a: &[Digit], b: &[Digit], nwords: usize) {
    let mut carry = 0u32;
    for i in 0..nwords {
        let (s, co) = addc(a[i], b[i], carry);
        c[i] = s;
        carry = co;
    }
}

/// `a += b` over `nwords` limbs (wrapping). In-place variant used where the
/// C code calls `mp_add(x, x, b, n)`.
#[inline]
pub fn mp_add_inplace(a: &mut [Digit], b: &[Digit], nwords: usize) {
    let mut carry = 0u32;
    for i in 0..nwords {
        let (s, co) = addc(a[i], b[i], carry);
        a[i] = s;
        carry = co;
    }
}

/// `c = a - b` over `nwords` limbs (wrapping). `c` must not alias `a` or `b`.
#[inline]
pub fn mp_sub(c: &mut [Digit], a: &[Digit], b: &[Digit], nwords: usize) {
    let mut borrow = 0u32;
    for i in 0..nwords {
        let (d, bo) = subc(a[i], b[i], borrow);
        c[i] = d;
        borrow = bo;
    }
}

/// In-place right shift by `shift ∈ 1..RADIX`. Returns the original LSB
/// (matching the C, which returns `x[0] & 1` regardless of `shift`).
#[inline]
pub fn mp_shiftr(x: &mut [Digit], shift: u32, nwords: usize) -> Digit {
    let bit_out = x[0] & 1;
    for i in 0..nwords - 1 {
        x[i] = shiftr(x[i + 1], x[i], shift);
    }
    x[nwords - 1] >>= shift;
    bit_out
}

/// In-place left shift by `shift ∈ 1..RADIX`.
#[inline]
pub fn mp_shiftl(x: &mut [Digit], shift: u32, nwords: usize) {
    for i in (1..nwords).rev() {
        x[i] = shiftl(x[i], x[i - 1], shift);
    }
    x[0] <<= shift;
}

/// Left shift by an arbitrary `shift` (may exceed `RADIX-1`). Shift of 0 is a
/// no-op (the C reference would hit UB via `>> 64` in `shiftl`).
#[inline]
pub fn multiple_mp_shiftl(x: &mut [Digit], shift: u32, nwords: usize) {
    if shift == 0 {
        return;
    }
    let mut t = shift as i64;
    while t > (RADIX - 1) as i64 {
        mp_shiftl(x, RADIX - 1, nwords);
        t -= (RADIX - 1) as i64;
    }
    mp_shiftl(x, t as u32, nwords);
}

/// Constant-time select: `c = (mask == 0) ? a : b` per limb, with
/// `mask ∈ {0, 0xFF..FF}`.
#[inline]
pub fn select_ct(c: &mut [Digit], a: &[Digit], b: &[Digit], mask: Digit, nwords: usize) {
    for i in 0..nwords {
        c[i] = ((a[i] ^ b[i]) & mask) ^ a[i];
    }
}

/// Constant-time conditional swap: if `option == 0xFF..FF` swap `a` and `b`,
/// if `option == 0` leave unchanged.
#[inline]
pub fn swap_ct(a: &mut [Digit], b: &mut [Digit], option: Digit, nwords: usize) {
    for i in 0..nwords {
        let temp = option & (a[i] ^ b[i]);
        a[i] ^= temp;
        b[i] ^= temp;
    }
}

/// Lexicographic compare (most-significant limb first).
/// Returns `1` if `a > b`, `0` if equal, `-1` if `a < b`. **Not** constant-time.
#[inline]
pub fn mp_compare(a: &[Digit], b: &[Digit], nwords: usize) -> i32 {
    for i in (0..nwords).rev() {
        if a[i] > b[i] {
            return 1;
        } else if a[i] < b[i] {
            return -1;
        }
    }
    0
}

/// Constant-time zero check over `nwords` limbs.
#[inline]
pub fn mp_is_zero(a: &[Digit], nwords: usize) -> bool {
    let mut r: Digit = 0;
    for i in 0..nwords {
        r |= a[i];
    }
    is_digit_zero_ct(r) != 0
}

/// 2-limb × 2-limb → 4-limb multiply.
///
/// NOTE: Faithful port of C `mp_mul2`, which omits the `a[1]*b[0]` partial
/// product. This matches the reference implementation; correctness review of
/// this routine against the spec is tracked separately.
#[inline]
pub fn mp_mul2(c: &mut [Digit], a: &[Digit], b: &[Digit]) {
    let mut carry = 0u32;
    let t0 = mul_pair(a[0], b[0]);
    let t1 = mul_pair(a[0], b[1]);
    let mut t0 = t0;
    let mut t1 = t1;
    let (s, co) = addc(t0[1], t1[0], carry);
    t0[1] = s;
    carry = co;
    let (s, co) = addc(0, t1[1], carry);
    t1[1] = s;
    carry = co;
    let mut t2 = mul_pair(a[1], b[1]);
    let (s, co) = addc(t2[0], t1[1], carry);
    t2[0] = s;
    carry = co;
    let (s, _co) = addc(0, t2[1], carry);
    t2[1] = s;
    c[0] = t0[0];
    c[1] = t0[1];
    c[2] = t2[0];
    c[3] = t2[1];
}

/// Debug print as a single big-endian hex literal.
pub fn mp_print(a: &[Digit], nwords: usize) -> String {
    let mut s = String::from("0x");
    for i in 0..nwords {
        let _ = write!(s, "{:016x}", a[nwords - 1 - i]);
    }
    s
}

/// `b = a` over `nwords` limbs.
#[inline]
pub fn mp_copy(b: &mut [Digit], a: &[Digit], nwords: usize) {
    b[..nwords].copy_from_slice(&a[..nwords]);
}

/// Low half of `a * b`: `nwords`-limb inputs, `nwords`-limb output.
/// Matches C `mp_mul`, which discards the high half.
pub fn mp_mul(c: &mut [Digit], a: &[Digit], b: &[Digit], nwords: usize) {
    let mut cc = vec![0 as Digit; nwords];
    let mut t = vec![0 as Digit; nwords];

    for i in 0..nwords {
        let p0 = mul_pair(a[i], b[0]);
        t[0] = p0[0];
        if nwords >= 2 {
            t[1] = p0[1];
        }

        for j in 1..nwords.saturating_sub(1) {
            let uv = mul_pair(a[i], b[j]);
            let (s, carry) = addc(t[j], uv[0], 0);
            t[j] = s;
            t[j + 1] = uv[1].wrapping_add(carry as Digit);
        }

        let j = nwords - 1;
        let uv = mul_pair(a[i], b[j]);
        let (s, _carry) = addc(t[j], uv[0], 0);
        t[j] = s;

        // C: mp_add(&cc[i], &cc[i], t, nwords - i)  — in-place on cc tail
        mp_add_inplace(&mut cc[i..], &t, nwords - i);
    }

    mp_copy(c, &cc, nwords);
}

/// In-place reduce `a` modulo `2^e` (zero out bits `>= e`).
#[inline]
pub fn mp_mod_2exp(a: &mut [Digit], e: u32, nwords: usize) {
    let q = (e >> LOG2RADIX) as usize;
    let r = e & (RADIX - 1);
    if q < nwords {
        a[q] &= ((1 as Digit) << r).wrapping_sub(1);
        for i in q + 1..nwords {
            a[i] = 0;
        }
    }
}

/// In-place two's-complement negate: `a = ~a; a[0] += 1`.
///
/// NOTE: Faithful port of C `mp_neg`, which does **not** propagate the +1
/// carry past limb 0. Correct for the callers that immediately reduce mod `2^e`.
#[inline]
pub fn mp_neg(a: &mut [Digit], nwords: usize) {
    for i in 0..nwords {
        a[i] ^= Digit::MAX;
    }
    a[0] = a[0].wrapping_add(1);
}

/// Returns `true` iff `x == 1` over `nwords` limbs. Not constant-time.
#[inline]
pub fn mp_is_one(x: &[Digit], nwords: usize) -> bool {
    if x[0] != 1 {
        return false;
    }
    for i in 1..nwords {
        if x[i] != 0 {
            return false;
        }
    }
    true
}

/// Returns `true` iff limb 0 is odd.
#[inline(always)]
pub fn mp_is_odd(x: &[Digit]) -> bool {
    (x[0] & 1) != 0
}

/// Returns `true` iff limb 0 is even.
#[inline(always)]
pub fn mp_is_even(x: &[Digit]) -> bool {
    !mp_is_odd(x)
}

/// Compute `b = a^{-1} mod 2^e` via Newton/Hensel lifting. `a` must be odd.
pub fn mp_inv_2e(b: &mut [Digit], a: &[Digit], e: i32, nwords: usize) {
    debug_assert!(a[0] & 1 == 1, "mp_inv_2e: input must be odd");

    let mut x = vec![0 as Digit; nwords];
    let mut y = vec![0 as Digit; nwords];
    let mut aa = vec![0 as Digit; nwords];
    let mut tmp = vec![0 as Digit; nwords];
    let mut mp_one = vec![0 as Digit; nwords];
    mp_one[0] = 1;

    mp_copy(&mut aa, a, nwords);

    let mut p: i32 = 1;
    while (1i32 << p) < e {
        p += 1;
    }
    p -= 2;
    let w = (1i32 << (p + 2)) as u32;

    mp_mod_2exp(&mut aa, w, nwords);
    mp_add(&mut x, &aa, &aa, nwords);
    mp_add_inplace(&mut x, &aa, nwords); // 3a
    x[0] ^= 1 << 1; // (3a) xor 2
    mp_mod_2exp(&mut x, w, nwords); // x*a == 1 mod 2^4

    mp_mul(&mut tmp, &aa, &x, nwords);
    mp_neg(&mut tmp, nwords);
    mp_add(&mut y, &mp_one, &tmp, nwords);

    let mut xtmp = vec![0 as Digit; nwords];
    let mut ytmp = vec![0 as Digit; nwords];
    for _ in 0..p {
        mp_add(&mut tmp, &mp_one, &y, nwords);
        mp_copy(&mut xtmp, &x, nwords);
        mp_mul(&mut x, &xtmp, &tmp, nwords);
        mp_copy(&mut ytmp, &y, nwords);
        mp_mul(&mut y, &ytmp, &ytmp, nwords);
    }

    mp_mod_2exp(&mut x, w, nwords);
    mp_copy(b, &x, nwords);

    #[cfg(debug_assertions)]
    {
        let mut check = vec![0 as Digit; nwords];
        mp_mul(&mut check, &x, &aa, nwords);
        mp_mod_2exp(&mut check, w, nwords);
        debug_assert!(mp_is_one(&check, nwords), "mp_inv_2e self-check failed");
    }
}

/// Invert the 2×2 matrix `((r1,r2),(s1,s2))` modulo `2^e` in place.
/// Determinant must be odd.
pub fn mp_invert_matrix(
    r1: &mut [Digit],
    r2: &mut [Digit],
    s1: &mut [Digit],
    s2: &mut [Digit],
    e: i32,
    nwords: usize,
) {
    let mut p: i32 = 1;
    while (1i32 << p) < e {
        p += 1;
    }
    let w = (1i32 << p) as u32;

    let mut det = vec![0 as Digit; nwords];
    let mut tmp = vec![0 as Digit; nwords];
    let mut resa = vec![0 as Digit; nwords];
    let mut resb = vec![0 as Digit; nwords];
    let mut resc = vec![0 as Digit; nwords];
    let mut resd = vec![0 as Digit; nwords];

    mp_mul(&mut tmp, r1, s2, nwords);
    mp_mul(&mut det, r2, s1, nwords);
    let det_copy = det.clone();
    mp_sub(&mut det, &tmp, &det_copy, nwords);
    let det_copy = det.clone();
    mp_inv_2e(&mut det, &det_copy, e, nwords);

    mp_mul(&mut resa, &det, s2, nwords);
    mp_mul(&mut resb, &det, r2, nwords);
    mp_mul(&mut resc, &det, s1, nwords);
    mp_mul(&mut resd, &det, r1, nwords);

    mp_neg(&mut resb, nwords);
    mp_neg(&mut resc, nwords);

    mp_mod_2exp(&mut resa, w, nwords);
    mp_mod_2exp(&mut resb, w, nwords);
    mp_mod_2exp(&mut resc, w, nwords);
    mp_mod_2exp(&mut resd, w, nwords);

    mp_copy(r1, &resa, nwords);
    mp_copy(r2, &resb, nwords);
    mp_copy(s1, &resc, nwords);
    mp_copy(s2, &resd, nwords);
}

/// Byte-swap a digit (matches C `BSWAP_DIGIT`).
#[inline(always)]
pub const fn bswap_digit(x: Digit) -> Digit {
    x.swap_bytes()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addc_basic() {
        assert_eq!(addc(1, 2, 0), (3, 0));
        assert_eq!(addc(Digit::MAX, 1, 0), (0, 1));
        assert_eq!(addc(Digit::MAX, 0, 1), (0, 1));
        assert_eq!(addc(Digit::MAX, Digit::MAX, 1), (Digit::MAX, 1));
    }

    #[test]
    fn subc_basic() {
        assert_eq!(subc(3, 2, 0), (1, 0));
        assert_eq!(subc(0, 1, 0), (Digit::MAX, 1));
        assert_eq!(subc(0, 0, 1), (Digit::MAX, 1));
        assert_eq!(subc(5, 3, 1), (1, 0));
    }

    #[test]
    fn mp_add_carry_propagates() {
        let a = [Digit::MAX, Digit::MAX, 0, 0];
        let b = [1, 0, 0, 0];
        let mut c = [0; 4];
        mp_add(&mut c, &a, &b, 4);
        assert_eq!(c, [0, 0, 1, 0]);
    }

    #[test]
    fn mp_sub_borrow_propagates() {
        let a = [0, 0, 1, 0];
        let b = [1, 0, 0, 0];
        let mut c = [0; 4];
        mp_sub(&mut c, &a, &b, 4);
        assert_eq!(c, [Digit::MAX, Digit::MAX, 0, 0]);
    }

    #[test]
    fn mul_full_width() {
        let mut out = [0; 2];
        mul(&mut out, Digit::MAX, Digit::MAX);
        // (2^64 - 1)^2 = 2^128 - 2^65 + 1
        assert_eq!(out, [1, Digit::MAX - 1]);
    }

    #[test]
    fn shifts() {
        let mut x = [0b1010u64, 0b11u64];
        mp_shiftl(&mut x, 1, 2);
        assert_eq!(x, [0b10100, 0b110]);
        let bit = mp_shiftr(&mut x, 1, 2);
        assert_eq!(bit, 0);
        assert_eq!(x, [0b1010, 0b11]);

        let mut y = [1u64, 0u64];
        let bit = mp_shiftr(&mut y, 1, 2);
        assert_eq!(bit, 1);
        assert_eq!(y, [0, 0]);
    }

    #[test]
    fn multiple_shiftl_crosses_limbs() {
        let mut x = [1u64, 0, 0];
        multiple_mp_shiftl(&mut x, 65, 3);
        assert_eq!(x, [0, 2, 0]);
    }

    #[test]
    fn compare_and_zero() {
        assert_eq!(mp_compare(&[2, 0], &[1, 0], 2), 1);
        assert_eq!(mp_compare(&[1, 0], &[2, 0], 2), -1);
        assert_eq!(mp_compare(&[5, 7], &[5, 7], 2), 0);
        assert_eq!(mp_compare(&[Digit::MAX, 0], &[0, 1], 2), -1);
        assert!(mp_is_zero(&[0, 0, 0], 3));
        assert!(!mp_is_zero(&[0, 1, 0], 3));
        assert!(mp_is_one(&[1, 0, 0], 3));
        assert!(!mp_is_one(&[1, 1, 0], 3));
    }

    #[test]
    fn select_and_swap_ct() {
        let a = [1u64, 2];
        let b = [3u64, 4];
        let mut c = [0u64; 2];
        select_ct(&mut c, &a, &b, 0, 2);
        assert_eq!(c, a);
        select_ct(&mut c, &a, &b, Digit::MAX, 2);
        assert_eq!(c, b);

        let mut x = [1u64, 2];
        let mut y = [3u64, 4];
        swap_ct(&mut x, &mut y, 0, 2);
        assert_eq!((x, y), ([1, 2], [3, 4]));
        swap_ct(&mut x, &mut y, Digit::MAX, 2);
        assert_eq!((x, y), ([3, 4], [1, 2]));
    }

    #[test]
    fn mp_mul_low_half_matches_u128() {
        // 2-word: low half of product fits in 128 bits, compare against u128.
        let a = [0x1234_5678_9abc_def0u64, 0x0fed_cba9_8765_4321];
        let b = [0xdead_beef_cafe_babeu64, 0x0011_2233_4455_6677];
        let mut c = [0u64; 2];
        mp_mul(&mut c, &a, &b, 2);
        let aa = (a[1] as u128) << 64 | a[0] as u128;
        let bb = (b[1] as u128) << 64 | b[0] as u128;
        let expect = aa.wrapping_mul(bb);
        assert_eq!(c, [expect as u64, (expect >> 64) as u64]);
    }

    #[test]
    fn mod_2exp() {
        let mut a = [Digit::MAX; 3];
        mp_mod_2exp(&mut a, 70, 3);
        assert_eq!(a, [Digit::MAX, 0x3F, 0]);
    }

    #[test]
    fn inv_2e_roundtrip() {
        const N: usize = 4;
        let a: [Digit; N] = [0xCAFEBABE_DEADBEEF | 1, 0x12345678, 0, 0];
        let mut inv = [0 as Digit; N];
        mp_inv_2e(&mut inv, &a, 200, N);
        let mut prod = [0 as Digit; N];
        mp_mul(&mut prod, &a, &inv, N);
        mp_mod_2exp(&mut prod, 200, N);
        assert!(mp_is_one(&prod, N), "a * a^-1 mod 2^200 != 1: {:?}", prod);
    }

    #[test]
    fn invert_matrix_roundtrip() {
        const N: usize = 4;
        let e = 120;
        let mut r1: [Digit; N] = [3, 0, 0, 0];
        let mut r2: [Digit; N] = [5, 0, 0, 0];
        let mut s1: [Digit; N] = [7, 0, 0, 0];
        let mut s2: [Digit; N] = [12, 0, 0, 0]; // det = 3*12 - 5*7 = 1 (odd)
        let (a, b, c, d) = (r1, r2, s1, s2);
        mp_invert_matrix(&mut r1, &mut r2, &mut s1, &mut s2, e, N);

        // Multiply original * inverse, expect identity mod 2^e (rounded up to 2^128).
        let mul2 = |x: &[Digit; N], y: &[Digit; N]| -> [Digit; N] {
            let mut o = [0; N];
            mp_mul(&mut o, x, y, N);
            o
        };
        let add2 = |x: &[Digit; N], y: &[Digit; N]| -> [Digit; N] {
            let mut o = [0; N];
            mp_add(&mut o, x, y, N);
            o
        };
        let mut m00 = add2(&mul2(&a, &r1), &mul2(&b, &s1));
        let mut m01 = add2(&mul2(&a, &r2), &mul2(&b, &s2));
        let mut m10 = add2(&mul2(&c, &r1), &mul2(&d, &s1));
        let mut m11 = add2(&mul2(&c, &r2), &mul2(&d, &s2));
        for m in [&mut m00, &mut m01, &mut m10, &mut m11] {
            mp_mod_2exp(m, 128, N);
        }
        assert!(mp_is_one(&m00, N));
        assert!(mp_is_zero(&m01, N));
        assert!(mp_is_zero(&m10, N));
        assert!(mp_is_one(&m11, N));
    }

    #[test]
    fn parity() {
        assert!(mp_is_even(&[0u64, 0]));
        assert!(!mp_is_odd(&[0u64, 0]));
        assert!(mp_is_odd(&[1u64, 0]));
        assert!(!mp_is_even(&[1u64, 0]));
        assert!(mp_is_odd(&[u64::MAX, u64::MAX]));
        assert!(mp_is_even(&[2u64, 7]));
    }

    #[test]
    fn digit_nonzero_ct_top_bit_only() {
        // Regression for x | -x vs x ^ -x: when x = 2^63, -x = 2^63, so x ^ -x = 0.
        assert_eq!(is_digit_nonzero_ct(1u64 << 63), 1);
        assert_eq!(is_digit_nonzero_ct(0), 0);
        assert_eq!(is_digit_zero_ct(1u64 << 63), 0);
    }

    #[test]
    fn multiple_shiftl_exact_radix() {
        // Shift of exactly 64 must enter the loop once (63) then finish with 1.
        let mut x = [1u64, 0, 0];
        multiple_mp_shiftl(&mut x, 64, 3);
        assert_eq!(x, [0, 1, 0]);
        let mut x = [1u64, 0, 0];
        multiple_mp_shiftl(&mut x, 127, 3);
        assert_eq!(x, [0, 1u64 << 63, 0]);
    }

    #[test]
    fn multiple_shiftl_zero_is_nop() {
        // Previously panicked in debug via `>> 64` in the funnel shift.
        let mut x = [3u64, 5, 7];
        multiple_mp_shiftl(&mut x, 0, 3);
        assert_eq!(x, [3, 5, 7]);
    }
}
