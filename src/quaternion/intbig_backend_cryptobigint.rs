// SPDX-License-Identifier: Apache-2.0
//! `crypto-bigint` backend for `Ibz`. Selected with `feature = "cryptobigint"`.
//!
//! Unlike the malachite/rug backends this is *fixed-precision*: every `Ibz`
//! is an `Int<IBZ_LIMBS>` and arithmetic wraps on overflow. `IBZ_LIMBS` is
//! sized per security level from an empirical audit (see PORTING.md), and
//! `ibz_mul`/`ibz_pow` carry `debug_assert!` overflow guards. SQIsign signing
//! is variable-time by algorithm design, so the constant-time properties of
//! `crypto-bigint` are not the goal here; this backend exists for users who
//! want a small, audited dependency surface.

use crate::mp::Digit;
use core::cmp::Ordering;
use crypto_bigint::modular::{FixedMontyForm, FixedMontyParams};
use crypto_bigint::{Gcd, Int, JacobiSymbol, NonZero, Odd, Uint};

/// Per-level fixed precision, in 64-bit limbs.
///
/// Derived empirically by instrumenting the malachite backend across all 100
/// KAT seeds (see `examples/_ibz_audit.rs`). The high-water mark is hit inside
/// `ibz_mat_4xn_hnf_mod_core` (HNF intermediate `c = u·a_k + v·a_j` before
/// modular reduction). Observed maxima:
///
/// | level | log₂ p | max bits seen | IBZ_LIMBS | bits | margin |
/// |-------|--------|---------------|-----------|------|--------|
/// | lvl1  | 251    | 12 746        | 256       | 16 384 | 1.29× |
/// | lvl3  | 381    | 19 182        | 384       | 24 576 | 1.28× |
/// | lvl5  | 505    | 25 532        | 512       | 32 768 | 1.28× |
///
/// One limb of headroom is reserved for the sign bit / add carry.
#[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
pub const IBZ_LIMBS: usize = 256;
#[cfg(all(feature = "lvl3", not(feature = "lvl5")))]
pub const IBZ_LIMBS: usize = 384;
#[cfg(feature = "lvl5")]
pub const IBZ_LIMBS: usize = 512;

pub type Ibz = Int<IBZ_LIMBS>;
type Ubz = Uint<IBZ_LIMBS>;

#[inline]
fn neg(x: &Ibz) -> bool {
    bool::from(x.is_negative())
}

#[inline]
fn from_abs_sign(mag: Ubz, negative: bool) -> Ibz {
    Int::new_from_abs_sign(mag, crypto_bigint::Choice::from(negative as u8))
        .expect("IBZ_LIMBS overflow: magnitude does not fit in Int")
}

/// Highest bit position that may carry value (everything above must be sign extension).
#[inline]
fn check_headroom(x: &Ibz) {
    debug_assert!(
        x.abs().bits_vartime() < (IBZ_LIMBS as u32 * 64) - 1,
        "Ibz overflow risk: {} bits used of {}",
        x.abs().bits_vartime(),
        IBZ_LIMBS * 64
    );
}

// ---------------------------------------------------------------------------
// Construction / assignment

#[inline]
pub fn ibz_new() -> Ibz {
    Int::ZERO
}
#[inline]
pub fn ibz_from_i64(v: i64) -> Ibz {
    Int::from_i64(v)
}
#[inline]
pub fn ibz_set(i: &mut Ibz, x: i32) {
    *i = Int::from_i64(i64::from(x));
}
#[inline]
pub fn ibz_copy(target: &mut Ibz, value: &Ibz) {
    *target = *value;
}
#[inline]
pub fn ibz_swap(a: &mut Ibz, b: &mut Ibz) {
    core::mem::swap(a, b);
}

// ---------------------------------------------------------------------------
// Basic arithmetic

#[inline]
pub fn ibz_add(sum: &mut Ibz, a: &Ibz, b: &Ibz) {
    *sum = a.wrapping_add(b);
    check_headroom(sum);
}
#[inline]
pub fn ibz_sub(diff: &mut Ibz, a: &Ibz, b: &Ibz) {
    *diff = a.wrapping_sub(b);
    check_headroom(diff);
}
#[inline]
pub fn ibz_mul(prod: &mut Ibz, a: &Ibz, b: &Ibz) {
    #[cfg(debug_assertions)]
    {
        let need = a.abs().bits_vartime() + b.abs().bits_vartime();
        debug_assert!(
            need < (IBZ_LIMBS as u32 * 64) - 1,
            "ibz_mul overflow: {}+{} bits",
            a.abs().bits_vartime(),
            b.abs().bits_vartime()
        );
    }
    *prod = a.wrapping_mul(b);
}
#[inline]
pub fn ibz_neg(out: &mut Ibz, a: &Ibz) {
    *out = a.wrapping_neg();
}
#[inline]
pub fn ibz_abs(abs: &mut Ibz, a: &Ibz) {
    *abs = from_abs_sign(a.abs(), false);
}
#[inline]
pub fn ibz_is_negative(a: &Ibz) -> bool {
    neg(a)
}
#[inline]
pub fn ibz_is_positive(a: &Ibz) -> bool {
    !neg(a) && !bool::from(a.is_zero())
}

/// Truncating division (matches `mpz_tdiv_qr`).
pub fn ibz_div(quotient: &mut Ibz, remainder: &mut Ibz, a: &Ibz, b: &Ibz) {
    let nz = NonZero::new(*b).expect("ibz_div: divisor is zero");
    let (q, r) = a.checked_div_rem_vartime(&nz);
    *quotient = q.expect("ibz_div: MIN / -1");
    *remainder = r;
}

/// Truncating right-shift of `|a|` with sign preserved (matches `mpz_tdiv_q_2exp`).
pub fn ibz_div_2exp(quotient: &mut Ibz, a: &Ibz, exp: u32) {
    let mag = a.abs().shr_vartime(exp);
    *quotient = from_abs_sign(mag, neg(a));
}

/// Floor division (matches `mpz_fdiv_qr`).
///
/// `crypto-bigint` 0.7's `checked_div_rem_floor` returns a remainder with the
/// wrong sign when `n` and `d` have opposing signs (the identity `q·d + r = n`
/// does not hold), so we derive floor from truncating div instead.
pub fn ibz_div_floor(q: &mut Ibz, r: &mut Ibz, n: &Ibz, d: &Ibz) {
    let nz = NonZero::new(*d).expect("ibz_div_floor: divisor is zero");
    let (qt, rt) = n.checked_div_rem_vartime(&nz);
    let mut qt = qt.expect("ibz_div_floor: MIN / -1");
    let mut rt = rt;
    if !bool::from(rt.is_zero()) && (neg(n) != neg(d)) {
        qt = qt.wrapping_sub(&Ibz::ONE);
        rt = rt.wrapping_add(d);
    }
    *q = qt;
    *r = rt;
}

/// `r = a mod |b|`, always non-negative (matches `mpz_mod`).
pub fn ibz_mod(r: &mut Ibz, a: &Ibz, b: &Ibz) {
    let bm = b.abs();
    let nz = NonZero::new(bm).expect("ibz_mod: modulus is zero");
    let rr = a.abs().rem_vartime(&nz);
    *r = if neg(a) && !bool::from(rr.is_zero()) {
        from_abs_sign(bm.wrapping_sub(&rr), false)
    } else {
        from_abs_sign(rr, false)
    };
}

/// Floor remainder by an unsigned long (matches `mpz_fdiv_ui`).
pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    debug_assert!(d != 0);
    let dd = Ubz::from_u64(d);
    let nz = NonZero::new(dd).unwrap();
    let r = n.abs().rem_vartime(&nz);
    let r0 = r.as_words()[0];
    if neg(n) && r0 != 0 {
        d - r0
    } else {
        r0
    }
}

#[inline]
pub fn ibz_divides(a: &Ibz, b: &Ibz) -> i32 {
    if bool::from(b.is_zero()) {
        return bool::from(a.is_zero()) as i32;
    }
    let nz = NonZero::new(*b).unwrap();
    let (_, r) = a.checked_div_rem_vartime(&nz);
    bool::from(r.is_zero()) as i32
}

/// Non-modular `out = x^e`. The only callers in the signing path raise
/// `2` to a small power (Tonelli–Shanks loop counter, DRBG range bounds),
/// which we special-case as a shift; everything else is repeated squaring.
pub fn ibz_pow(out: &mut Ibz, x: &Ibz, e: u32) {
    if x.abs() == Ubz::from_u64(2) {
        debug_assert!(e < (IBZ_LIMBS as u32 * 64) - 1);
        let mag = Ubz::ONE.shl_vartime(e);
        *out = from_abs_sign(mag, neg(x) && e % 2 == 1);
        return;
    }
    let mut acc = Ibz::ONE;
    let mut base = *x;
    let mut ee = e;
    while ee > 0 {
        if ee & 1 == 1 {
            let a = acc;
            ibz_mul(&mut acc, &a, &base);
        }
        ee >>= 1;
        if ee > 0 {
            let b = base;
            ibz_mul(&mut base, &b, &b);
        }
    }
    *out = acc;
}

/// `out = x^e mod m` (matches `mpz_powm` for `e ≥ 0`). The only modular-exp
/// callers (Tonelli–Shanks, BPSW Miller–Rabin) always pass an odd modulus,
/// so we use `FixedMontyForm`; even moduli are handled by a slow fallback.
pub fn ibz_pow_mod(out: &mut Ibz, x: &Ibz, e: &Ibz, m: &Ibz) {
    debug_assert!(!neg(e), "negative exponents not used");
    let mm = m.abs();
    if mm == Ubz::ONE {
        *out = Ibz::ZERO;
        return;
    }
    let exp = e.abs();
    let exp_bits = exp.bits_vartime().max(1);
    let nz = NonZero::new(mm).expect("ibz_pow_mod: modulus is zero");
    let base = {
        // Reduce x into [0, m).
        let r = x.abs().rem_vartime(&nz);
        if neg(x) && !bool::from(r.is_zero()) {
            mm.wrapping_sub(&r)
        } else {
            r
        }
    };
    if let Some(odd) = Option::<Odd<Ubz>>::from(Odd::new(mm)) {
        let params = FixedMontyParams::new_vartime(odd);
        let mf = FixedMontyForm::new(&base, &params);
        let r = mf.pow_bounded_exp(&exp, exp_bits);
        *out = from_abs_sign(r.retrieve(), false);
    } else {
        // Even modulus: square-and-multiply with explicit reduction. Never on
        // the KAT path, so the cost is irrelevant.
        let mut acc = Ubz::ONE.rem_vartime(&nz);
        let mut b = base;
        for i in 0..exp_bits {
            if exp.bit_vartime(i) {
                acc = acc.wrapping_mul(&b).rem_vartime(&nz);
            }
            b = b.wrapping_mul(&b).rem_vartime(&nz);
        }
        *out = from_abs_sign(acc, false);
    }
}

/// Position of the lowest set bit (`mpz_scan1(x, 0)`).
#[inline]
pub fn ibz_two_adic(x: &Ibz) -> i32 {
    if bool::from(x.is_zero()) {
        -1
    } else {
        x.abs().trailing_zeros_vartime() as i32
    }
}

#[inline]
pub fn ibz_significant_bits(a: &Ibz) -> u64 {
    u64::from(a.abs().bits_vartime())
}
#[inline]
pub fn ibz_get_bit(a: &Ibz, i: u64) -> bool {
    // Match GMP/malachite: negatives are infinite-precision two's-complement,
    // i.e. bit_i(a) = bit_i(!( |a| − 1 )) for a < 0.
    if neg(a) {
        let m = a.abs().wrapping_sub(&Ubz::ONE);
        !m.bit_vartime(i as u32)
    } else {
        a.abs().bit_vartime(i as u32)
    }
}

#[inline]
pub fn ibz_cmp(a: &Ibz, b: &Ibz) -> i32 {
    match a.cmp_vartime(b) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}
#[inline]
pub fn ibz_is_zero(x: &Ibz) -> i32 {
    bool::from(x.is_zero()) as i32
}
#[inline]
pub fn ibz_is_one(x: &Ibz) -> i32 {
    (x.cmp_vartime(&Ibz::ONE) == Ordering::Equal) as i32
}
#[inline]
pub fn ibz_cmp_int32(x: &Ibz, y: i32) -> i32 {
    ibz_cmp(x, &Int::from_i64(i64::from(y)))
}
#[inline]
pub fn ibz_is_even(x: &Ibz) -> i32 {
    (!x.abs().bit_vartime(0)) as i32
}
#[inline]
pub fn ibz_is_odd(x: &Ibz) -> i32 {
    x.abs().bit_vartime(0) as i32
}

pub fn ibz_convert_to_str(i: &Ibz, base: i32) -> Option<String> {
    if base != 10 && base != 16 {
        return None;
    }
    let mag = i.abs();
    let s = if base == 16 {
        // crypto-bigint formats with leading zeros; strip them.
        let raw = format!("{:x}", mag);
        let trimmed = raw.trim_start_matches('0');
        if trimmed.is_empty() {
            "0".to_string()
        } else {
            trimmed.to_string()
        }
    } else {
        // Decimal via repeated div-by-1e18.
        if bool::from(mag.is_zero()) {
            "0".to_string()
        } else {
            let mut n = mag;
            let ten18 = Ubz::from_u64(1_000_000_000_000_000_000);
            let nz = NonZero::new(ten18).unwrap();
            let mut chunks: Vec<u64> = Vec::new();
            while !bool::from(n.is_zero()) {
                let (q, r) = n.div_rem_vartime(&nz);
                chunks.push(r.as_words()[0]);
                n = q;
            }
            let mut s = format!("{}", chunks.pop().unwrap());
            while let Some(c) = chunks.pop() {
                s.push_str(&format!("{:018}", c));
            }
            s
        }
    };
    Some(if neg(i) && s != "0" {
        format!("-{s}")
    } else {
        s
    })
}

#[inline]
pub fn ibz_print(num: &Ibz, base: i32) {
    print!("{}", ibz_convert_to_str(num, base).unwrap());
}

pub fn ibz_set_from_str(i: &mut Ibz, s: &str, base: i32) -> i32 {
    let (negative, body) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s),
    };
    if body.is_empty() {
        return 0;
    }
    let mut acc = Ubz::ZERO;
    let bb = Ubz::from_u64(base as u64);
    for c in body.chars() {
        let d = match c.to_digit(base as u32) {
            Some(d) => d as u64,
            None => return 0,
        };
        acc = acc.wrapping_mul(&bb).wrapping_add(&Ubz::from_u64(d));
    }
    *i = from_abs_sign(acc, negative);
    1
}

#[inline]
pub fn ibz_bitsize(a: &Ibz) -> i32 {
    a.abs().bits_vartime() as i32
}
pub fn ibz_size_in_base(a: &Ibz, base: i32) -> i32 {
    if base == 2 {
        a.abs().bits_vartime().max(1) as i32
    } else {
        ibz_convert_to_str(a, base)
            .unwrap()
            .trim_start_matches('-')
            .len()
            .max(1) as i32
    }
}

// ---------------------------------------------------------------------------
// Digit / byte conversion

/// Read little-endian u64 limbs into `target` (matches `mpz_import`).
pub fn ibz_copy_digits(target: &mut Ibz, dig: &[Digit]) {
    debug_assert!(dig.len() <= IBZ_LIMBS);
    let mut words = [0u64; IBZ_LIMBS];
    words[..dig.len()].copy_from_slice(dig);
    *target = from_abs_sign(Ubz::from_words(words), false);
}

/// Overwrite the value with zeros. `Int<L>` is plain stack data.
pub fn ibz_secure_clear(x: &mut Ibz) {
    let words: &mut [u64; IBZ_LIMBS] = x.as_mut();
    zeroize::Zeroize::zeroize(words.as_mut_slice());
}

/// Write `ibz` to `target` as little-endian u64 limbs, zero-padding to len.
pub fn ibz_to_digits(target: &mut [Digit], ibz: &Ibz) {
    debug_assert!(!neg(ibz));
    let words = ibz.abs().to_words();
    debug_assert!(
        words[target.len()..].iter().all(|&w| w == 0),
        "target too small for ibz"
    );
    target.copy_from_slice(&words[..target.len()]);
}

/// Exact `mpz_get_si` → C `ibz_get` semantics. See the malachite backend for
/// the derivation; we extract the low limb of `|x|`, mask to 63 bits, apply
/// the sign, then fold to 32 bits as the C wrapper does.
pub fn ibz_get(i: &Ibz) -> i32 {
    let mag = i.abs();
    let low_limb = mag.as_words()[0];
    let abs63 = (low_limb & (i64::MAX as u64)) as i64;
    let t: i64 = if neg(i) { -abs63 } else { abs63 };
    let sign_bit = ((t >> 32) as i32) & i32::MIN;
    let low = (t as i32) & i32::MAX;
    sign_bit | low
}

/// `mpz_get_d_2exp` semantics: `(d, exp)` with `0.5 ≤ |d| < 1` and `d`
/// truncated (toward zero) to 53 bits.
pub fn ibz_get_d_2exp(z: &Ibz) -> (f64, i64) {
    if bool::from(z.is_zero()) {
        return (0.0, 0);
    }
    let mag = z.abs();
    let bits = mag.bits_vartime() as i64;
    let top = if bits > 53 {
        mag.shr_vartime((bits - 53) as u32)
    } else {
        mag
    };
    let top_u = top.as_words()[0];
    let m = top_u as f64;
    let scale = if bits > 53 { 53 } else { bits };
    let d = m / (1u64 << scale) as f64;
    (if neg(z) { -d } else { d }, bits)
}

pub fn ibz_set_from_mantissa_shift(out: &mut Ibz, int_m: i64, shift: i64) {
    let mag = Ubz::from_u64(int_m.unsigned_abs());
    let mag = if shift >= 0 {
        mag.shl_vartime(shift as u32)
    } else {
        mag.shr_vartime((-shift) as u32)
    };
    *out = from_abs_sign(mag, int_m < 0);
}

// ---------------------------------------------------------------------------
// Number theory

pub fn ibz_gcd(gcd: &mut Ibz, a: &Ibz, b: &Ibz) {
    *gcd = from_abs_sign(a.gcd(b), false);
}

/// Extended GCD with the `mpz_gcdext` cofactor convention (`|u| ≤ |b/(2g)|`,
/// `|v| ≤ |a/(2g)|`, plus sign rules). `crypto-bigint`'s binary xgcd does not
/// guarantee this convention, so we normalise afterwards. The HNF caller
/// already re-normalises in `ibz_xgcd_with_u_not_0`, but other callers
/// (e.g. tests) expect the GMP convention directly.
pub fn ibz_xgcd(gcd: &mut Ibz, u: &mut Ibz, v: &mut Ibz, a: &Ibz, b: &Ibz) {
    let out = a.xgcd(b);
    *gcd = from_abs_sign(out.gcd, false);
    // GMP edge cases.
    let az = bool::from(a.is_zero());
    let bz = bool::from(b.is_zero());
    if az && bz {
        *u = Ibz::ZERO;
        *v = Ibz::ZERO;
        return;
    }
    if bz {
        *u = if neg(a) { Ibz::MINUS_ONE } else { Ibz::ONE };
        *v = Ibz::ZERO;
        return;
    }
    if az || a.abs() == b.abs() {
        *u = Ibz::ZERO;
        *v = if neg(b) { Ibz::MINUS_ONE } else { Ibz::ONE };
        return;
    }
    // General case: GMP picks the (u,v) with |u| ≤ |b|/(2g). All Bézout pairs
    // differ by multiples of (b/g, −a/g), so reduce crypto-bigint's u to the
    // centred residue mod |b/g|, then derive v from the identity.
    let bog = out.rhs_on_gcd.abs(); // |b/g|, ≥ 2 here.
    let bog_i = from_abs_sign(bog, false);
    let mut q_ = Ibz::ZERO;
    let mut r = Ibz::ZERO;
    ibz_div_floor(&mut q_, &mut r, &out.x, &bog_i);
    // r ∈ [0, |b/g|); centre to (−|b|/2g, |b|/2g].
    let half = from_abs_sign(bog.shr_vartime(1), false);
    if r.cmp_vartime(&half) == Ordering::Greater {
        r = r.wrapping_sub(&bog_i);
    }
    // GMP edge case: when |b| = 2g (so |b/g| = 2, half = 1, r ∈ {0,1}),
    // mpz_gcdext picks u = sgn(a). Our centring can only yield 0 or 1 here.
    if bog == Ubz::from_u64(2) {
        r = if neg(a) { Ibz::MINUS_ONE } else { Ibz::ONE };
    }
    *u = r;
    // v = (g − u·a) / b  (exact since u·a + v·b = g).
    let num = gcd.wrapping_sub(&u.wrapping_mul(a));
    let nz_b = NonZero::new(*b).unwrap();
    let (vv, rr) = num.checked_div_rem_vartime(&nz_b);
    debug_assert!(bool::from(rr.is_zero()));
    *v = vv.expect("xgcd norm");
}

pub fn ibz_invmod(inv: &mut Ibz, a: &Ibz, m: &Ibz) -> i32 {
    let mm = m.abs();
    if bool::from(mm.is_zero()) {
        return 0;
    }
    let mut ar = Ibz::ZERO;
    ibz_mod(&mut ar, a, m);
    let out = ar.xgcd(&from_abs_sign(mm, false));
    if out.gcd != Ubz::ONE {
        return 0;
    }
    // x·a' + y·m ≡ 1 → x ≡ a⁻¹ (mod m). Reduce into [0, m).
    ibz_mod(inv, &out.x, m);
    1
}

#[inline]
fn jacobi_to_i32(j: JacobiSymbol) -> i32 {
    match j {
        JacobiSymbol::MinusOne => -1,
        JacobiSymbol::Zero => 0,
        JacobiSymbol::One => 1,
    }
}

#[inline]
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    let odd = Odd::new(p.abs()).expect("ibz_legendre: even modulus");
    jacobi_to_i32(a.jacobi_symbol_vartime(&odd))
}

pub fn ibz_sqrt(sqrt: &mut Ibz, a: &Ibz) -> i32 {
    if neg(a) {
        return 0;
    }
    let r = a.abs().floor_sqrt_vartime();
    if r.wrapping_mul(&r) == a.abs() {
        *sqrt = from_abs_sign(r, false);
        1
    } else {
        0
    }
}

#[inline]
pub fn ibz_sqrt_floor(sqrt: &mut Ibz, a: &Ibz) {
    // rug/malachite panic on negative input; align rather than silently
    // returning sqrt(|a|).
    debug_assert!(!neg(a), "ibz_sqrt_floor of negative");
    *sqrt = from_abs_sign(a.abs().floor_sqrt_vartime(), false);
}

// ---------------------------------------------------------------------------
// Primality: Baillie–PSW (deterministic; matches GMP's first phase).
// Expressed entirely via crypto-bigint primitives.

const SMALL_PRIMES: &[u64] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

fn miller_rabin_round(
    n: &Ubz,
    n_minus_1: &Ubz,
    d: &Ubz,
    s: u32,
    a: &Ubz,
    params: &FixedMontyParams<IBZ_LIMBS>,
) -> bool {
    let af = FixedMontyForm::new(a, params);
    let dbits = d.bits_vartime().max(1);
    let mut x = af.pow_bounded_exp(d, dbits).retrieve();
    if x == Ubz::ONE || x == *n_minus_1 {
        return true;
    }
    let nz = NonZero::new(*n).unwrap();
    for _ in 1..s {
        x = x.wrapping_mul(&x).rem_vartime(&nz);
        if x == *n_minus_1 {
            return true;
        }
        if x == Ubz::ONE {
            return false;
        }
    }
    false
}

fn strong_lucas_prp(n: &Ubz) -> bool {
    let n_int = from_abs_sign(*n, false);
    let n_odd = Odd::new(*n).expect("n is odd");
    // Selfridge: D = 5, -7, 9, -11, ... with Jacobi(D/n) = -1.
    let mut d: i64 = 5;
    loop {
        let dd: Ibz = Int::from_i64(d);
        let j = jacobi_to_i32(dd.jacobi_symbol_vartime(&n_odd));
        if j == -1 {
            break;
        }
        if j == 0 {
            return n_int.abs() == Int::<IBZ_LIMBS>::from_i64(d).abs();
        }
        if d.abs() > 1_000_000 {
            return false;
        }
        d = if d > 0 { -(d + 2) } else { -(d - 2) };
    }
    let q = (1 - d) / 4;
    let q_int = Ibz::from_i64(q);

    let n_plus_1 = n.wrapping_add(&Ubz::ONE);
    let s = n_plus_1.trailing_zeros_vartime();
    let k = n_plus_1.shr_vartime(s);
    let bits = k.bits_vartime();

    // Work mod n with values in [0, n) stored as Ibz (always non-negative).
    let modn = |x: &Ibz| -> Ibz {
        let mut r = Ibz::ZERO;
        ibz_mod(&mut r, x, &n_int);
        r
    };
    let mut two_inv = Ibz::ZERO;
    let _ = ibz_invmod(&mut two_inv, &Ibz::from_i64(2), &n_int);
    let d_int = Ibz::from_i64(d);

    let mut u = Ibz::ONE;
    let mut v = Ibz::ONE;
    let mut qk = modn(&q_int);

    for i in (0..bits.saturating_sub(1)).rev() {
        let u2 = modn(&u.wrapping_mul(&v));
        let v2 = modn(&v.wrapping_mul(&v).wrapping_sub(&qk.wrapping_add(&qk)));
        u = u2;
        v = v2;
        qk = modn(&qk.wrapping_mul(&qk));
        if k.bit_vartime(i) {
            let u1 = modn(&u.wrapping_add(&v).wrapping_mul(&two_inv));
            let v1 = modn(
                &d_int
                    .wrapping_mul(&u)
                    .wrapping_add(&v)
                    .wrapping_mul(&two_inv),
            );
            u = u1;
            v = v1;
            qk = modn(&qk.wrapping_mul(&q_int));
        }
    }

    if bool::from(u.is_zero()) || bool::from(v.is_zero()) {
        return true;
    }
    for _ in 1..s {
        v = modn(&v.wrapping_mul(&v).wrapping_sub(&qk.wrapping_add(&qk)));
        qk = modn(&qk.wrapping_mul(&qk));
        if bool::from(v.is_zero()) {
            return true;
        }
    }
    false
}

pub fn ibz_probab_prime(n: &Ibz, _reps: i32) -> i32 {
    // GMP tests |n|.
    let nn = n.abs();
    if nn.cmp_vartime(&Ubz::from_u64(2)) == Ordering::Less {
        return 0;
    }
    for &p in SMALL_PRIMES {
        let pp = Ubz::from_u64(p);
        if nn == pp {
            return 2;
        }
        let nz = NonZero::new(pp).unwrap();
        if bool::from(nn.rem_vartime(&nz).is_zero()) {
            return 0;
        }
    }
    // Perfect-square reject (breaks Selfridge D-search otherwise).
    let r = nn.floor_sqrt_vartime();
    if r.wrapping_mul(&r) == nn {
        return 0;
    }
    let n_minus_1 = nn.wrapping_sub(&Ubz::ONE);
    let s = n_minus_1.trailing_zeros_vartime();
    let d = n_minus_1.shr_vartime(s);
    let params = FixedMontyParams::new_vartime(Odd::new(nn).unwrap());
    if !miller_rabin_round(&nn, &n_minus_1, &d, s, &Ubz::from_u64(2), &params) {
        return 0;
    }
    if !strong_lucas_prp(&nn) {
        return 0;
    }
    1
}
