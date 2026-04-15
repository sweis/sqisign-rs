//! Malachite backend for `Ibz` (the default, pure-Rust bigint).
//!
//! Exposes the same `ibz_*` surface as the rug backend. The only nontrivial
//! piece is `ibz_probab_prime`, which malachite lacks for arbitrary-precision
//! integers; a deterministic Baillie–PSW test is implemented here so that
//! prime/composite classification matches GMP's `mpz_probab_prime_p` (which
//! also begins with BPSW) on every input the KAT path encounters.

use crate::mp::Digit;
use core::cmp::Ordering;
use malachite_base::num::arithmetic::traits::{
    Abs, DivMod, DivRem, DivisibleBy, ExtendedGcd, FloorSqrt, JacobiSymbol, LegendreSymbol, Mod,
    ModPow, Parity, Pow, Sign, UnsignedAbs,
};
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_base::num::conversion::traits::{ExactFrom, FromStringBase, ToStringBase};
use malachite_base::num::logic::traits::{BitAccess, SignificantBits};
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

pub type Ibz = Integer;

#[cfg(ibz_audit)]
mod audit {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering as AO};
    static MAX_BITS: AtomicU64 = AtomicU64::new(0);
    pub fn track(x: &Ibz) {
        MAX_BITS.fetch_max(x.significant_bits(), AO::Relaxed);
    }
    pub fn report() {
        eprintln!("[ibz_audit] max_bits = {}", MAX_BITS.load(AO::Relaxed));
    }
}
#[cfg(ibz_audit)]
pub use audit::report as ibz_audit_report;
#[cfg(ibz_audit)]
#[inline]
fn track(x: &Ibz) {
    audit::track(x);
}
#[cfg(not(ibz_audit))]
#[inline]
fn track(_: &Ibz) {}
#[inline]
fn track_mul(x: &Ibz) {
    track(x);
}

// ---------------------------------------------------------------------------
// Construction / assignment

#[inline]
pub fn ibz_new() -> Ibz {
    Integer::ZERO
}
#[inline]
pub fn ibz_from_i64(v: i64) -> Ibz {
    Integer::from(v)
}
#[inline]
pub fn ibz_set(i: &mut Ibz, x: i32) {
    *i = Integer::from(x);
}
#[inline]
pub fn ibz_copy(target: &mut Ibz, value: &Ibz) {
    target.clone_from(value);
}
#[inline]
pub fn ibz_swap(a: &mut Ibz, b: &mut Ibz) {
    core::mem::swap(a, b);
}

// ---------------------------------------------------------------------------
// Basic arithmetic

#[inline]
pub fn ibz_add(sum: &mut Ibz, a: &Ibz, b: &Ibz) {
    *sum = a + b;
    track(sum);
}
#[inline]
pub fn ibz_sub(diff: &mut Ibz, a: &Ibz, b: &Ibz) {
    *diff = a - b;
    track(diff);
}
#[inline]
pub fn ibz_mul(prod: &mut Ibz, a: &Ibz, b: &Ibz) {
    *prod = a * b;
    track_mul(prod);
}
#[inline]
pub fn ibz_neg(neg: &mut Ibz, a: &Ibz) {
    *neg = -a;
}
#[inline]
pub fn ibz_abs(abs: &mut Ibz, a: &Ibz) {
    *abs = a.abs();
}
#[inline]
pub fn ibz_is_negative(a: &Ibz) -> bool {
    a.sign() == Ordering::Less
}
#[inline]
pub fn ibz_is_positive(a: &Ibz) -> bool {
    a.sign() == Ordering::Greater
}

/// Truncating division (matches `mpz_tdiv_qr`).
pub fn ibz_div(quotient: &mut Ibz, remainder: &mut Ibz, a: &Ibz, b: &Ibz) {
    let (q, r) = a.div_rem(b);
    *quotient = q;
    *remainder = r;
}

/// Truncating right-shift of |a| with sign preserved (matches `mpz_tdiv_q_2exp`).
pub fn ibz_div_2exp(quotient: &mut Ibz, a: &Ibz, exp: u32) {
    let mag: Natural = a.unsigned_abs_ref() >> u64::from(exp);
    *quotient = if ibz_is_negative(a) {
        -Integer::from(mag)
    } else {
        Integer::from(mag)
    };
}

/// Floor division (matches `mpz_fdiv_qr`).
pub fn ibz_div_floor(q: &mut Ibz, r: &mut Ibz, n: &Ibz, d: &Ibz) {
    let (qq, rr) = n.div_mod(d);
    *q = qq;
    *r = rr;
}

/// `r = a mod |b|`, always non-negative (matches `mpz_mod`).
pub fn ibz_mod(r: &mut Ibz, a: &Ibz, b: &Ibz) {
    let m: Natural = b.unsigned_abs_ref().clone();
    let rr = a.mod_op(&Integer::from(m));
    *r = rr;
}

/// Floor remainder by an unsigned long (matches `mpz_fdiv_ui`).
pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    debug_assert!(d != 0);
    let m = Integer::from(d);
    let r = n.mod_op(&m);
    u64::exact_from(&r)
}

#[inline]
pub fn ibz_divides(a: &Ibz, b: &Ibz) -> bool {
    a.divisible_by(b)
}

#[inline]
pub fn ibz_pow(out: &mut Ibz, x: &Ibz, e: u32) {
    *out = x.pow(u64::from(e));
    track(out);
}

pub fn ibz_pow_mod(out: &mut Ibz, x: &Ibz, e: &Ibz, m: &Ibz) {
    debug_assert!(!ibz_is_negative(e), "negative exponents not used");
    let mm: &Natural = m.unsigned_abs_ref();
    let xm: Integer = x.mod_op(&Integer::from(mm.clone()));
    let base: Natural = xm.unsigned_abs();
    let exp: &Natural = e.unsigned_abs_ref();
    *out = Integer::from(base.mod_pow(exp, mm));
}

/// Position of the lowest set bit (`mpz_scan1(x, 0)`).
#[inline]
pub fn ibz_two_adic(x: &Ibz) -> i32 {
    x.trailing_zeros().map_or(-1, |b| b as i32)
}

#[inline]
pub fn ibz_significant_bits(a: &Ibz) -> u64 {
    a.significant_bits()
}

#[inline]
pub fn ibz_get_bit(a: &Ibz, i: u64) -> bool {
    a.get_bit(i)
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
pub fn ibz_is_zero(x: &Ibz) -> bool {
    *x == Integer::ZERO
}
#[inline]
pub fn ibz_is_one(x: &Ibz) -> bool {
    *x == Integer::ONE
}
#[inline]
pub fn ibz_cmp_int32(x: &Ibz, y: i32) -> i32 {
    match x.partial_cmp(&i64::from(y)).unwrap() {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}
#[inline]
pub fn ibz_is_even(x: &Ibz) -> bool {
    x.even()
}
#[inline]
pub fn ibz_is_odd(x: &Ibz) -> bool {
    x.odd()
}

pub fn ibz_convert_to_str(i: &Ibz, base: i32) -> Option<String> {
    if base != 10 && base != 16 {
        return None;
    }
    Some(i.to_string_base(base as u8))
}

#[inline]
pub fn ibz_print(num: &Ibz, base: i32) {
    print!("{}", num.to_string_base(base as u8));
}

pub fn ibz_set_from_str(i: &mut Ibz, s: &str, base: i32) -> i32 {
    let (neg, body) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s),
    };
    match Natural::from_string_base(base as u8, body) {
        Some(n) => {
            *i = if neg {
                -Integer::from(n)
            } else {
                Integer::from(n)
            };
            1
        }
        None => 0,
    }
}

#[inline]
pub fn ibz_bitsize(a: &Ibz) -> i32 {
    a.significant_bits() as i32
}

pub fn ibz_size_in_base(a: &Ibz, base: i32) -> i32 {
    if base == 2 {
        a.significant_bits().max(1) as i32
    } else {
        a.to_string_base(base as u8)
            .trim_start_matches('-')
            .len()
            .max(1) as i32
    }
}

// ---------------------------------------------------------------------------
// Digit / byte conversion

/// Read little-endian u64 limbs into `target` (matches `mpz_import`).
pub fn ibz_copy_digits(target: &mut Ibz, dig: &[Digit]) {
    *target = Integer::from(Natural::from_limbs_asc(dig));
    track(target);
}

/// Best-effort secure clear. Malachite does not expose mutable access to its
/// internal limb `Vec`, so this overwrites the value bit-by-bit (forcing the
/// existing buffer to all-zero bits) before reassigning. Intermediate buffers
/// from prior reallocations are not covered.
pub fn ibz_secure_clear(x: &mut Ibz) {
    use malachite_base::num::logic::traits::BitAccess;
    let bits = ibz_significant_bits(x);
    for i in 0..bits {
        x.clear_bit(i);
    }
    *x = Integer::from(0);
}

/// Write `ibz` to `target` as little-endian u64 limbs, zero-padding to len.
pub fn ibz_to_digits(target: &mut [Digit], ibz: &Ibz) {
    debug_assert!(!ibz_is_negative(ibz));
    for d in target.iter_mut() {
        *d = 0;
    }
    let limbs = ibz.unsigned_abs_ref().to_limbs_asc();
    debug_assert!(limbs.len() <= target.len(), "target too small for ibz");
    target[..limbs.len()].copy_from_slice(&limbs);
}

/// Exact `mpz_get_si` → C `ibz_get` semantics. GMP's `mpz_get_si` returns
/// `sign(x) · (|x|.limb[0] & i64::MAX)` as an `i64`; the C wrapper then folds
/// that into 32 bits as `(t>>32 & i32::MIN) | (t & i32::MAX)`.
pub fn ibz_get(i: &Ibz) -> i32 {
    let mag = i.unsigned_abs_ref();
    let low_limb = if mag.limb_count() == 0 {
        0u64
    } else {
        mag.to_limbs_asc()[0]
    };
    let abs63 = (low_limb & (i64::MAX as u64)) as i64;
    let t: i64 = if ibz_is_negative(i) { -abs63 } else { abs63 };
    let sign_bit = ((t >> 32) as i32) & i32::MIN;
    let low = (t as i32) & i32::MAX;
    sign_bit | low
}

/// `mpz_get_d_2exp` semantics: `(d, exp)` with `0.5 <= |d| < 1` and `d`
/// truncated (toward zero) to 53 bits, such that `z = trunc(d · 2^exp)`.
pub fn ibz_get_d_2exp(z: &Ibz) -> (f64, i64) {
    if *z == Integer::ZERO {
        return (0.0, 0);
    }
    let bits = z.significant_bits() as i64;
    let mag = z.unsigned_abs_ref();
    let top: Natural = if bits > 53 {
        mag >> ((bits - 53) as u64)
    } else {
        mag.clone()
    };
    let top_u = u64::exact_from(&top);
    // top_u has ≤53 bits → conversion to f64 is exact.
    let m = top_u as f64;
    let scale = if bits > 53 { 53 } else { bits };
    let d = m / (1u64 << scale) as f64; // now in [0.5, 1)
    (if ibz_is_negative(z) { -d } else { d }, bits)
}

/// Set `out = z` from a 53-bit signed integer mantissa and a power-of-two
/// shift (positive = left, negative = right truncating).
pub fn ibz_set_from_mantissa_shift(out: &mut Ibz, int_m: i64, shift: i64) {
    *out = Integer::from(int_m);
    if shift >= 0 {
        *out <<= shift as u64;
        track(out);
    } else {
        // Truncating division by power of two (toward zero).
        let mag: Natural = out.unsigned_abs_ref() >> ((-shift) as u64);
        *out = if int_m < 0 {
            -Integer::from(mag)
        } else {
            Integer::from(mag)
        };
    }
}

// ---------------------------------------------------------------------------
// Number theory

pub fn ibz_gcd(gcd: &mut Ibz, a: &Ibz, b: &Ibz) {
    use malachite_base::num::arithmetic::traits::Gcd;
    let g = a.unsigned_abs_ref().gcd(b.unsigned_abs_ref());
    *gcd = Integer::from(g);
}

/// Extended GCD: `gcd = u·a + v·b` with `mpz_gcdext` cofactor convention
/// (`|u| ≤ |b/(2g)|`, `|v| ≤ |a/(2g)|`, and the documented edge cases).
pub fn ibz_xgcd(gcd: &mut Ibz, u: &mut Ibz, v: &mut Ibz, a: &Ibz, b: &Ibz) {
    // Malachite's extended_gcd already follows the GMP convention.
    let (g, su, sv) = a.clone().extended_gcd(b.clone());
    *gcd = Integer::from(g);
    *u = su;
    *v = sv;
}

pub fn ibz_invmod(inv: &mut Ibz, a: &Ibz, m: &Ibz) -> i32 {
    let mm = Integer::from(m.unsigned_abs_ref().clone());
    let aa = a.mod_op(&mm);
    let (g, x, _y) = aa.extended_gcd(mm.clone());
    if g != Natural::ONE {
        return 0;
    }
    *inv = x.mod_op(&mm);
    1
}

#[inline]
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    a.legendre_symbol(p) as i32
}

pub fn ibz_sqrt(sqrt: &mut Ibz, a: &Ibz) -> bool {
    if ibz_is_negative(a) {
        return false;
    }
    let r = a.unsigned_abs_ref().floor_sqrt();
    let exact = &r * &r == *a.unsigned_abs_ref();
    if exact {
        *sqrt = Integer::from(r);
    }
    exact
}

#[inline]
pub fn ibz_sqrt_floor(sqrt: &mut Ibz, a: &Ibz) {
    *sqrt = a.floor_sqrt();
}

// ---------------------------------------------------------------------------
// Primality: Baillie–PSW (deterministic; matches GMP's mpz_probab_prime_p
// on all inputs with no known counterexample).

const SMALL_PRIMES: &[u32] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

fn nat_mod_pow(base: &Natural, exp: &Natural, m: &Natural) -> Natural {
    (base.mod_op(m)).mod_pow(exp, m)
}

/// One strong probable-prime (Miller–Rabin) round with base `a` against `n`.
fn miller_rabin_round(n: &Natural, n_minus_1: &Natural, d: &Natural, s: u64, a: &Natural) -> bool {
    let mut x = nat_mod_pow(a, d, n);
    if x == Natural::ONE || x == *n_minus_1 {
        return true;
    }
    for _ in 1..s {
        x = (&x * &x).mod_op(n);
        if x == *n_minus_1 {
            return true;
        }
        if x == Natural::ONE {
            return false;
        }
    }
    false
}

/// Strong Lucas probable-prime test (Selfridge parameters), matching FIPS 186-4
/// Appendix C.3.3 / GMP's `mpz_stronglucas`.
fn strong_lucas_prp(n: &Natural) -> bool {
    // Selfridge: find first D in 5, -7, 9, -11, ... with Jacobi(D/n) = -1.
    let mut d: i64 = 5;
    let n_int = Integer::from(n.clone());
    loop {
        let j = Integer::from(d).jacobi_symbol(&n_int);
        if j == -1 {
            break;
        }
        if j == 0 {
            // n shares a factor with |D|; composite unless n == |D|.
            return n_int == Integer::from(d).abs();
        }
        if d.abs() > 1_000_000 {
            // Perfect squares never reach -1; caller already excluded them.
            return false;
        }
        d = if d > 0 { -(d + 2) } else { -(d - 2) };
    }
    // P=1, Q=(1-D)/4.
    let q = (1 - d) / 4;
    let q_int = Integer::from(q);

    // Factor n+1 = 2^s · k with k odd.
    let n_plus_1: Natural = n + Natural::ONE;
    let s = n_plus_1.trailing_zeros().unwrap();
    let k: Natural = &n_plus_1 >> s;

    // Lucas sequence by binary method: walk bits of k from msb-1 to 0.
    let bits = k.significant_bits();
    let modn = |x: &Integer| x.mod_op(&n_int);
    let mut u = Integer::ONE;
    let mut v = Integer::ONE;
    let mut qk = modn(&q_int);
    let two_inv = {
        // n is odd ⇒ 2 is invertible.
        let mut t = Integer::ZERO;
        let _ = ibz_invmod(&mut t, &Integer::TWO, &n_int);
        t
    };
    let d_int = Integer::from(d);

    for i in (0..bits.saturating_sub(1)).rev() {
        // Double: U_{2m}=U_m V_m,  V_{2m}=V_m^2 - 2 Q^m.
        let u2 = modn(&(&u * &v));
        let v2 = modn(&(&v * &v - (&qk + &qk)));
        u = u2;
        v = v2;
        qk = modn(&(&qk * &qk));
        if k.get_bit(i) {
            // Add-one: U_{m+1}=(P U_m + V_m)/2, V_{m+1}=(D U_m + P V_m)/2. P=1.
            let u1 = modn(&(&(&u + &v) * &two_inv));
            let v1 = modn(&(&(&d_int * &u + &v) * &two_inv));
            u = u1;
            v = v1;
            qk = modn(&(&qk * &q_int));
        }
    }

    // Strong test: U_k ≡ 0, or V_{k·2^r} ≡ 0 for some 0 ≤ r < s.
    if u == Integer::ZERO || v == Integer::ZERO {
        return true;
    }
    for _ in 1..s {
        v = modn(&(&v * &v - (&qk + &qk)));
        qk = modn(&(&qk * &qk));
        if v == Integer::ZERO {
            return true;
        }
    }
    false
}

/// `mpz_probab_prime_p`-compatible primality. Returns 2 for small certain
/// primes, 1 for probable primes, 0 for composites. GMP tests |n|, so we do too.
pub fn ibz_probab_prime(n: &Ibz, _reps: i32) -> i32 {
    let nn: &Natural = n.unsigned_abs_ref();
    if *nn < Natural::TWO {
        return 0;
    }
    // Trial division.
    for &p in SMALL_PRIMES {
        let pp = Natural::from(p);
        if *nn == pp {
            return 2;
        }
        if (nn.mod_op(&pp)) == Natural::ZERO {
            return 0;
        }
    }
    // Perfect squares are never prime and break the Lucas D-search.
    {
        let r = nn.floor_sqrt();
        if &r * &r == *nn {
            return 0;
        }
    }
    // Miller–Rabin base 2.
    let n_minus_1: Natural = nn - Natural::ONE;
    let s = n_minus_1.trailing_zeros().unwrap();
    let d: Natural = &n_minus_1 >> s;
    if !miller_rabin_round(nn, &n_minus_1, &d, s, &Natural::TWO) {
        return 0;
    }
    // Strong Lucas.
    if !strong_lucas_prp(nn) {
        return 0;
    }
    1
}
