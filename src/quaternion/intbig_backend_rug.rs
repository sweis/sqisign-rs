//! GMP backend for `Ibz`, via the `rug` crate. Selected with `feature = "gmp"`.

use crate::mp::Digit;
use core::cmp::Ordering;
use rug::integer::{IsPrime, Order};
use rug::ops::{NegAssign, Pow, RemRounding};
use rug::{Assign, Complete, Integer};

pub type Ibz = Integer;

// Construction / assignment

#[inline]
pub fn ibz_new() -> Ibz {
    Integer::new()
}
#[inline]
pub fn ibz_from_i64(v: i64) -> Ibz {
    Integer::from(v)
}
#[inline]
pub fn ibz_set(i: &mut Ibz, x: i32) {
    i.assign(x);
}
#[inline]
pub fn ibz_copy(target: &mut Ibz, value: &Ibz) {
    target.assign(value);
}
#[inline]
pub fn ibz_swap(a: &mut Ibz, b: &mut Ibz) {
    core::mem::swap(a, b);
}

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
#[inline]
pub fn ibz_is_negative(a: &Ibz) -> bool {
    a.is_negative()
}
#[inline]
pub fn ibz_is_positive(a: &Ibz) -> bool {
    a.is_positive()
}

pub fn ibz_div(quotient: &mut Ibz, remainder: &mut Ibz, a: &Ibz, b: &Ibz) {
    let (q, r) = a.div_rem_ref(b).complete();
    quotient.assign(q);
    remainder.assign(r);
}

pub fn ibz_div_2exp(quotient: &mut Ibz, a: &Ibz, exp: u32) {
    let neg = a.is_negative();
    quotient.assign(a.abs_ref());
    *quotient >>= exp;
    if neg {
        quotient.neg_assign();
    }
}

pub fn ibz_div_floor(q: &mut Ibz, r: &mut Ibz, n: &Ibz, d: &Ibz) {
    let (qq, rr) = n.div_rem_floor_ref(d).complete();
    q.assign(qq);
    r.assign(rr);
}

pub fn ibz_mod(r: &mut Ibz, a: &Ibz, b: &Ibz) {
    let abs_b = b.clone().abs();
    r.assign(a.rem_floor(&abs_b));
}

pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    debug_assert!(d != 0);
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

#[inline]
pub fn ibz_two_adic(x: &Ibz) -> i32 {
    x.find_one(0).map(|b| b as i32).unwrap_or(-1)
}
#[inline]
pub fn ibz_significant_bits(a: &Ibz) -> u64 {
    u64::from(a.significant_bits())
}
#[inline]
pub fn ibz_get_bit(a: &Ibz, i: u64) -> bool {
    a.get_bit(i as u32)
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

pub fn ibz_convert_to_str(i: &Ibz, base: i32) -> Option<String> {
    if base != 10 && base != 16 {
        return None;
    }
    Some(i.to_string_radix(base))
}
#[inline]
pub fn ibz_print(num: &Ibz, base: i32) {
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
pub fn ibz_bitsize(a: &Ibz) -> i32 {
    a.significant_bits() as i32
}
pub fn ibz_size_in_base(a: &Ibz, base: i32) -> i32 {
    if base == 2 {
        a.significant_bits().max(1) as i32
    } else {
        a.to_string_radix(base).trim_start_matches('-').len().max(1) as i32
    }
}

// Digit conversion

pub fn ibz_copy_digits(target: &mut Ibz, dig: &[Digit]) {
    target.assign_digits(dig, Order::Lsf);
}

/// Overwrite the limb buffer with zeros and set the value to 0.
#[allow(unsafe_code)]
pub fn ibz_secure_clear(x: &mut Ibz) {
    // SAFETY: rug's Integer is a transparent wrapper over mpz_t. We zero the
    // allocated limb buffer (`_mp_d[0.._mp_alloc]`) before setting size to 0.
    unsafe {
        let raw = x.as_raw_mut();
        let alloc = (*raw).alloc as usize;
        if alloc > 0 {
            let limbs = core::slice::from_raw_parts_mut((*raw).d.as_ptr(), alloc);
            zeroize::Zeroize::zeroize(limbs);
        }
        (*raw).size = 0;
    }
}
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

#[allow(unsafe_code)]
pub fn ibz_get(i: &Ibz) -> i32 {
    // SAFETY: rug's Integer is a transparent wrapper over mpz_t.
    let t: i64 = unsafe { gmp_mpfr_sys::gmp::mpz_get_si(i.as_raw()) };
    let sign_bit = ((t >> 32) as i32) & i32::MIN;
    let low = (t as i32) & i32::MAX;
    sign_bit | low
}

pub fn ibz_get_d_2exp(z: &Ibz) -> (f64, i64) {
    if z.is_zero() {
        return (0.0, 0);
    }
    let bits = z.significant_bits() as i64;
    let shift = (bits - 53).max(0);
    let top: Ibz = z.clone().abs() >> (shift as u32);
    let top_u = top.to_u64().expect("≤53 bits");
    let m = top_u as f64;
    let scale = if bits > 53 { 53 } else { bits };
    let d = m / (1u64 << scale) as f64;
    (if z.is_negative() { -d } else { d }, bits)
}

pub fn ibz_set_from_mantissa_shift(out: &mut Ibz, int_m: i64, shift: i64) {
    out.assign(int_m);
    if shift >= 0 {
        *out <<= shift as u32;
    } else {
        let neg = out.is_negative();
        if neg {
            out.neg_assign();
        }
        *out >>= (-shift) as u32;
        if neg {
            out.neg_assign();
        }
    }
}

// Number theory

pub fn ibz_gcd(gcd: &mut Ibz, a: &Ibz, b: &Ibz) {
    gcd.assign(a.gcd_ref(b));
}

/// Extended GCD: `gcd = u·a + v·b` (matches `mpz_gcdext`).
pub fn ibz_xgcd(gcd: &mut Ibz, u: &mut Ibz, v: &mut Ibz, a: &Ibz, b: &Ibz) {
    let (g, su, sv) = a.clone().extended_gcd(b.clone(), Integer::new());
    *gcd = g;
    *u = su;
    *v = sv;
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
