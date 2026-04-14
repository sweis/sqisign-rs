// SPDX-License-Identifier: Apache-2.0
//! Double-plus-exponent: an `f64` mantissa with a separate `i64` exponent,
//! allowing values far outside the normal `f64` range.
//!
//! Clean-room implementation of the subset of operations used by `lll/l2.c`;
//! the upstream `dpe.h` is LGPL and is not used or referenced here.

use super::intbig::Ibz;
use rug::Assign;

/// Value represented is `m · 2^e`. After normalization either `m == 0` or
/// `0.5 <= |m| < 1`.
#[derive(Clone, Copy, Debug, Default)]
pub struct Dpe {
    m: f64,
    e: i64,
}

/// Split a finite `f64` into `(mantissa, exponent)` with `0.5 <= |mantissa| < 1`,
/// or `(0.0, 0)` for zero. Mirrors C `frexp`.
fn frexp(x: f64) -> (f64, i64) {
    if x == 0.0 {
        return (0.0, 0);
    }
    debug_assert!(x.is_finite());
    let bits = x.to_bits();
    let sign = bits & 0x8000_0000_0000_0000;
    let raw_exp = ((bits >> 52) & 0x7FF) as i64;
    let frac = bits & 0x000F_FFFF_FFFF_FFFF;
    if raw_exp == 0 {
        // Subnormal: scale up by 2^54 to normalize, then recurse.
        let (m, e) = frexp(x * (1u64 << 54) as f64);
        return (m, e - 54);
    }
    // Replace biased exponent with 1022 (so value lies in [0.5, 1)).
    let m = f64::from_bits(sign | (1022u64 << 52) | frac);
    (m, raw_exp - 1022)
}

/// `m · 2^e` as `f64` (may overflow to ±inf or underflow to 0).
fn ldexp(m: f64, e: i64) -> f64 {
    if m == 0.0 {
        return 0.0;
    }
    // Clamp to a range that f64 can represent via repeated scaling.
    let mut m = m;
    let mut e = e;
    while e > 1000 {
        m *= f64::from_bits((1023 + 1000) << 52);
        e -= 1000;
    }
    while e < -1000 {
        m *= f64::from_bits((1023 - 1000) << 52);
        e += 1000;
    }
    m * f64::from_bits(((1023 + e) as u64) << 52)
}

impl Dpe {
    pub const ZERO: Dpe = Dpe { m: 0.0, e: 0 };

    #[inline]
    fn normalize(&mut self) {
        let (m, de) = frexp(self.m);
        self.m = m;
        self.e = if m == 0.0 { 0 } else { self.e + de };
    }

    pub fn set_d(&mut self, d: f64) {
        self.m = d;
        self.e = 0;
        self.normalize();
    }

    pub fn set_ui(&mut self, u: u64) {
        self.set_d(u as f64);
    }

    /// Load from a big integer with `mpz_get_d_2exp` semantics: mantissa in
    /// `[0.5, 1)` truncated (not rounded) to 53 bits, plus the exact exponent.
    pub fn set_z(&mut self, z: &Ibz) {
        if z.is_zero() {
            *self = Dpe::ZERO;
            return;
        }
        let bits = z.significant_bits() as i64;
        // Extract exactly the top 53 bits of |z| (truncating any lower bits).
        let shift = (bits - 53).max(0);
        let top: Ibz = z.clone().abs() >> (shift as u32);
        let top_u = top.to_u64().expect("≤53 bits");
        // top_u has ≤53 bits → conversion to f64 is exact.
        let m = top_u as f64;
        let (mn, me) = frexp(m);
        self.m = if z.is_negative() { -mn } else { mn };
        self.e = shift + me;
    }

    /// Convert to a big integer (truncating any fractional part).
    pub fn get_z(&self, out: &mut Ibz) {
        if self.m == 0.0 {
            out.assign(0);
            return;
        }
        // Mantissa has 53 significant bits; scale to an integer then shift.
        let scaled = self.m * (1u64 << 53) as f64;
        let int_m = scaled.trunc() as i64;
        out.assign(int_m);
        let shift = self.e - 53;
        if shift >= 0 {
            *out <<= shift as u32;
        } else {
            // Truncating division by power of two.
            let neg = out.is_negative();
            if neg {
                let t = out.clone();
                out.assign(-t);
            }
            *out >>= (-shift) as u32;
            if neg {
                let t = out.clone();
                out.assign(-t);
            }
        }
    }

    pub fn set(&mut self, src: &Dpe) {
        *self = *src;
    }

    pub fn abs(&mut self, src: &Dpe) {
        self.m = src.m.abs();
        self.e = src.e;
    }

    pub fn cmp(&self, other: &Dpe) -> i32 {
        let a = self.m.signum();
        let b = other.m.signum();
        if a != b {
            return if a < b { -1 } else { 1 };
        }
        if self.m == 0.0 {
            return 0;
        }
        // Same sign, both nonzero, both normalized: compare exponents first.
        if self.e != other.e {
            let r = if self.e < other.e { -1 } else { 1 };
            return if a < 0.0 { -r } else { r };
        }
        if self.m < other.m {
            -1
        } else if self.m > other.m {
            1
        } else {
            0
        }
    }

    pub fn cmp_d(&self, d: f64) -> i32 {
        let mut tmp = Dpe::ZERO;
        tmp.set_d(d);
        self.cmp(&tmp)
    }

    pub fn mul(&mut self, a: &Dpe, b: &Dpe) {
        self.m = a.m * b.m;
        self.e = a.e + b.e;
        self.normalize();
    }

    pub fn sub(&mut self, a: &Dpe, b: &Dpe) {
        if b.m == 0.0 {
            *self = *a;
            return;
        }
        if a.m == 0.0 {
            self.m = -b.m;
            self.e = b.e;
            return;
        }
        // Align to the larger exponent and subtract.
        if a.e >= b.e {
            let d = (a.e - b.e).min(120);
            self.m = a.m - ldexp(b.m, -d);
            self.e = a.e;
        } else {
            let d = (b.e - a.e).min(120);
            self.m = ldexp(a.m, -d) - b.m;
            self.e = b.e;
        }
        self.normalize();
    }

    pub fn div(&mut self, a: &Dpe, b: &Dpe) {
        debug_assert!(b.m != 0.0);
        self.m = a.m / b.m;
        self.e = a.e - b.e;
        self.normalize();
    }

    /// Round to nearest integer (in-place semantics matching C usage:
    /// `dpe_round(x, x)`).
    pub fn round(&mut self, src: &Dpe) {
        if src.m == 0.0 {
            *self = Dpe::ZERO;
            return;
        }
        if src.e >= 53 {
            // Already an integer at this precision.
            *self = *src;
            return;
        }
        let v = ldexp(src.m, src.e);
        // The reference DPE rounds with C `round()` (half-away-from-zero),
        // which matches Rust's `f64::round()`.
        self.set_d(v.round());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-10 * b.abs().max(1.0)
    }

    fn val(d: &Dpe) -> f64 {
        ldexp(d.m, d.e)
    }

    #[test]
    fn frexp_ldexp_roundtrip() {
        for &x in &[0.0, 1.0, 2.5, -3.75, 1e-300, 1e300, f64::MIN_POSITIVE / 4.0] {
            let (m, e) = frexp(x);
            if x != 0.0 {
                assert!(m.abs() >= 0.5 && m.abs() < 1.0, "x={x} m={m}");
            }
            assert!(approx_eq(ldexp(m, e), x), "x={x}");
        }
    }

    #[test]
    fn arithmetic() {
        let mut a = Dpe::ZERO;
        a.set_d(3.0);
        let mut b = Dpe::ZERO;
        b.set_d(7.0);
        let mut r = Dpe::ZERO;
        r.mul(&a, &b);
        assert!(approx_eq(val(&r), 21.0));
        r.sub(&a, &b);
        assert!(approx_eq(val(&r), -4.0));
        r.div(&a, &b);
        assert!(approx_eq(val(&r), 3.0 / 7.0));
        assert!(a.cmp(&b) < 0);
        assert!(b.cmp(&a) > 0);
        assert!(a.cmp_d(3.0) == 0);
    }

    #[test]
    fn round_and_get_z() {
        let mut x = Dpe::ZERO;
        x.set_d(2.5);
        let mut r = Dpe::ZERO;
        r.round(&x);
        let mut z = Ibz::new();
        r.get_z(&mut z);
        assert_eq!(z, 3); // half-away-from-zero
        x.set_d(4.5);
        r.round(&x);
        r.get_z(&mut z);
        assert_eq!(z, 5);
        x.set_d(-2.5);
        r.round(&x);
        r.get_z(&mut z);
        assert_eq!(z, -3);
        x.set_d(-2.7);
        r.round(&x);
        r.get_z(&mut z);
        assert_eq!(z, -3);
    }

    #[test]
    fn set_z_get_z_roundtrip() {
        let big = Integer::from_str_radix("123456789012345678901234567890", 10).unwrap();
        let mut d = Dpe::ZERO;
        d.set_z(&big);
        let mut back = Ibz::new();
        d.get_z(&mut back);
        // Only ~53 bits of precision survive.
        let diff: Ibz = (big.clone() - &back).abs();
        let allow: Ibz = big.clone() >> 50;
        assert!(diff <= allow, "diff={diff}");
        // Small ints are exact.
        d.set_z(&Integer::from(-12345));
        d.get_z(&mut back);
        assert_eq!(back, -12345);
    }

    #[test]
    fn huge_exponents() {
        // Values way outside f64 range.
        let mut a = Dpe::ZERO;
        let mut b = Dpe::ZERO;
        a.set_z(&(Ibz::from(1) << 5000));
        b.set_z(&(Ibz::from(3) << 4999));
        let mut r = Dpe::ZERO;
        r.div(&a, &b);
        assert!(approx_eq(val(&r), 2.0 / 3.0));
        assert!(a.cmp(&b) < 0);
    }
}
