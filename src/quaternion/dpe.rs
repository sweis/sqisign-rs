// SPDX-License-Identifier: Apache-2.0
//! Double-plus-exponent: an `f64` mantissa with a separate `i64` exponent,
//! allowing values far outside the normal `f64` range.
//!
//! Clean-room implementation of the subset of operations used by `lll/l2.c`;
//! the upstream `dpe.h` is LGPL and is not used or referenced here.

use super::intbig::{ibz_get_d_2exp, ibz_set, ibz_set_from_mantissa_shift, Ibz};

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
        let (d, e) = ibz_get_d_2exp(z);
        self.m = d;
        self.e = e;
    }

    /// Convert to a big integer (truncating any fractional part).
    pub fn get_z(&self, out: &mut Ibz) {
        if self.m == 0.0 {
            ibz_set(out, 0);
            return;
        }
        // Mantissa has 53 significant bits; scale to an integer then shift.
        let scaled = self.m * (1u64 << 53) as f64;
        let int_m = scaled.trunc() as i64;
        ibz_set_from_mantissa_shift(out, int_m, self.e - 53);
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
        } else {
            i32::from(self.m > other.m)
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
    use super::super::intbig::{ibz_abs, ibz_div_2exp, ibz_from_i64, ibz_pow, ibz_set_from_str};
    use super::*;

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
        let mut r = Dpe::ZERO;
        let mut z = Ibz::default();
        for &(v, expect) in &[(2.5, 3), (4.5, 5), (-2.5, -3), (-2.7, -3)] {
            x.set_d(v);
            r.round(&x);
            r.get_z(&mut z);
            assert_eq!(z, ibz_from_i64(expect), "round({v})");
        }
    }

    #[test]
    fn set_z_get_z_roundtrip() {
        let mut big = Ibz::default();
        ibz_set_from_str(&mut big, "123456789012345678901234567890", 10);
        let mut d = Dpe::ZERO;
        d.set_z(&big);
        let mut back = Ibz::default();
        d.get_z(&mut back);
        let mut diff = Ibz::default();
        let mut delta = Ibz::default();
        super::super::intbig::ibz_sub(&mut delta, &big, &back);
        ibz_abs(&mut diff, &delta);
        let mut allow = Ibz::default();
        ibz_div_2exp(&mut allow, &big, 50);
        assert!(diff <= allow, "diff={diff}");
        d.set_z(&ibz_from_i64(-12345));
        d.get_z(&mut back);
        assert_eq!(back, ibz_from_i64(-12345));
    }

    #[test]
    fn huge_exponents() {
        let mut a = Dpe::ZERO;
        let mut b = Dpe::ZERO;
        let mut p = Ibz::default();
        ibz_pow(&mut p, &ibz_from_i64(2), 5000);
        a.set_z(&p);
        let mut p3 = Ibz::default();
        super::super::intbig::ibz_mul(&mut p3, &p, &ibz_from_i64(3));
        ibz_div_2exp(&mut p, &p3, 1);
        b.set_z(&p);
        let mut r = Dpe::ZERO;
        r.div(&a, &b);
        assert!(approx_eq(val(&r), 2.0 / 3.0));
        assert!(a.cmp(&b) < 0);
    }
}
