// SPDX-License-Identifier: Apache-2.0
//! Big-integer wrapper, ported from `intbig.c`.
//!
//! All `ibz_*` functions mirror the C signatures (out-params first) so the rest
//! of the quaternion module can be ported mechanically. The actual arithmetic
//! is delegated to a backend module: `malachite` by default (pure Rust), or
//! `rug`/GMP under `feature = "gmp"`. Only this file and the two backend files
//! reference the underlying bigint crate; everything else uses `Ibz` opaquely.

#[cfg(feature = "gmp")]
#[path = "intbig_backend_rug.rs"]
mod backend;
#[cfg(all(not(feature = "gmp"), feature = "cryptobigint"))]
#[path = "intbig_backend_cryptobigint.rs"]
mod backend;
#[cfg(all(not(feature = "gmp"), not(feature = "cryptobigint")))]
#[path = "intbig_backend_malachite.rs"]
mod backend;

pub use backend::*;

use crate::common::ctrdrbg::randombytes;
use crate::mp::Digit;
use std::sync::OnceLock;

// ---------------------------------------------------------------------------
// Constants

fn const_int(v: i32) -> &'static Ibz {
    static CELLS: [OnceLock<Ibz>; 4] = [
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
    ];
    CELLS[v as usize].get_or_init(|| ibz_from_i64(i64::from(v)))
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

#[inline]
pub fn ibz_init(x: &mut Ibz) {
    *x = ibz_new();
}
#[inline]
pub fn ibz_finalize(_x: &mut Ibz) {}

// ---------------------------------------------------------------------------
// Random sampling (deterministic via the KAT DRBG). Backend-agnostic; all
// big-integer access goes through `ibz_*` so the byte-consumption pattern is
// identical across backends.

/// Uniform sample in `[a, b]` by rejection. Matches the C byte-sampling loop.
pub fn ibz_rand_interval(rand: &mut Ibz, a: &Ibz, b: &Ibz) -> i32 {
    let mut bmina = ibz_new();
    ibz_sub(&mut bmina, b, a);
    if ibz_is_zero(&bmina) != 0 {
        ibz_copy(rand, a);
        return 1;
    }
    debug_assert!(ibz_is_positive(&bmina));

    let len_bits = ibz_significant_bits(&bmina) as usize;
    let len_bytes = (len_bits + 7) / 8;
    const LIMB: usize = core::mem::size_of::<Digit>();
    const LIMB_BITS: usize = LIMB * 8;
    let len_limbs = (len_bytes + LIMB - 1) / LIMB;

    let shift = (LIMB_BITS - len_bits % LIMB_BITS) % LIMB_BITS;
    let mask: Digit = (!0u64) >> shift;

    let mut buf = vec![0u8; len_limbs * LIMB];
    let mut tmp = ibz_new();
    let mut limbs = vec![0u64; len_limbs];
    loop {
        // C draws exactly len_bytes; the upper bytes of the top limb are
        // discarded by the mask, so their initial value is irrelevant.
        randombytes(&mut buf[..len_bytes]);
        for (i, chunk) in buf.chunks(LIMB).enumerate() {
            let mut arr = [0u8; LIMB];
            arr[..chunk.len()].copy_from_slice(chunk);
            limbs[i] = u64::from_le_bytes(arr);
        }
        limbs[len_limbs - 1] &= mask;
        ibz_copy_digits(&mut tmp, &limbs);
        if ibz_cmp(&tmp, &bmina) <= 0 {
            break;
        }
    }
    ibz_add(rand, &tmp, a);
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
    let two_m = ibz_from_i64(2 * i64::from(m));
    let ret = ibz_rand_interval(rand, ibz_const_zero(), &two_m);
    let r = rand.clone();
    ibz_sub(rand, &r, &ibz_from_i64(i64::from(m)));
    ret
}

/// Uniform sample in `[-2^m, 2^m]`.
pub fn ibz_rand_interval_bits(rand: &mut Ibz, m: u32) -> i32 {
    let mut hi = ibz_new();
    ibz_pow(&mut hi, ibz_const_two(), m);
    let mut lo = ibz_new();
    ibz_neg(&mut lo, &hi);
    let ret = ibz_rand_interval(rand, &lo, &hi);
    if ret != 1 {
        return ret;
    }
    // Match the C code's trailing `mpz_sub_ui(*rand, *rand, m)` exactly.
    let r = rand.clone();
    ibz_sub(rand, &r, &ibz_from_i64(i64::from(m)));
    ret
}

// ---------------------------------------------------------------------------
// Tonelli–Shanks modular square root. Backend-agnostic.

pub fn ibz_sqrt_mod_p(sqrt: &mut Ibz, a: &Ibz, p: &Ibz) -> i32 {
    debug_assert!(ibz_probab_prime(p, 30) != 0);

    let mut amod = ibz_new();
    ibz_mod(&mut amod, a, p);
    if ibz_is_zero(&amod) != 0 {
        ibz_set(sqrt, 0);
        // C falls through; legendre(0,p)=0 so it returns 0.
    }
    if ibz_legendre(&amod, p) != 1 {
        return 0;
    }

    let mut pm1 = ibz_new();
    ibz_sub(&mut pm1, p, ibz_const_one());
    let mut tmp = ibz_new();

    let p_mod_8 = ibz_mod_ui(p, 8);
    if p_mod_8 % 4 == 3 {
        // sqrt = a^{(p+1)/4}
        let mut e = ibz_new();
        ibz_add(&mut e, p, ibz_const_one());
        let ee = e.clone();
        ibz_div_2exp(&mut e, &ee, 2);
        ibz_pow_mod(sqrt, &amod, &e, p);
    } else if p_mod_8 == 5 {
        // Atkin's trick.
        let mut e = ibz_new();
        ibz_div_2exp(&mut e, &pm1, 2);
        let mut t = ibz_new();
        ibz_pow_mod(&mut t, &amod, &e, p);
        if ibz_is_one(&t) != 0 {
            ibz_add(&mut e, p, ibz_const_three());
            let ee = e.clone();
            ibz_div_2exp(&mut e, &ee, 3);
            ibz_pow_mod(sqrt, &amod, &e, p);
        } else {
            let mut e2 = ibz_new();
            ibz_sub(&mut e2, p, &ibz_from_i64(5));
            let ee = e2.clone();
            ibz_div_2exp(&mut e2, &ee, 3);
            let mut a4 = ibz_new();
            ibz_mul(&mut a4, &amod, &ibz_from_i64(4));
            let mut r = ibz_new();
            ibz_pow_mod(&mut r, &a4, &e2, p);
            let mut a2 = ibz_new();
            ibz_mul(&mut a2, &amod, ibz_const_two());
            ibz_mul(&mut tmp, &r, &a2);
            ibz_mod(sqrt, &tmp, p);
        }
    } else {
        // p ≡ 1 (mod 8): full Tonelli–Shanks.
        let e = ibz_two_adic(&pm1) as u32;
        let mut q = ibz_new();
        ibz_div_2exp(&mut q, &pm1, e);

        // Find a non-residue.
        let mut qnr = ibz_from_i64(2);
        while ibz_legendre(&qnr, p) != -1 {
            let q1 = qnr.clone();
            ibz_add(&mut qnr, &q1, ibz_const_one());
        }
        let mut z = ibz_new();
        ibz_pow_mod(&mut z, &qnr, &q, p);

        let mut y = ibz_new();
        ibz_pow_mod(&mut y, &amod, &q, p);
        let mut half_q1 = ibz_new();
        ibz_add(&mut half_q1, &q, ibz_const_one());
        let h = half_q1.clone();
        ibz_div_2exp(&mut half_q1, &h, 1);
        let mut x = ibz_new();
        ibz_pow_mod(&mut x, &amod, &half_q1, p);

        let mut exp = ibz_new();
        ibz_pow(&mut exp, ibz_const_two(), e - 2);
        for _ in 0..e {
            let mut b = ibz_new();
            ibz_pow_mod(&mut b, &y, &exp, p);
            if ibz_cmp(&b, &pm1) == 0 {
                ibz_mul(&mut tmp, &x, &z);
                ibz_mod(&mut x, &tmp, p);
                ibz_mul(&mut tmp, &y, &z);
                let t = tmp.clone();
                ibz_mul(&mut tmp, &t, &z);
                ibz_mod(&mut y, &tmp, p);
            }
            // z = z² mod p
            let zz = z.clone();
            ibz_mul(&mut tmp, &zz, &zz);
            ibz_mod(&mut z, &tmp, p);
            let ee = exp.clone();
            ibz_div_2exp(&mut exp, &ee, 1);
        }
        ibz_copy(sqrt, &x);
    }
    1
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    fn z(v: i64) -> Ibz {
        ibz_from_i64(v)
    }

    #[test]
    fn arithmetic_basics() {
        let mut r = ibz_new();
        ibz_add(&mut r, &z(7), &z(-3));
        assert_eq!(ibz_cmp(&r, &z(4)), 0);
        ibz_sub(&mut r, &z(7), &z(-3));
        assert_eq!(ibz_cmp(&r, &z(10)), 0);
        ibz_mul(&mut r, &z(7), &z(-3));
        assert_eq!(ibz_cmp(&r, &z(-21)), 0);
        ibz_neg(&mut r, &z(-5));
        assert_eq!(ibz_cmp(&r, &z(5)), 0);
        ibz_abs(&mut r, &z(-9));
        assert_eq!(ibz_cmp(&r, &z(9)), 0);
    }

    #[test]
    fn div_tdiv_semantics() {
        let (mut q, mut r) = (ibz_new(), ibz_new());
        ibz_div(&mut q, &mut r, &z(-7), &z(3));
        assert_eq!((ibz_get(&q), ibz_get(&r)), (-2, -1));
        ibz_div(&mut q, &mut r, &z(7), &z(-3));
        assert_eq!((ibz_get(&q), ibz_get(&r)), (-2, 1));
    }

    #[test]
    fn div_floor_semantics() {
        let (mut q, mut r) = (ibz_new(), ibz_new());
        ibz_div_floor(&mut q, &mut r, &z(-7), &z(3));
        assert_eq!((ibz_get(&q), ibz_get(&r)), (-3, 2));
    }

    #[test]
    fn div_2exp_trunc() {
        let mut q = ibz_new();
        ibz_div_2exp(&mut q, &z(-7), 1);
        assert_eq!(ibz_get(&q), -3); // tdiv: -7/2 → -3
        ibz_div_2exp(&mut q, &z(7), 1);
        assert_eq!(ibz_get(&q), 3);
    }

    #[test]
    fn mod_nonneg() {
        let mut r = ibz_new();
        ibz_mod(&mut r, &z(-7), &z(3));
        assert_eq!(ibz_get(&r), 2);
        ibz_mod(&mut r, &z(-7), &z(-3));
        assert_eq!(ibz_get(&r), 2);
        ibz_mod(&mut r, &z(7), &z(3));
        assert_eq!(ibz_get(&r), 1);
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
        assert_eq!(ibz_two_adic(&z(12)), 2);
        assert_eq!(ibz_two_adic(&z(1)), 0);
        assert_eq!(ibz_bitsize(&z(255)), 8);
        assert_eq!(ibz_bitsize(&z(256)), 9);
    }

    #[test]
    fn ibz_get_low_bits() {
        let mut big = ibz_new();
        ibz_set_from_str(&mut big, "1234567890ABCDEF1234567890ABCDEF", 16);
        let g = ibz_get(&big);
        assert_eq!(g & 1, ibz_get_bit(&big, 0) as i32);
        assert_eq!((g & 3) as u64, ibz_mod_ui(&big, 4));
        assert_eq!(ibz_get(&z(42)), 42);
        assert_eq!(ibz_get(&z(-42)), -42);
        // Edge cases for mpz_get_si emulation: sign+magnitude truncation.
        // |2^40| → low 31 bits are 0.
        assert_eq!(ibz_get(&ibz_from_i64(1i64 << 40)), 0);
        assert_eq!(ibz_get(&ibz_from_i64(-(1i64 << 40))), i32::MIN);
        // 2^31 + 5 → low 31 bits = 5.
        assert_eq!(ibz_get(&ibz_from_i64((1i64 << 31) + 5)), 5);
    }

    #[test]
    fn digits_roundtrip() {
        let limbs: [u64; 3] = [0xDEAD_BEEF, 0xCAFE_BABE, 0x1];
        let mut x = ibz_new();
        ibz_copy_digits(&mut x, &limbs);
        let mut out = [0u64; 4];
        ibz_to_digits(&mut out, &x);
        assert_eq!(&out[..3], &limbs);
        assert_eq!(out[3], 0);
        let zero = ibz_new();
        let mut out2 = [99u64; 2];
        ibz_to_digits(&mut out2, &zero);
        assert_eq!(out2, [0, 0]);
    }

    #[test]
    fn gcd_invmod() {
        let mut g = ibz_new();
        ibz_gcd(&mut g, &z(462), &z(1071));
        assert_eq!(ibz_get(&g), 21);
        let mut inv = ibz_new();
        assert_eq!(ibz_invmod(&mut inv, &z(3), &z(7)), 1);
        assert_eq!(ibz_get(&inv), 5);
        assert_eq!(ibz_invmod(&mut inv, &z(4), &z(8)), 0);
    }

    #[test]
    fn sqrt_perfect() {
        let mut s = ibz_new();
        assert_eq!(ibz_sqrt(&mut s, &z(144)), 1);
        assert_eq!(ibz_get(&s), 12);
        assert_eq!(ibz_sqrt(&mut s, &z(145)), 0);
        ibz_sqrt_floor(&mut s, &z(145));
        assert_eq!(ibz_get(&s), 12);
    }

    #[test]
    fn sqrt_mod_p_all_branches() {
        let mut s = ibz_new();
        for &(a, p) in &[(2i64, 7i64), (3, 13), (10, 13), (2, 17), (8, 41)] {
            assert_eq!(ibz_sqrt_mod_p(&mut s, &z(a), &z(p)), 1, "a={a} p={p}");
            let mut s2 = ibz_new();
            ibz_mul(&mut s2, &s, &s);
            assert_eq!(ibz_mod_ui(&s2, p as u64), a as u64, "a={a} p={p}");
        }
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
        let mut r = ibz_new();
        ibz_pow(&mut r, &z(3), 5);
        assert_eq!(ibz_get(&r), 243);
        ibz_pow_mod(&mut r, &z(3), &z(5), &z(7));
        assert_eq!(ibz_get(&r), 5);
    }

    #[test]
    fn str_roundtrip() {
        let mut x = ibz_new();
        assert_eq!(ibz_set_from_str(&mut x, "-12345678901234567890", 10), 1);
        assert_eq!(ibz_convert_to_str(&x, 10).unwrap(), "-12345678901234567890");
    }

    #[test]
    fn d_2exp_semantics() {
        // Golden values: matches mpz_get_d_2exp truncation.
        let (d, e) = ibz_get_d_2exp(&z(1));
        assert_eq!((d, e), (0.5, 1));
        let (d, e) = ibz_get_d_2exp(&z(-3));
        assert_eq!((d, e), (-0.75, 2));
        // 2^60 - 1: top 53 bits truncated, not rounded.
        let (d, e) = ibz_get_d_2exp(&z((1i64 << 60) - 1));
        assert_eq!(e, 60);
        // mantissa is (2^53 - 1)/2^53 (truncated), not 1.0 (rounded).
        assert!((d - (1.0 - 1.0 / (1u64 << 53) as f64)).abs() < 1e-18);
    }

    #[test]
    fn primality_bpsw() {
        // Exercise the BPSW path with composites that pass single-base MR.
        // 2047 = 23·89 is a strong pseudoprime to base 2.
        assert_eq!(ibz_probab_prime(&z(2047), 30), 0);
        // Carmichael number 561 = 3·11·17.
        assert_eq!(ibz_probab_prime(&z(561), 30), 0);
        // Large prime.
        let mut p = ibz_new();
        ibz_set_from_str(
            &mut p,
            "170141183460469231731687303715884105727", // 2^127 - 1
            10,
        );
        assert!(ibz_probab_prime(&p, 30) > 0);
        // Large composite (product of two ~64-bit primes).
        let mut c = ibz_new();
        ibz_set_from_str(&mut c, "340282366920938463463374607431768211457", 10);
        // = (2^64+13)(2^64-?) ... actually just check it's not the Mersenne.
        ibz_mul(
            &mut c,
            &ibz_from_i64(1_000_000_007),
            &ibz_from_i64(998_244_353),
        );
        assert_eq!(ibz_probab_prime(&c, 30), 0);
    }

    // ---- Backend-convergence regressions (must agree with GMP semantics). ----

    #[test]
    fn probab_prime_tests_abs() {
        // GMP tests |n|; previous malachite/cryptobigint backends returned 0.
        assert!(ibz_probab_prime(&ibz_from_i64(-7), 20) > 0);
        assert!(ibz_probab_prime(&ibz_from_i64(-9), 20) == 0);
        assert!(ibz_probab_prime(&ibz_from_i64(-1), 20) == 0);
    }

    #[test]
    fn get_bit_twos_complement_on_negatives() {
        // -3 = ...11101₂. GMP/malachite use infinite two's-complement.
        let m3 = ibz_from_i64(-3);
        assert!(ibz_get_bit(&m3, 0));
        assert!(!ibz_get_bit(&m3, 1));
        assert!(ibz_get_bit(&m3, 2));
        assert!(ibz_get_bit(&m3, 5));
        assert!(ibz_get_bit(&m3, 63));
    }

    #[test]
    fn xgcd_matches_gmp_on_2g_edge() {
        // GMP convention: when |b| = 2g, u = sgn(a).
        let mut g = ibz_new();
        let mut u = ibz_new();
        let mut v = ibz_new();
        for &(a, b, eg, eu, ev) in &[
            (-3i64, 2i64, 1i64, -1i64, -1i64),
            (-6, 4, 2, -1, -1),
            (3, 2, 1, 1, -1),
            (6, -4, 2, 1, 1),
        ] {
            ibz_xgcd(&mut g, &mut u, &mut v, &ibz_from_i64(a), &ibz_from_i64(b));
            assert_eq!(ibz_cmp(&g, &ibz_from_i64(eg)), 0, "g({a},{b})");
            assert_eq!(ibz_cmp(&u, &ibz_from_i64(eu)), 0, "u({a},{b})");
            assert_eq!(ibz_cmp(&v, &ibz_from_i64(ev)), 0, "v({a},{b})");
        }
    }

    #[test]
    fn set_from_str_rejects_empty() {
        let mut i = ibz_new();
        assert_eq!(ibz_set_from_str(&mut i, "", 10), 0);
        assert_eq!(ibz_set_from_str(&mut i, "-", 10), 0);
    }

    #[test]
    fn secure_clear_yields_zero() {
        let mut x = ibz_new();
        ibz_set_from_str(&mut x, "123456789012345678901234567890", 10);
        ibz_secure_clear(&mut x);
        assert_eq!(ibz_cmp(&x, ibz_const_zero()), 0);
    }
}
