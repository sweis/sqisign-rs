// SPDX-License-Identifier: Apache-2.0
//! GF(p) arithmetic for the SQIsign lvl1 prime p = 5·2²⁴⁸ - 1.
//!
//! Ported from the auto-generated `fp_p5248_64.c` (modarith / Montgomery,
//! unsaturated radix-2⁵¹ representation, 5 limbs).
//!
//! All values are kept in Montgomery form internally; `encode`/`decode`
//! convert between canonical little-endian bytes and the internal form.

#![allow(clippy::needless_range_loop)]

use core::fmt;

// ---------------------------------------------------------------------------
// Parameters (lvl1, RADIX_64, ref impl)
// ---------------------------------------------------------------------------

/// Word type used by the C reference (`digit_t` for RADIX_64).
pub type Digit = u64;

/// Number of 51-bit limbs in an Fp element (`NWORDS_FIELD` for ref/64-bit).
pub const NWORDS_FIELD: usize = 5;
/// Number of 64-bit words to hold a scalar of order ~2²⁴⁸ (`NWORDS_ORDER`).
pub const NWORDS_ORDER: usize = 4;
/// `RADIX` from `tutil.h`.
pub const RADIX: u32 = 64;
/// `BITS` from `fp_constants.h`.
pub const BITS: u32 = 256;
/// `LOG2P` from `fp_constants.h`.
pub const LOG2P: u32 = 8;
/// `FP_ENCODED_BYTES` (= `Nbytes` in modarith output).
pub const FP_ENCODED_BYTES: usize = 32;

// modarith internal parameters
const LIMB_BITS: u32 = 51;
const NLIMBS: usize = 5;
const NBYTES: usize = 32;
const MASK: u64 = (1u64 << LIMB_BITS) - 1;
// p = 5·2²⁴⁸ - 1; in radix-2⁵¹ the top limb of p is 5·2⁴⁴ = 0x5000_0000_0000.
const P4: u64 = 0x5000_0000_0000;

/// An element of GF(p), in Montgomery form, unsaturated 5×51-bit limbs.
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
            write!(f, "{:02x}", b)?;
        }
        write!(f, ")")
    }
}

// ---------------------------------------------------------------------------
// Precomputed constants (Montgomery form)
// ---------------------------------------------------------------------------

/// Montgomery representation of 0.
pub const ZERO: Fp = Fp([0, 0, 0, 0, 0]);

/// Montgomery representation of 1.
pub const ONE: Fp = Fp([
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
]);

/// Montgomery representation of 2⁻¹.
const TWO_INV: Fp = Fp([
    0x0000_0000_0000_000c,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_4000_0000_0000,
]);

/// Montgomery representation of 3⁻¹.
const THREE_INV: Fp = Fp([
    0x0005_5555_5555_555d,
    0x0002_aaaa_aaaa_aaaa,
    0x0005_5555_5555_5555,
    0x0002_aaaa_aaaa_aaaa,
    0x0000_4555_5555_5555,
]);

/// Montgomery representation of 2²⁵⁶ (i.e. R²).
const R2: Fp = Fp([
    0x0001_9999_9999_9eb8,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0004_cccc_cccc_cccc,
    0x0000_1999_9999_9999,
]);

/// nres constant: R² mod p in radix form (for `nres`).
const NRES_C: [u64; NLIMBS] = [
    0x4_cccc_cccc_cf5c,
    0x1_9999_9999_9999,
    0x3_3333_3333_3333,
    0x6_6666_6666_6666,
    0x0_0ccc_cccc_cccc,
];

// ===========================================================================
// Internal modarith primitives (port of generated code in fp_p5248_64.c)
// ===========================================================================

/// Propagate carries through the unsaturated limbs.
/// Returns an all-ones mask if the result is negative (top bit of limb 4 set).
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
    // -((n[4] >> 1) >> 62)  →  all-ones if bit 63 set, else 0
    ((n[NLIMBS - 1] >> 1) >> 62).wrapping_neg()
}

/// Propagate carries and conditionally add p if negative.
/// Returns 1 if the value was negative (correction applied), else 0.
#[inline]
fn flatten(n: &mut [u64; NLIMBS]) -> i32 {
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(1u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(P4 & carry);
    let _ = prop(n);
    (carry & 1) as i32
}

/// Montgomery final subtract: n ← n - p, then add p back if it underflowed.
/// Returns 1 if n < p (i.e. was already canonical).
#[inline]
fn modfsb(n: &mut [u64; NLIMBS]) -> i32 {
    n[0] = n[0].wrapping_add(1);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_sub(P4);
    flatten(n)
}

/// n = a + b mod 2p.
#[inline]
fn modadd(a: &[u64; NLIMBS], b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = a[i].wrapping_add(b[i]);
    }
    // subtract 2p
    n[0] = n[0].wrapping_add(2);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_sub(0xa000_0000_0000);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(0xa000_0000_0000 & carry);
    let _ = prop(n);
}

/// n = a - b mod 2p.
#[inline]
fn modsub(a: &[u64; NLIMBS], b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = a[i].wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(0xa000_0000_0000 & carry);
    let _ = prop(n);
}

/// n = -b mod 2p.
#[inline]
fn modneg(b: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    for i in 0..NLIMBS {
        n[i] = 0u64.wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2u64 & carry);
    n[NLIMBS - 1] = n[NLIMBS - 1].wrapping_add(0xa000_0000_0000 & carry);
    let _ = prop(n);
}

/// c = a * b mod 2p (Montgomery).
#[inline]
fn modmul(a: &[u64; NLIMBS], b: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p4 = P4 as u128;
    let mut t: u128 = 0;

    t += (a[0] as u128) * (b[0] as u128);
    let v0 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[1] as u128);
    t += (a[1] as u128) * (b[0] as u128);
    let v1 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[2] as u128);
    t += (a[1] as u128) * (b[1] as u128);
    t += (a[2] as u128) * (b[0] as u128);
    let v2 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[3] as u128);
    t += (a[1] as u128) * (b[2] as u128);
    t += (a[2] as u128) * (b[1] as u128);
    t += (a[3] as u128) * (b[0] as u128);
    let v3 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[4] as u128);
    t += (a[1] as u128) * (b[3] as u128);
    t += (a[2] as u128) * (b[2] as u128);
    t += (a[3] as u128) * (b[1] as u128);
    t += (a[4] as u128) * (b[0] as u128);
    t += (v0 as u128) * p4;
    let v4 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[1] as u128) * (b[4] as u128);
    t += (a[2] as u128) * (b[3] as u128);
    t += (a[3] as u128) * (b[2] as u128);
    t += (a[4] as u128) * (b[1] as u128);
    t += (v1 as u128) * p4;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[2] as u128) * (b[4] as u128);
    t += (a[3] as u128) * (b[3] as u128);
    t += (a[4] as u128) * (b[2] as u128);
    t += (v2 as u128) * p4;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[3] as u128) * (b[4] as u128);
    t += (a[4] as u128) * (b[3] as u128);
    t += (v3 as u128) * p4;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[4] as u128) * (b[4] as u128);
    t += (v4 as u128) * p4;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[4] = t as u64;
}

/// c = a² mod 2p (Montgomery).
#[inline]
fn modsqr(a: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p4 = P4 as u128;
    let mut t: u128;
    let mut tot: u128;

    tot = (a[0] as u128) * (a[0] as u128);
    t = tot;
    let v0 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[1] as u128);
    tot *= 2;
    t += tot;
    let v1 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[2] as u128);
    tot *= 2;
    tot += (a[1] as u128) * (a[1] as u128);
    t += tot;
    let v2 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[3] as u128);
    tot += (a[1] as u128) * (a[2] as u128);
    tot *= 2;
    t += tot;
    let v3 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[4] as u128);
    tot += (a[1] as u128) * (a[3] as u128);
    tot *= 2;
    tot += (a[2] as u128) * (a[2] as u128);
    t += tot;
    t += (v0 as u128) * p4;
    let v4 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[1] as u128) * (a[4] as u128);
    tot += (a[2] as u128) * (a[3] as u128);
    tot *= 2;
    t += tot;
    t += (v1 as u128) * p4;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[2] as u128) * (a[4] as u128);
    tot *= 2;
    tot += (a[3] as u128) * (a[3] as u128);
    t += tot;
    t += (v2 as u128) * p4;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[3] as u128) * (a[4] as u128);
    tot *= 2;
    t += tot;
    t += (v3 as u128) * p4;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[4] as u128) * (a[4] as u128);
    t += tot;
    t += (v4 as u128) * p4;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[4] = t as u64;
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

/// `c = a * b mod 2p`, with aliasing allowed (matches C semantics).
#[inline]
fn modmul_ip(a: [u64; NLIMBS], b: [u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    modmul(&a, &b, c);
}

/// Progenitor: z = w^((p-3)/4). Addition chain from modarith.
fn modpro(w: &[u64; NLIMBS], z: &mut [u64; NLIMBS]) {
    let mut x = [0u64; NLIMBS];
    let mut t0 = [0u64; NLIMBS];
    let mut t1 = [0u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4 = [0u64; NLIMBS];
    modcpy(w, &mut x);
    modsqr(&x, z);
    modmul_ip(x, *z, &mut t0);
    modsqr(&t0, z);
    modmul_ip(x, *z, z);
    modsqr(z, &mut t1);
    modsqr(&t1, &mut t3);
    modsqr(&t3, &mut t2);
    modcpy(&t2, &mut t4);
    modnsqr(&mut t4, 3);
    modmul_ip(t2, t4, &mut t2);
    modcpy(&t2, &mut t4);
    modnsqr(&mut t4, 6);
    modmul_ip(t2, t4, &mut t2);
    modcpy(&t2, &mut t4);
    modnsqr(&mut t4, 2);
    modmul_ip(t3, t4, &mut t3);
    modnsqr(&mut t3, 13);
    modmul_ip(t2, t3, &mut t2);
    modcpy(&t2, &mut t3);
    modnsqr(&mut t3, 27);
    modmul_ip(t2, t3, &mut t2);
    modmul_ip(*z, t2, z);
    modcpy(z, &mut t2);
    modnsqr(&mut t2, 4);
    modmul_ip(t1, t2, &mut t1);
    modmul_ip(t0, t1, &mut t0);
    modmul_ip(t1, t0, &mut t1);
    modmul_ip(t0, t1, &mut t0);
    modmul_ip(t1, t0, &mut t2);
    modmul_ip(t0, t2, &mut t0);
    modmul_ip(t1, t0, &mut t1);
    modnsqr(&mut t1, 63);
    modmul_ip(t0, t1, &mut t1);
    modnsqr(&mut t1, 64);
    modmul_ip(t0, t1, &mut t0);
    modnsqr(&mut t0, 57);
    modmul_ip(*z, t0, z);
}

/// z = x⁻¹ (using optional precomputed progenitor h).
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

/// Convert m to Montgomery form: n = m·R mod p.
#[inline]
fn nres(m: &[u64; NLIMBS], n: &mut [u64; NLIMBS]) {
    modmul(m, &NRES_C, n);
}

/// Convert n out of Montgomery form: m = n·R⁻¹ mod p, fully reduced.
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
    *a = [1, 0, 0, 0, 0];
    let t = *a;
    nres(&t, a);
}

#[inline]
fn modint(x: i32, a: &mut [u64; NLIMBS]) {
    *a = [x as u64, 0, 0, 0, 0];
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

/// Conditional move: if b==1, f ← g.
#[inline(never)]
fn modcmv(b: i32, g: &[u64; NLIMBS], f: &mut [u64; NLIMBS]) {
    let r: u64 = 0x3cc3_c33c_5aa5_a55a;
    let c0 = (1u64.wrapping_sub(b as u64)).wrapping_add(r);
    let c1 = (b as u64).wrapping_add(r);
    for i in 0..NLIMBS {
        let s = g[i];
        let t = f[i];
        let v = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s));
        f[i] = v.wrapping_sub(r.wrapping_mul(t.wrapping_add(s)));
    }
}

/// Conditional swap: if b==1, swap f and g.
#[inline(never)]
fn modcsw(b: i32, g: &mut [u64; NLIMBS], f: &mut [u64; NLIMBS]) {
    let r: u64 = 0x3cc3_c33c_5aa5_a55a;
    let c0 = (1u64.wrapping_sub(b as u64)).wrapping_add(r);
    let c1 = (b as u64).wrapping_add(r);
    for i in 0..NLIMBS {
        let s = g[i];
        let t = f[i];
        let w = r.wrapping_mul(t.wrapping_add(s));
        f[i] = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s)).wrapping_sub(w);
        g[i] = c0.wrapping_mul(s).wrapping_add(c1.wrapping_mul(t)).wrapping_sub(w);
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

#[allow(dead_code)]
fn mod2r(r: u32, a: &mut [u64; NLIMBS]) {
    let n = (r / LIMB_BITS) as usize;
    let m = r % LIMB_BITS;
    modzer(a);
    if r >= 32 * 8 {
        return;
    }
    a[n] = 1u64 << m;
    let t = *a;
    nres(&t, a);
}

#[allow(dead_code)]
fn modsign(a: &[u64; NLIMBS]) -> i32 {
    let mut c = [0u64; NLIMBS];
    redc(a, &mut c);
    (c[0] % 2) as i32
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

/// Returns 0xFFFFFFFF if a == b, 0 otherwise (constant time).
#[inline]
pub fn fp_is_equal(a: &Fp, b: &Fp) -> u32 {
    (modcmp(&a.0, &b.0) as u32).wrapping_neg()
}

/// Returns 0xFFFFFFFF if a == 0, 0 otherwise (constant time).
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

#[inline]
pub fn fp_inv(x: &mut Fp) {
    let a = x.0;
    modinv(&a, None, &mut x.0);
}

/// Returns 0xFFFFFFFF if a is a square in GF(p), 0 otherwise.
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

/// out = a^((p-3)/4).
#[inline]
pub fn fp_exp3div4(out: &mut Fp, a: &Fp) {
    modpro(&a.0, &mut out.0);
}

#[inline]
pub fn fp_div3(out: &mut Fp, a: &Fp) {
    modmul(&THREE_INV.0, &a.0, &mut out.0);
}

/// `ctl` MUST be 0 or 0xFFFFFFFF. d = ctl ? a1 : a0.
#[inline]
pub fn fp_select(d: &mut Fp, a0: &Fp, a1: &Fp, ctl: u32) {
    // Sign-extend the 32-bit mask to a 64-bit word, matching `(int32_t)ctl` in C.
    let cw = (ctl as i32) as i64 as u64;
    for i in 0..NWORDS_FIELD {
        d.0[i] = a0.0[i] ^ (cw & (a0.0[i] ^ a1.0[i]));
    }
}

/// Encode a (in Montgomery form) into 32 little-endian canonical bytes.
pub fn fp_encode(dst: &mut [u8], a: &Fp) {
    debug_assert!(dst.len() >= FP_ENCODED_BYTES);
    let mut c = [0u64; NLIMBS];
    redc(&a.0, &mut c);
    for i in 0..NBYTES {
        dst[i] = (c[0] & 0xff) as u8;
        let _ = modshr(8, &mut c);
    }
}

/// Decode 32 little-endian bytes into Montgomery form.
/// Returns 0xFFFFFFFF if the input was a canonical encoding (< p), else 0
/// and the output is forced to zero.
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

#[inline]
fn add_carry(cc: u8, a: u64, b: u64, d: &mut u64) -> u8 {
    let t = (a as u128) + (b as u128) + (cc as u128);
    *d = t as u64;
    (t >> 64) as u8
}

/// Partial reduce a 4×64-bit packed value modulo p using 5·2²⁴⁸ ≡ 1.
fn partial_reduce(out: &mut [u64; 4], src: &[u64; 4]) {
    let h = src[3] >> 56;
    let l = src[3] & 0x00FF_FFFF_FFFF_FFFF;
    // 5·2²⁴⁸ ≡ 1 ⇒ add ⌊h/5⌋ + (h mod 5)·2²⁴⁸ to the low part.
    let quo = (h * 0xCD) >> 10;
    let rem = h - 5 * quo;
    let mut o0 = 0u64;
    let mut o1 = 0u64;
    let mut o2 = 0u64;
    let mut o3 = 0u64;
    let mut cc = add_carry(0, src[0], quo, &mut o0);
    cc = add_carry(cc, src[1], 0, &mut o1);
    cc = add_carry(cc, src[2], 0, &mut o2);
    let _ = add_carry(cc, l, rem << 56, &mut o3);
    out[0] = o0;
    out[1] = o1;
    out[2] = o2;
    out[3] = o3;
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
        for b in tmp[rem..].iter_mut() {
            *b = 0;
        }
        fp_decode(d, &tmp);
        len = k;
    }

    while len > 0 {
        let prev = *d;
        fp_mul(d, &prev, &R2);
        len -= NBYTES;
        let mut t = [0u64; 4];
        for j in 0..4 {
            t[j] = u64::from_le_bytes(src[len + 8 * j..len + 8 * j + 8].try_into().unwrap());
        }
        let tin = t;
        partial_reduce(&mut t, &tin);
        for j in 0..4 {
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

    /// Deterministic xorshift64* RNG for tests (independent of randombytes).
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

    fn fp_random(prng: &mut Prng) -> Fp {
        let mut buf = [0u8; FP_ENCODED_BYTES];
        prng.fill(&mut buf);
        let mut a = Fp::default();
        fp_decode_reduce(&mut a, &buf);
        a
    }

    const ITERS: usize = 500;

    #[test]
    fn one_constant_matches_set_one() {
        let mut a = Fp::default();
        fp_set_one(&mut a);
        // The internal limbs should equal the precomputed ONE constant.
        // Both are in Montgomery form; modone computes nres([1,0,0,0,0]).
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
            let val = (prng.next() as u32) & 0x7FFF_FFFF; // keep positive as i32
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
        // All 0xFF bytes is far above p.
        let buf = [0xFFu8; FP_ENCODED_BYTES];
        let mut d = Fp::default();
        let ok = fp_decode(&mut d, &buf);
        assert_eq!(ok, 0);
        assert_ne!(fp_is_zero(&d), 0);
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
        let mut enc = [0u8; FP_ENCODED_BYTES];

        let mut a = Fp::default();
        fp_set_one(&mut a);
        fp_encode(&mut enc, &a);
        assert_hex(&enc, "0100000000000000000000000000000000000000000000000000000000000000");

        let mut b = Fp::default();
        let mut c = Fp::default();
        fp_set_small(&mut a, 5);
        fp_set_small(&mut b, 7);
        fp_mul(&mut c, &a, &b);
        fp_encode(&mut enc, &c);
        assert_hex(&enc, "2300000000000000000000000000000000000000000000000000000000000000");

        fp_inv(&mut c);
        fp_encode(&mut enc, &c);
        assert_hex(&enc, "9224499224499224499224499224499224499224499224499224499224499201");

        fp_set_small(&mut a, 12345);
        fp_sqr(&mut b, &a);
        fp_sqrt(&mut b);
        fp_encode(&mut enc, &b);
        assert_hex(&enc, "c6cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffff04");

        let mut src = [0u8; 64];
        for i in 0..64 {
            src[i] = (i as u8).wrapping_mul(7).wrapping_add(1);
        }
        fp_decode_reduce(&mut a, &src);
        fp_encode(&mut enc, &a);
        assert_hex(&enc, "a43cd712e8565f34a3ab4d895ecdd577b388f7ffa1ddb2212acc07dd4b54f602");

        fp_decode_reduce(&mut a, &src[..47]);
        fp_encode(&mut enc, &a);
        assert_hex(&enc, "2c35d712e8565f34a3ab4d895ecdd57771787f868d949ba2a9b0b7bec5ccd303");
    }

    #[test]
    fn precomp_fp_constants_match() {
        // The precomp module's internal-form constants must equal what gf computes.
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
