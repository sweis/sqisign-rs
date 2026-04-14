// SPDX-License-Identifier: Apache-2.0
//! GF(p) backend for the lvl1 prime p = 5·2²⁴⁸ − 1.
//! Ported from `fp_p5248_64.c` (modarith / Montgomery, unsaturated radix-2⁵¹).

#![allow(clippy::unreadable_literal)]

use super::fp::Fp;

pub const NWORDS_FIELD: usize = 5;
pub const NWORDS_ORDER: usize = 4;
pub const BITS: u32 = 256;
pub const LOG2P: u32 = 8;
pub const FP_ENCODED_BYTES: usize = 32;

pub const NLIMBS: usize = 5;
pub const LIMB_BITS: u32 = 51;
pub const NBYTES: usize = 32;
pub const MASK: u64 = (1u64 << LIMB_BITS) - 1;

/// Top limb of p in radix-2⁵¹ form.
pub const P_TOP: u64 = 0x5000_0000_0000;
/// Top limb of 2p.
pub const TWO_P_TOP: u64 = 0xa000_0000_0000;

/// Number of saturated 64-bit words in `partial_reduce` input/output.
pub const PR_WORDS: usize = 4;

pub const ZERO: Fp = Fp([0; NLIMBS]);
pub const ONE: Fp = Fp([
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
]);
pub const MINUS_ONE: Fp = Fp([
    0x7ffffffffffe6,
    0x7ffffffffffff,
    0x7ffffffffffff,
    0x7ffffffffffff,
    0x1fffffffffff,
]);
pub const TWO_INV: Fp = Fp([
    0x0000_0000_0000_000c,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_4000_0000_0000,
]);
pub const THREE_INV: Fp = Fp([
    0x0005_5555_5555_555d,
    0x0002_aaaa_aaaa_aaaa,
    0x0005_5555_5555_5555,
    0x0002_aaaa_aaaa_aaaa,
    0x0000_4555_5555_5555,
]);
pub const R2: Fp = Fp([
    0x0001_9999_9999_9eb8,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0004_cccc_cccc_cccc,
    0x0000_1999_9999_9999,
]);
pub const NRES_C: [u64; NLIMBS] = [
    0x4_cccc_cccc_cf5c,
    0x1_9999_9999_9999,
    0x3_3333_3333_3333,
    0x6_6666_6666_6666,
    0x0_0ccc_cccc_cccc,
];

#[inline]
pub fn modmul(a: &[u64; NLIMBS], b: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p4 = P_TOP as u128;
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

#[inline]
pub fn modsqr(a: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p4 = P_TOP as u128;
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
fn modnsqr(a: &mut [u64; NLIMBS], n: u32) {
    for _ in 0..n {
        let t = *a;
        modsqr(&t, a);
    }
}
#[inline]
fn modmul_ip(a: [u64; NLIMBS], b: [u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    modmul(&a, &b, c);
}

/// Progenitor: z = w^((p-3)/4). Addition chain from modarith.
pub fn modpro(w: &[u64; NLIMBS], z: &mut [u64; NLIMBS]) {
    let x = *w;
    let mut t0 = [0u64; NLIMBS];
    let mut t1 = [0u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4;
    modsqr(&x, z);
    modmul_ip(x, *z, &mut t0);
    modsqr(&t0, z);
    modmul_ip(x, *z, z);
    modsqr(z, &mut t1);
    modsqr(&t1, &mut t3);
    modsqr(&t3, &mut t2);
    t4 = t2;
    modnsqr(&mut t4, 3);
    modmul_ip(t2, t4, &mut t2);
    t4 = t2;
    modnsqr(&mut t4, 6);
    modmul_ip(t2, t4, &mut t2);
    t4 = t2;
    modnsqr(&mut t4, 2);
    modmul_ip(t3, t4, &mut t3);
    modnsqr(&mut t3, 13);
    modmul_ip(t2, t3, &mut t2);
    t3 = t2;
    modnsqr(&mut t3, 27);
    modmul_ip(t2, t3, &mut t2);
    modmul_ip(*z, t2, z);
    t2 = *z;
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

/// Partial reduce a saturated 4×64-bit value modulo p using 5·2²⁴⁸ ≡ 1.
pub fn partial_reduce(out: &mut [u64; PR_WORDS], src: &[u64; PR_WORDS]) {
    let h = src[3] >> 56;
    let l = src[3] & 0x00FF_FFFF_FFFF_FFFF;
    let quo = (h * 0xCD) >> 10;
    let rem = h - 5 * quo;
    let mut cc = adc(0, src[0], quo, &mut out[0]);
    cc = adc(cc, src[1], 0, &mut out[1]);
    cc = adc(cc, src[2], 0, &mut out[2]);
    let _ = adc(cc, l, rem << 56, &mut out[3]);
}

#[inline]
fn adc(cc: u8, a: u64, b: u64, d: &mut u64) -> u8 {
    let t = (a as u128) + (b as u128) + (cc as u128);
    *d = t as u64;
    (t >> 64) as u8
}
