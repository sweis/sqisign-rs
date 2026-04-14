// SPDX-License-Identifier: Apache-2.0
//! GF(p) backend for the lvl3 prime p = 65·2³⁷⁶ − 1.
//! Ported from `fp_p65376_64.c` (modarith / Montgomery, unsaturated radix-2⁵⁵).

#![allow(clippy::unreadable_literal)]

use super::fp::Fp;

pub const NWORDS_FIELD: usize = 7;
pub const NWORDS_ORDER: usize = 6;
pub const BITS: u32 = 384;
pub const LOG2P: u32 = 9;
pub const FP_ENCODED_BYTES: usize = 48;

pub const NLIMBS: usize = 7;
pub const LIMB_BITS: u32 = 55;
pub const NBYTES: usize = 48;
pub const MASK: u64 = (1u64 << LIMB_BITS) - 1;

pub const P_TOP: u64 = 0x10400000000000;
pub const TWO_P_TOP: u64 = 0x20800000000000;

pub const PR_WORDS: usize = 6;

pub const ZERO: Fp = Fp([0; NLIMBS]);
pub const ONE: Fp = Fp([0x0000000000000007, 0, 0, 0, 0, 0, 0x000e400000000000]);
pub const MINUS_ONE: Fp = Fp([
    0x7ffffffffffff8,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x1ffffffffffff,
]);
pub const TWO_INV: Fp = Fp([0x0000000000000003, 0, 0, 0, 0, 0, 0x000f400000000000]);
pub const THREE_INV: Fp = Fp([
    0x0055555555555557,
    0x002aaaaaaaaaaaaa,
    0x0055555555555555,
    0x002aaaaaaaaaaaaa,
    0x0055555555555555,
    0x002aaaaaaaaaaaaa,
    0x000f955555555555,
]);
pub const R2: Fp = Fp([
    0x0007e07e07e07e26,
    0x007c0fc0fc0fc0fc,
    0x0001f81f81f81f81,
    0x003f03f03f03f03f,
    0x00607e07e07e07e0,
    0x000fc0fc0fc0fc0f,
    0x000e9f81f81f81f8,
]);
pub const NRES_C: [u64; NLIMBS] = [
    0xfc0fc0fc0fc4d,
    0x781f81f81f81f8,
    0x3f03f03f03f03,
    0x7e07e07e07e07e,
    0x40fc0fc0fc0fc0,
    0x1f81f81f81f81f,
    0xcff03f03f03f0,
];

#[inline]
pub fn modmul(a: &[u64; NLIMBS], b: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p6 = P_TOP as u128;
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
    let v4 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[5] as u128);
    t += (a[1] as u128) * (b[4] as u128);
    t += (a[2] as u128) * (b[3] as u128);
    t += (a[3] as u128) * (b[2] as u128);
    t += (a[4] as u128) * (b[1] as u128);
    t += (a[5] as u128) * (b[0] as u128);
    let v5 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[6] as u128);
    t += (a[1] as u128) * (b[5] as u128);
    t += (a[2] as u128) * (b[4] as u128);
    t += (a[3] as u128) * (b[3] as u128);
    t += (a[4] as u128) * (b[2] as u128);
    t += (a[5] as u128) * (b[1] as u128);
    t += (a[6] as u128) * (b[0] as u128);
    t += (v0 as u128) * p6;
    let v6 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[1] as u128) * (b[6] as u128);
    t += (a[2] as u128) * (b[5] as u128);
    t += (a[3] as u128) * (b[4] as u128);
    t += (a[4] as u128) * (b[3] as u128);
    t += (a[5] as u128) * (b[2] as u128);
    t += (a[6] as u128) * (b[1] as u128);
    t += (v1 as u128) * p6;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[2] as u128) * (b[6] as u128);
    t += (a[3] as u128) * (b[5] as u128);
    t += (a[4] as u128) * (b[4] as u128);
    t += (a[5] as u128) * (b[3] as u128);
    t += (a[6] as u128) * (b[2] as u128);
    t += (v2 as u128) * p6;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[3] as u128) * (b[6] as u128);
    t += (a[4] as u128) * (b[5] as u128);
    t += (a[5] as u128) * (b[4] as u128);
    t += (a[6] as u128) * (b[3] as u128);
    t += (v3 as u128) * p6;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[4] as u128) * (b[6] as u128);
    t += (a[5] as u128) * (b[5] as u128);
    t += (a[6] as u128) * (b[4] as u128);
    t += (v4 as u128) * p6;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[5] as u128) * (b[6] as u128);
    t += (a[6] as u128) * (b[5] as u128);
    t += (v5 as u128) * p6;
    c[4] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[6] as u128) * (b[6] as u128);
    t += (v6 as u128) * p6;
    c[5] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[6] = t as u64;
}

#[inline]
pub fn modsqr(a: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p6 = P_TOP as u128;
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
    let v4 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[5] as u128);
    tot += (a[1] as u128) * (a[4] as u128);
    tot += (a[2] as u128) * (a[3] as u128);
    tot *= 2;
    t += tot;
    let v5 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[6] as u128);
    tot += (a[1] as u128) * (a[5] as u128);
    tot += (a[2] as u128) * (a[4] as u128);
    tot *= 2;
    tot += (a[3] as u128) * (a[3] as u128);
    t += tot;
    t += (v0 as u128) * p6;
    let v6 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[1] as u128) * (a[6] as u128);
    tot += (a[2] as u128) * (a[5] as u128);
    tot += (a[3] as u128) * (a[4] as u128);
    tot *= 2;
    t += tot;
    t += (v1 as u128) * p6;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[2] as u128) * (a[6] as u128);
    tot += (a[3] as u128) * (a[5] as u128);
    tot *= 2;
    tot += (a[4] as u128) * (a[4] as u128);
    t += tot;
    t += (v2 as u128) * p6;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[3] as u128) * (a[6] as u128);
    tot += (a[4] as u128) * (a[5] as u128);
    tot *= 2;
    t += tot;
    t += (v3 as u128) * p6;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[4] as u128) * (a[6] as u128);
    tot *= 2;
    tot += (a[5] as u128) * (a[5] as u128);
    t += tot;
    t += (v4 as u128) * p6;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[5] as u128) * (a[6] as u128);
    tot *= 2;
    t += tot;
    t += (v5 as u128) * p6;
    c[4] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[6] as u128) * (a[6] as u128);
    t += tot;
    t += (v6 as u128) * p6;
    c[5] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[6] = t as u64;
}

#[inline(always)]
fn modnsqr(a: &mut [u64; NLIMBS], n: u32) {
    for _ in 0..n {
        let t = *a;
        modsqr(&t, a);
    }
}
#[inline(always)]
fn modmul_ip(a: [u64; NLIMBS], b: [u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    modmul(&a, &b, c);
}
#[inline(always)]
fn modsqr_ip(a: [u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    modsqr(&a, c);
}

/// Progenitor: z = w^((p-3)/4). Addition chain from modarith.
pub fn modpro(w: &[u64; NLIMBS], z: &mut [u64; NLIMBS]) {
    let x = *w;
    let mut t0 = [0u64; NLIMBS];
    let mut t1 = [0u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4 = [0u64; NLIMBS];
    let mut t5;
    modsqr(&x, z);
    modsqr(z, &mut t0);
    modmul_ip(x, t0, &mut t1);
    modmul_ip(*z, t1, z);
    modsqr(z, &mut t0);
    modsqr(&t0, &mut t3);
    modsqr(&t3, &mut t4);
    modsqr(&t4, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 3);
    modmul_ip(t2, t5, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 6);
    modmul_ip(t2, t5, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 2);
    modmul_ip(t4, t5, &mut t5);
    modnsqr(&mut t5, 13);
    modmul_ip(t2, t5, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 2);
    modmul_ip(t4, t5, &mut t4);
    modnsqr(&mut t4, 28);
    modmul_ip(t2, t4, &mut t2);
    modsqr(&t2, &mut t4);
    modmul_ip(t3, t4, &mut t3);
    modnsqr(&mut t3, 59);
    modmul_ip(t2, t3, &mut t2);
    modmul_ip(t1, t2, &mut t1);
    modmul_ip(*z, t1, z);
    modmul_ip(t0, *z, &mut t0);
    modmul_ip(t1, t0, &mut t1);
    modsqr(&t1, &mut t2);
    modmul_ip(t1, t2, &mut t2);
    modsqr_ip(t2, &mut t2);
    modmul_ip(t1, t2, &mut t2);
    modmul_ip(t0, t2, &mut t0);
    modmul_ip(*z, t0, z);
    modsqr(z, &mut t2);
    modmul_ip(*z, t2, &mut t2);
    modmul_ip(t0, t2, &mut t0);
    modmul_ip(t1, t0, &mut t1);
    t2 = t1;
    modnsqr(&mut t2, 128);
    modmul_ip(t1, t2, &mut t1);
    modmul_ip(t0, t1, &mut t0);
    modnsqr(&mut t0, 125);
    modmul_ip(*z, t0, z);
}

/// Partial reduce a saturated 6×64-bit value modulo p using 65·2³⁷⁶ ≡ 1.
pub fn partial_reduce(out: &mut [u64; PR_WORDS], src: &[u64; PR_WORDS]) {
    let h = src[5] >> 56;
    let l = src[5] & 0x00FF_FFFF_FFFF_FFFF;
    let quo = (h * 0xFC1) >> 18;
    let rem = h - 65 * quo;
    let mut cc = adc(0, src[0], quo, &mut out[0]);
    cc = adc(cc, src[1], 0, &mut out[1]);
    cc = adc(cc, src[2], 0, &mut out[2]);
    cc = adc(cc, src[3], 0, &mut out[3]);
    cc = adc(cc, src[4], 0, &mut out[4]);
    let _ = adc(cc, l, rem << 56, &mut out[5]);
}

#[inline]
fn adc(cc: u8, a: u64, b: u64, d: &mut u64) -> u8 {
    let t = (a as u128) + (b as u128) + (cc as u128);
    *d = t as u64;
    (t >> 64) as u8
}
