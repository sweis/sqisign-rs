//! GF(p) backend for the lvl5 prime p = 27·2⁵⁰⁰ − 1.
//! Ported from `fp_p27500_64.c` (modarith / Montgomery, unsaturated radix-2⁵⁷).

#![allow(clippy::unreadable_literal)]

use super::fp::Fp;

pub const NWORDS_FIELD: usize = 9;
pub const NWORDS_ORDER: usize = 8;
pub const BITS: u32 = 512;
pub const LOG2P: u32 = 9;
pub const FP_ENCODED_BYTES: usize = 64;

pub const NLIMBS: usize = 9;
pub const LIMB_BITS: u32 = 57;
pub const NBYTES: usize = 64;
pub const MASK: u64 = (1u64 << LIMB_BITS) - 1;

pub const P_TOP: u64 = 0x1b00000000000;
pub const TWO_P_TOP: u64 = 0x3600000000000;

pub const PR_WORDS: usize = 8;

pub const ZERO: Fp = Fp([0; NLIMBS]);
pub const ONE: Fp = Fp([0x000000000000012f, 0, 0, 0, 0, 0, 0, 0, 0x0000b00000000000]);
pub const MINUS_ONE: Fp = Fp([
    0x1fffffffffffed0,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0xffffffffffff,
]);
pub const TWO_INV: Fp = Fp([0x0000000000000097, 0, 0, 0, 0, 0, 0, 0, 0x0001300000000000]);
pub const THREE_INV: Fp = Fp([
    0x00aaaaaaaaaaab0f,
    0x0155555555555555,
    0x00aaaaaaaaaaaaaa,
    0x0155555555555555,
    0x00aaaaaaaaaaaaaa,
    0x0155555555555555,
    0x00aaaaaaaaaaaaaa,
    0x0155555555555555,
    0x00015aaaaaaaaaaa,
]);
pub const R2: Fp = Fp([
    0x0012f684bda1e334,
    0x01425ed097b425ed,
    0x01684bda12f684bd,
    0x01ed097b425ed097,
    0x00bda12f684bda12,
    0x0097b425ed097b42,
    0x0012f684bda12f68,
    0x01425ed097b425ed,
    0x00008bda12f684bd,
]);
pub const NRES_C: [u64; NLIMBS] = [
    0x25ed097b43c668,
    0x84bda12f684bda,
    0xd097b425ed097b,
    0x1da12f684bda12f,
    0x17b425ed097b425,
    0x12f684bda12f684,
    0x25ed097b425ed0,
    0x84bda12f684bda,
    0x117b425ed097b,
];

#[inline]
pub fn modmul(a: &[u64; NLIMBS], b: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p8 = P_TOP as u128;
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
    let v6 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[7] as u128);
    t += (a[1] as u128) * (b[6] as u128);
    t += (a[2] as u128) * (b[5] as u128);
    t += (a[3] as u128) * (b[4] as u128);
    t += (a[4] as u128) * (b[3] as u128);
    t += (a[5] as u128) * (b[2] as u128);
    t += (a[6] as u128) * (b[1] as u128);
    t += (a[7] as u128) * (b[0] as u128);
    let v7 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[0] as u128) * (b[8] as u128);
    t += (a[1] as u128) * (b[7] as u128);
    t += (a[2] as u128) * (b[6] as u128);
    t += (a[3] as u128) * (b[5] as u128);
    t += (a[4] as u128) * (b[4] as u128);
    t += (a[5] as u128) * (b[3] as u128);
    t += (a[6] as u128) * (b[2] as u128);
    t += (a[7] as u128) * (b[1] as u128);
    t += (a[8] as u128) * (b[0] as u128);
    t += (v0 as u128) * p8;
    let v8 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[1] as u128) * (b[8] as u128);
    t += (a[2] as u128) * (b[7] as u128);
    t += (a[3] as u128) * (b[6] as u128);
    t += (a[4] as u128) * (b[5] as u128);
    t += (a[5] as u128) * (b[4] as u128);
    t += (a[6] as u128) * (b[3] as u128);
    t += (a[7] as u128) * (b[2] as u128);
    t += (a[8] as u128) * (b[1] as u128);
    t += (v1 as u128) * p8;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[2] as u128) * (b[8] as u128);
    t += (a[3] as u128) * (b[7] as u128);
    t += (a[4] as u128) * (b[6] as u128);
    t += (a[5] as u128) * (b[5] as u128);
    t += (a[6] as u128) * (b[4] as u128);
    t += (a[7] as u128) * (b[3] as u128);
    t += (a[8] as u128) * (b[2] as u128);
    t += (v2 as u128) * p8;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[3] as u128) * (b[8] as u128);
    t += (a[4] as u128) * (b[7] as u128);
    t += (a[5] as u128) * (b[6] as u128);
    t += (a[6] as u128) * (b[5] as u128);
    t += (a[7] as u128) * (b[4] as u128);
    t += (a[8] as u128) * (b[3] as u128);
    t += (v3 as u128) * p8;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[4] as u128) * (b[8] as u128);
    t += (a[5] as u128) * (b[7] as u128);
    t += (a[6] as u128) * (b[6] as u128);
    t += (a[7] as u128) * (b[5] as u128);
    t += (a[8] as u128) * (b[4] as u128);
    t += (v4 as u128) * p8;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[5] as u128) * (b[8] as u128);
    t += (a[6] as u128) * (b[7] as u128);
    t += (a[7] as u128) * (b[6] as u128);
    t += (a[8] as u128) * (b[5] as u128);
    t += (v5 as u128) * p8;
    c[4] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[6] as u128) * (b[8] as u128);
    t += (a[7] as u128) * (b[7] as u128);
    t += (a[8] as u128) * (b[6] as u128);
    t += (v6 as u128) * p8;
    c[5] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[7] as u128) * (b[8] as u128);
    t += (a[8] as u128) * (b[7] as u128);
    t += (v7 as u128) * p8;
    c[6] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    t += (a[8] as u128) * (b[8] as u128);
    t += (v8 as u128) * p8;
    c[7] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[8] = t as u64;
}

#[inline]
pub fn modsqr(a: &[u64; NLIMBS], c: &mut [u64; NLIMBS]) {
    let p8 = P_TOP as u128;
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
    let v6 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[7] as u128);
    tot += (a[1] as u128) * (a[6] as u128);
    tot += (a[2] as u128) * (a[5] as u128);
    tot += (a[3] as u128) * (a[4] as u128);
    tot *= 2;
    t += tot;
    let v7 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[0] as u128) * (a[8] as u128);
    tot += (a[1] as u128) * (a[7] as u128);
    tot += (a[2] as u128) * (a[6] as u128);
    tot += (a[3] as u128) * (a[5] as u128);
    tot *= 2;
    tot += (a[4] as u128) * (a[4] as u128);
    t += tot;
    t += (v0 as u128) * p8;
    let v8 = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[1] as u128) * (a[8] as u128);
    tot += (a[2] as u128) * (a[7] as u128);
    tot += (a[3] as u128) * (a[6] as u128);
    tot += (a[4] as u128) * (a[5] as u128);
    tot *= 2;
    t += tot;
    t += (v1 as u128) * p8;
    c[0] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[2] as u128) * (a[8] as u128);
    tot += (a[3] as u128) * (a[7] as u128);
    tot += (a[4] as u128) * (a[6] as u128);
    tot *= 2;
    tot += (a[5] as u128) * (a[5] as u128);
    t += tot;
    t += (v2 as u128) * p8;
    c[1] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[3] as u128) * (a[8] as u128);
    tot += (a[4] as u128) * (a[7] as u128);
    tot += (a[5] as u128) * (a[6] as u128);
    tot *= 2;
    t += tot;
    t += (v3 as u128) * p8;
    c[2] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[4] as u128) * (a[8] as u128);
    tot += (a[5] as u128) * (a[7] as u128);
    tot *= 2;
    tot += (a[6] as u128) * (a[6] as u128);
    t += tot;
    t += (v4 as u128) * p8;
    c[3] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[5] as u128) * (a[8] as u128);
    tot += (a[6] as u128) * (a[7] as u128);
    tot *= 2;
    t += tot;
    t += (v5 as u128) * p8;
    c[4] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[6] as u128) * (a[8] as u128);
    tot *= 2;
    tot += (a[7] as u128) * (a[7] as u128);
    t += tot;
    t += (v6 as u128) * p8;
    c[5] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[7] as u128) * (a[8] as u128);
    tot *= 2;
    t += tot;
    t += (v7 as u128) * p8;
    c[6] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    tot = (a[8] as u128) * (a[8] as u128);
    t += tot;
    t += (v8 as u128) * p8;
    c[7] = (t as u64) & MASK;
    t >>= LIMB_BITS;
    c[8] = t as u64;
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
/// Progenitor: z = w^((p-3)/4). Addition chain from modarith.
pub fn modpro(w: &[u64; NLIMBS], z: &mut [u64; NLIMBS]) {
    let x = *w;
    let mut t0 = [0u64; NLIMBS];
    let mut t1 = [0u64; NLIMBS];
    let mut t2 = [0u64; NLIMBS];
    let mut t3 = [0u64; NLIMBS];
    let mut t4 = [0u64; NLIMBS];
    let mut t5 = [0u64; NLIMBS];
    let mut t6 = [0u64; NLIMBS];
    *z = x;
    modnsqr(z, 2);
    modmul_ip(x, *z, &mut t0);
    modmul_ip(x, t0, z);
    modsqr(z, &mut t1);
    modmul_ip(x, t1, &mut t1);
    modsqr(&t1, &mut t3);
    modsqr(&t3, &mut t2);
    modmul_ip(t3, t2, &mut t4);
    modsqr(&t4, &mut t5);
    t2 = t5;
    modnsqr(&mut t2, 2);
    modsqr(&t2, &mut t6);
    modmul_ip(t2, t6, &mut t6);
    modmul_ip(t5, t6, &mut t5);
    modnsqr(&mut t5, 5);
    modmul_ip(t2, t5, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 12);
    modmul_ip(t2, t5, &mut t2);
    t5 = t2;
    modnsqr(&mut t5, 2);
    modmul_ip(t2, t5, &mut t5);
    modmul_ip(t4, t5, &mut t4);
    modsqr(&t4, &mut t5);
    modmul_ip(t2, t5, &mut t2);
    modmul_ip(t4, t2, &mut t4);
    modnsqr(&mut t4, 27);
    modmul_ip(t2, t4, &mut t2);
    modmul_ip(t1, t2, &mut t2);
    t4 = t2;
    modnsqr(&mut t4, 2);
    modmul_ip(t3, t4, &mut t3);
    modnsqr(&mut t3, 58);
    modmul_ip(t2, t3, &mut t2);
    modmul_ip(*z, t2, z);
    t2 = *z;
    modnsqr(&mut t2, 4);
    modmul_ip(t1, t2, &mut t1);
    modmul_ip(t0, t1, &mut t0);
    modmul_ip(t1, t0, &mut t1);
    modsqr(&t1, &mut t2);
    modmul_ip(t0, t2, &mut t0);
    t2 = t0;
    modnsqr(&mut t2, 2);
    modmul_ip(t0, t2, &mut t2);
    modmul_ip(t1, t2, &mut t1);
    modmul_ip(t0, t1, &mut t0);
    modnsqr(&mut t1, 128);
    modmul_ip(t0, t1, &mut t1);
    modnsqr(&mut t1, 128);
    modmul_ip(t0, t1, &mut t0);
    modnsqr(&mut t0, 119);
    modmul_ip(*z, t0, z);
}

/// Partial reduce a saturated 8×64-bit value modulo p using 27·2⁵⁰⁰ ≡ 1.
pub fn partial_reduce(out: &mut [u64; PR_WORDS], src: &[u64; PR_WORDS]) {
    let h = src[7] >> 52;
    let l = src[7] & 0x000F_FFFF_FFFF_FFFF;
    let quo = (h * 0x12F7) >> 17;
    let rem = h - 27 * quo;
    let mut cc = adc(0, src[0], quo, &mut out[0]);
    cc = adc(cc, src[1], 0, &mut out[1]);
    cc = adc(cc, src[2], 0, &mut out[2]);
    cc = adc(cc, src[3], 0, &mut out[3]);
    cc = adc(cc, src[4], 0, &mut out[4]);
    cc = adc(cc, src[5], 0, &mut out[5]);
    cc = adc(cc, src[6], 0, &mut out[6]);
    let _ = adc(cc, l, rem << 52, &mut out[7]);
}

#[inline]
fn adc(cc: u8, a: u64, b: u64, d: &mut u64) -> u8 {
    let t = (a as u128) + (b as u128) + (cc as u128);
    *d = t as u64;
    (t >> 64) as u8
}
