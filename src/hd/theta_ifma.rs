//! AVX-512 IFMA batched kernels for `ThetaPoint`.
//!
//! `ThetaSoa` keeps a `ThetaPoint`'s eight `Fp` coordinates in radix-2⁵² SoA
//! across an entire `double_iter` / `theta_isogeny_eval` so the IFMA pipeline
//! is not interrupted by per-step pack/unpack. Values stay in the asm `Fp`'s
//! Montgomery domain (R₆₄ = 2²⁵⁶) by post-shifting every IFMA product by 4
//! (since `fp_mul_ifma8` reduces wrt R₅₂ = 2²⁶⁰ = 16·R₆₄), so the boundary
//! conversion is bit-shuffle only.
//!
//! Lane order is `[x.re, y.re, z.re, t.re, x.im, y.im, z.im, t.im]`: the
//! Hadamard transform acts identically on the re-half and im-half, and each
//! half is a self-contained 4-lane butterfly.

#![cfg(all(
    feature = "lvl1",
    not(feature = "lvl3"),
    not(feature = "lvl5"),
    not(feature = "gf-portable"),
    target_arch = "x86_64",
    target_feature = "avx512ifma",
    target_feature = "avx512f"
))]
#![allow(unsafe_code)]

use core::arch::x86_64::*;

use super::ThetaPoint;
use crate::gf::fp_lvl1_ifma::{fp_mul8_r64, normalize_signed8, pack8, permute8, unpack8, FpIfma8};
#[cfg(test)]
use crate::gf::fp_lvl1_ifma::{fp_mul8_u64x4, fp_mul8_u64x4_shared_a};
#[cfg(test)]
use crate::gf::{Fp, FpInner};

// ThetaPoint memory order is [x.re, x.im, y.re, y.im, z.re, z.im, t.re, t.im];
// these tables convert to/from the SoA lane order [re×4 | im×4].
const DEINTERLEAVE: [u64; 8] = [0, 2, 4, 6, 1, 3, 5, 7];
const INTERLEAVE: [u64; 8] = [0, 4, 1, 5, 2, 6, 3, 7];
const SWAP_HALVES: [u64; 8] = [4, 5, 6, 7, 0, 1, 2, 3];

#[inline(always)]
fn tp_as_u64x4x8(p: &ThetaPoint) -> &[[u64; 4]; 8] {
    const _: () = assert!(core::mem::size_of::<ThetaPoint>() == 256);
    // SAFETY: `ThetaPoint` is `repr(C)` of 4 `Fp2`, each `repr(C)` of 2 `Fp`,
    // each `repr(transparent)` over `GF5_248` `repr(C)` over `[u64; 4]`.
    unsafe { &*core::ptr::from_ref(p).cast() }
}
#[inline(always)]
fn tp_from_u64x4x8(w: [[u64; 4]; 8]) -> ThetaPoint {
    unsafe { core::mem::transmute(w) }
}

#[cfg(test)]
#[inline(always)]
fn limbs(x: &Fp) -> [u64; 4] {
    x.0 .0
}
#[cfg(test)]
#[inline(always)]
fn from_limbs(w: [u64; 4]) -> Fp {
    Fp(FpInner(w))
}

/// One `ThetaPoint` in radix-2⁵² SoA (lanes `[x.re,y.re,z.re,t.re,x.im,…,t.im]`).
#[allow(missing_debug_implementations)]
#[derive(Clone, Copy)]
pub struct ThetaSoa(pub FpIfma8);

impl ThetaSoa {
    #[inline]
    pub fn from_tp(p: &ThetaPoint) -> Self {
        Self(permute8(&pack8(tp_as_u64x4x8(p)), &DEINTERLEAVE))
    }

    #[inline]
    pub fn to_tp(mut self) -> ThetaPoint {
        normalize_signed8(&mut self.0);
        tp_from_u64x4x8(unpack8(permute8(&self.0, &INTERLEAVE)))
    }

    /// Coordinate-wise Fp2 squaring, one IFMA batch.
    /// `(re+im)(re−im)` for the re-half, `(2re)·im` for the im-half.
    #[inline]
    pub fn pointwise_square_ip(&mut self) {
        unsafe {
            let v = &self.0.limbs;
            let s = permute8(&self.0, &SWAP_HALVES); // [im×4 | re×4]
            let mut a = FpIfma8::default(); // [re+im ×4 | 2re ×4]
            let mut b = FpIfma8::default(); // [re−im ×4 | im ×4]
            for i in 0..5 {
                a.limbs[i] = _mm512_add_epi64(v[i], s.limbs[i]); // re+im | im+re
                                                                 // im-half of `a` should be 2re; replace lanes 4..7 with re+re.
                a.limbs[i] = blend_hi(a.limbs[i], _mm512_add_epi64(s.limbs[i], s.limbs[i]));
                let diff = _mm512_sub_epi64(v[i], s.limbs[i]); // re−im (signed) | im−re
                b.limbs[i] = blend_hi(diff, v[i]); // [re−im | im]
            }
            normalize_signed8(&mut a);
            normalize_signed8(&mut b);
            self.0 = fp_mul8_r64(&a, &b);
        }
    }

    /// `self[k] *= c[k]` for each Fp2 coordinate, two IFMA batches.
    /// `c` must be strict (limbs < 2⁵²); self may be lazy.
    #[inline]
    pub fn pointwise_mul_ip(&mut self, c: &ThetaSoa) {
        unsafe {
            // self = [a×4|b×4], c = [c×4|d×4].
            //   r1 = self·c            → [ac×4 | bd×4]
            //   r2 = self·swap(c)      → [ad×4 | bc×4]
            //   re = ac − bd, im = ad + bc.
            normalize_signed8(&mut self.0);
            let r1 = fp_mul8_r64(&self.0, &c.0);
            let r2 = fp_mul8_r64(&self.0, &permute8(&c.0, &SWAP_HALVES));
            let r1s = permute8(&r1, &SWAP_HALVES); // [bd | ac]
            let r2s = permute8(&r2, &SWAP_HALVES); // [bc | ad]
            let mut out = FpIfma8::default();
            for i in 0..5 {
                let re = _mm512_sub_epi64(r1.limbs[i], r1s.limbs[i]); // ac−bd | bd−ac
                let im = _mm512_add_epi64(r2.limbs[i], r2s.limbs[i]); // ad+bc | bc+ad
                out.limbs[i] = blend_hi(re, im);
            }
            self.0 = out;
        }
    }

    /// Hadamard transform in SoA (signed radix-2⁵² adds; output lazy).
    #[inline]
    pub fn hadamard_ip(&mut self) {
        unsafe {
            // Lanes [x,y,z,t | X,Y,Z,T]. Two butterfly levels per limb:
            //   v ← [x+y, x−y, z+t, z−t | …]  via swap-within-pairs
            //   v ← [·+·, ·+·, ·−·, ·−· | …]  via swap-pairs
            static P1: [u64; 8] = [1, 0, 3, 2, 5, 4, 7, 6];
            static P2: [u64; 8] = [2, 3, 0, 1, 6, 7, 4, 5];
            let i1 = _mm512_loadu_si512(P1.as_ptr().cast());
            let i2 = _mm512_loadu_si512(P2.as_ptr().cast());
            for i in 0..5 {
                let v = self.0.limbs[i];
                let p1 = _mm512_permutexvar_epi64(i1, v);
                // lanes 0,2,4,6 = v+p1; lanes 1,3,5,7 = p1−v (= partner − self).
                let l1 = _mm512_mask_sub_epi64(_mm512_add_epi64(v, p1), 0xAA, p1, v);
                let p2 = _mm512_permutexvar_epi64(i2, l1);
                // lanes 0,1,4,5 = l1+p2; lanes 2,3,6,7 = p2−l1.
                self.0.limbs[i] = _mm512_mask_sub_epi64(_mm512_add_epi64(l1, p2), 0xCC, p2, l1);
            }
        }
    }

    #[inline]
    #[allow(clippy::wrong_self_convention)]
    pub fn to_squared_theta_ip(&mut self) {
        self.pointwise_square_ip();
        self.hadamard_ip();
    }
}

/// SoA `theta_isogeny_eval`. `precomp` is `phi.precomputation` already in SoA
/// (the chain driver converts it once per step).
#[inline]
pub fn theta_isogeny_eval_soa(
    out: &mut ThetaPoint,
    p: &ThetaPoint,
    precomp: &ThetaSoa,
    h1: bool,
    h2: bool,
) {
    let mut s = ThetaSoa::from_tp(p);
    if h1 {
        s.hadamard_ip();
    }
    s.to_squared_theta_ip();
    s.pointwise_mul_ip(precomp);
    if h2 {
        s.hadamard_ip();
    }
    *out = s.to_tp();
}

#[inline(always)]
unsafe fn blend_hi(lo_half: __m512i, hi_half: __m512i) -> __m512i {
    _mm512_mask_blend_epi64(0xF0, lo_half, hi_half)
}

// ---------------------------------------------------------------------------
// Per-call kernels (conversion at the boundary). Kept for tests/benches; the
// SoA-resident path above is what `hd` actually uses.
// ---------------------------------------------------------------------------

/// `p[k].square_ip()` for each Fp2 coordinate, in one IFMA batch.
#[cfg(test)]
fn pointwise_square_ifma(p: &mut ThetaPoint) {
    let q = [&p.x, &p.y, &p.z, &p.t];
    let mut a = [[0u64; 4]; 8];
    let mut b = [[0u64; 4]; 8];
    for k in 0..4 {
        let re = q[k].re;
        let im = q[k].im;
        a[k] = limbs(&Fp::add_noreduce(&re, &im));
        b[k] = limbs(&Fp::sub_2p_noreduce(&re, &im));
        a[k + 4] = limbs(&Fp::add_noreduce(&re, &re));
        b[k + 4] = limbs(&im);
    }
    let r = fp_mul8_u64x4(&a, &b);
    p.x.re = from_limbs(r[0]);
    p.y.re = from_limbs(r[1]);
    p.z.re = from_limbs(r[2]);
    p.t.re = from_limbs(r[3]);
    p.x.im = from_limbs(r[4]);
    p.y.im = from_limbs(r[5]);
    p.z.im = from_limbs(r[6]);
    p.t.im = from_limbs(r[7]);
}

/// `p[k] *= c[k]` for each Fp2 coordinate, in two IFMA batches.
#[cfg(test)]
fn pointwise_mul_ifma(p: &mut ThetaPoint, c: &ThetaPoint) {
    let pa = [&p.x, &p.y, &p.z, &p.t];
    let pc = [&c.x, &c.y, &c.z, &c.t];
    let mut a = [[0u64; 4]; 8];
    let mut b1 = [[0u64; 4]; 8];
    let mut b2 = [[0u64; 4]; 8];
    for k in 0..4 {
        a[k] = limbs(&pa[k].re);
        a[k + 4] = limbs(&pa[k].im);
        b1[k] = limbs(&pc[k].re);
        b1[k + 4] = limbs(&pc[k].im);
        b2[k] = limbs(&pc[k].im);
        b2[k + 4] = limbs(&pc[k].re);
    }
    let (r1, r2) = fp_mul8_u64x4_shared_a(&a, &b1, &b2);
    for k in 0..4 {
        let ac = from_limbs(r1[k]);
        let bd = from_limbs(r1[k + 4]);
        let ad = from_limbs(r2[k]);
        let bc = from_limbs(r2[k + 4]);
        let re = ac - bd;
        let im = ad + bc;
        match k {
            0 => p.x = crate::gf::Fp2 { re, im },
            1 => p.y = crate::gf::Fp2 { re, im },
            2 => p.z = crate::gf::Fp2 { re, im },
            _ => p.t = crate::gf::Fp2 { re, im },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{hadamard_ip, ThetaPoint, ThetaStructure};
    use super::*;
    use crate::test_util::Prng;

    fn rand_tp(prng: &mut Prng) -> ThetaPoint {
        ThetaPoint {
            x: prng.fp2(),
            y: prng.fp2(),
            z: prng.fp2(),
            t: prng.fp2(),
        }
    }

    fn assert_tp_eq(a: &ThetaPoint, b: &ThetaPoint) {
        assert_eq!(a.x, b.x);
        assert_eq!(a.y, b.y);
        assert_eq!(a.z, b.z);
        assert_eq!(a.t, b.t);
    }

    #[test]
    fn soa_roundtrip() {
        let mut prng = Prng(0x71F40);
        for _ in 0..2_000 {
            let p = rand_tp(&mut prng);
            assert_tp_eq(&ThetaSoa::from_tp(&p).to_tp(), &p);
        }
    }

    #[test]
    fn soa_pointwise_square() {
        let mut prng = Prng(0x71F4A);
        for _ in 0..5_000 {
            let p = rand_tp(&mut prng);
            let mut s = ThetaSoa::from_tp(&p);
            s.pointwise_square_ip();
            let got = s.to_tp();
            let want = ThetaPoint {
                x: p.x.square(),
                y: p.y.square(),
                z: p.z.square(),
                t: p.t.square(),
            };
            assert_tp_eq(&got, &want);
        }
    }

    #[test]
    fn soa_pointwise_mul() {
        let mut prng = Prng(0x71F4B);
        for _ in 0..5_000 {
            let p = rand_tp(&mut prng);
            let c = rand_tp(&mut prng);
            let mut s = ThetaSoa::from_tp(&p);
            s.pointwise_mul_ip(&ThetaSoa::from_tp(&c));
            let got = s.to_tp();
            let want = ThetaPoint {
                x: p.x * c.x,
                y: p.y * c.y,
                z: p.z * c.z,
                t: p.t * c.t,
            };
            assert_tp_eq(&got, &want);
        }
    }

    #[test]
    fn soa_hadamard() {
        let mut prng = Prng(0x71F4C);
        for _ in 0..5_000 {
            let p = rand_tp(&mut prng);
            let mut s = ThetaSoa::from_tp(&p);
            s.hadamard_ip();
            let got = s.to_tp();
            let mut want = p;
            hadamard_ip(&mut want);
            assert_tp_eq(&got, &want);
        }
    }

    #[test]
    fn soa_long_chain() {
        // Stress the lazy-limb invariant across many iterations.
        use super::super::theta_precomputation;
        let mut prng = Prng(0x71F4E);
        for _ in 0..200 {
            let np = rand_tp(&mut prng);
            let mut a = ThetaStructure {
                null_point: np,
                ..Default::default()
            };
            theta_precomputation(&mut a);
            let pd = ThetaSoa::from_tp(&a.precomp_d());
            let pp = ThetaSoa::from_tp(&a.precomp());
            let p = rand_tp(&mut prng);
            let mut s = ThetaSoa::from_tp(&p);
            for _ in 0..60 {
                s.to_squared_theta_ip();
                s.pointwise_square_ip();
                s.pointwise_mul_ip(&pd);
                s.hadamard_ip();
                s.pointwise_mul_ip(&pp);
            }
            let got = s.to_tp();
            let mut want = ThetaPoint::default();
            super::super::double_iter_scalar_ref(&mut want, &mut a, &p, 60);
            assert_tp_eq(&got, &want);
        }
    }

    #[test]
    fn soa_chain_matches_scalar_double_point() {
        use super::super::theta_precomputation;
        let mut prng = Prng(0x71F4D);
        for _ in 0..1_000 {
            let np = rand_tp(&mut prng);
            let mut a = ThetaStructure {
                null_point: np,
                ..Default::default()
            };
            theta_precomputation(&mut a);
            let pd = ThetaSoa::from_tp(&a.precomp_d());
            let pp = ThetaSoa::from_tp(&a.precomp());
            let p = rand_tp(&mut prng);
            // SoA: 5 doublings.
            let mut s = ThetaSoa::from_tp(&p);
            for _ in 0..5 {
                s.to_squared_theta_ip();
                s.pointwise_square_ip();
                s.pointwise_mul_ip(&pd);
                s.hadamard_ip();
                s.pointwise_mul_ip(&pp);
            }
            let got = s.to_tp();
            // Scalar reference.
            let mut want = ThetaPoint::default();
            super::super::double_iter_scalar_ref(&mut want, &mut a, &p, 5);
            assert_tp_eq(&got, &want);
        }
    }

    #[test]
    fn pointwise_square_matches_scalar() {
        let mut prng = Prng(0x71F4A);
        for _ in 0..2_000 {
            let p = rand_tp(&mut prng);
            let mut a = p;
            pointwise_square_ifma(&mut a);
            let want = ThetaPoint {
                x: p.x.square(),
                y: p.y.square(),
                z: p.z.square(),
                t: p.t.square(),
            };
            assert_tp_eq(&a, &want);
        }
    }

    #[test]
    fn pointwise_mul_matches_scalar() {
        let mut prng = Prng(0x71F4B);
        for _ in 0..2_000 {
            let p = rand_tp(&mut prng);
            let c = rand_tp(&mut prng);
            let mut a = p;
            pointwise_mul_ifma(&mut a, &c);
            let want = ThetaPoint {
                x: p.x * c.x,
                y: p.y * c.y,
                z: p.z * c.z,
                t: p.t * c.t,
            };
            assert_tp_eq(&a, &want);
        }
    }
}
