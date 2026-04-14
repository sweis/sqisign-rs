// SPDX-License-Identifier: Apache-2.0
//! 2ⁿ-isogeny chain evaluation and curve isomorphisms.
//! Port of `ec/ref/lvlx/isog_chains.c`.

use super::*;
use crate::gf::*;

/// Strategy-based 2ⁿ-isogeny walk using degree-4 steps.
fn ec_eval_even_strategy(
    curve: &mut EcCurve,
    points: &mut [EcPoint],
    kernel: &EcPoint,
    isog_len: i32,
) -> u32 {
    ec_curve_normalize_a24(curve);
    let mut a24 = curve.a24;

    let mut space = 1usize;
    let mut i = 1;
    while i < isog_len {
        space += 1;
        i *= 2;
    }

    let mut splits: Vec<EcPoint> = vec![EcPoint::default(); space];
    let mut todo: Vec<u16> = vec![0; space];
    splits[0] = *kernel;
    todo[0] = isog_len as u16;

    let mut current: isize = 0;

    for j in 0..isog_len / 2 {
        debug_assert!(current >= 0);
        debug_assert!(todo[current as usize] >= 1);
        while todo[current as usize] != 2 {
            debug_assert!(todo[current as usize] >= 3);
            current += 1;
            debug_assert!((current as usize) < space);
            splits[current as usize] = splits[(current - 1) as usize];
            let prev = todo[(current - 1) as usize];
            let mut num_dbls = prev / 4 * 2 + prev % 2;
            todo[current as usize] = prev - num_dbls;
            while num_dbls > 0 {
                let s = splits[current as usize];
                xdbl_a24(&mut splits[current as usize], &s, &a24, false);
                num_dbls -= 1;
            }
        }

        if j == 0 {
            debug_assert!(fp2_is_one(&a24.z) != 0);
            if ec_is_four_torsion(&splits[current as usize], curve) == 0 {
                return u32::MAX;
            }
            let mut t = EcPoint::default();
            xdbl_a24(&mut t, &splits[current as usize], &a24, false);
            if fp2_is_zero(&t.x) != 0 {
                return u32::MAX;
            }
        } else {
            debug_assert_eq!(todo[current as usize], 2);
            #[cfg(debug_assertions)]
            {
                if fp2_is_zero(&splits[current as usize].z) != 0 {
                    eprintln!("splitting point z coordinate is unexpectedly zero");
                }
                let mut test = EcPoint::default();
                xdbl_a24(&mut test, &splits[current as usize], &a24, false);
                if fp2_is_zero(&test.z) != 0 {
                    eprintln!("z coordinate is unexpectedly zero before doubling");
                }
                let s = test;
                xdbl_a24(&mut test, &s, &a24, false);
                if fp2_is_zero(&test.z) == 0 {
                    eprintln!("z coordinate is unexpectedly not zero after doubling");
                }
            }
        }

        let mut kps4 = EcKps4::default();
        let kp = splits[current as usize];
        xisog_4(&mut kps4, &mut a24, &kp);
        xeval_4_inplace(&mut splits[..current as usize], &kps4);
        for i in 0..current as usize {
            todo[i] -= 2;
        }
        xeval_4_inplace(points, &kps4);

        current -= 1;
    }
    debug_assert!(if isog_len % 2 != 0 { current == 0 } else { current == -1 });

    if isog_len % 2 != 0 {
        #[cfg(debug_assertions)]
        {
            if fp2_is_zero(&splits[0].z) != 0 {
                eprintln!("splitting point z coordinate is unexpectedly zero");
            }
            let mut test = splits[0];
            let s = test;
            xdbl_a24(&mut test, &s, &a24, false);
            if fp2_is_zero(&test.z) == 0 {
                eprintln!("z coordinate is unexpectedly not zero after doubling");
            }
        }

        if isog_len == 1 && ec_is_two_torsion(&splits[0], curve) == 0 {
            return u32::MAX;
        }
        if fp2_is_zero(&splits[0].x) != 0 {
            return u32::MAX;
        }

        let mut kps2 = EcKps2::default();
        let kp = splits[0];
        xisog_2(&mut kps2, &mut a24, &kp);
        xeval_2_inplace(points, &kps2);
    }

    a24_to_ac(curve, &a24);
    curve.is_a24_computed_and_normalized = false;
    0
}

/// Evaluate a 2ⁿ-isogeny on a list of points; image curve written to `image`.
/// Returns 0 on success, 0xFFFFFFFF if the kernel is malformed.
pub fn ec_eval_even(image: &mut EcCurve, phi: &EcIsogEven, points: &mut [EcPoint]) -> u32 {
    copy_curve(image, &phi.curve);
    ec_eval_even_strategy(image, points, &phi.kernel, phi.length as i32)
}

/// Naive (multiplicative-strategy) 2ⁿ-isogeny chain. Returns 0 on success.
pub fn ec_eval_small_chain(
    curve: &mut EcCurve,
    kernel: &EcPoint,
    len: i32,
    points: &mut [EcPoint],
    special: bool,
) -> u32 {
    let mut a24 = EcPoint::default();
    ac_to_a24(&mut a24, curve);

    let mut kps = EcKps2::default();
    let mut big_k = *kernel;

    for i in 0..len {
        let mut small_k = big_k;
        for _ in 0..len - i - 1 {
            let s = small_k;
            xdbl_a24(&mut small_k, &s, &a24, false);
        }
        if i == 0 && ec_is_two_torsion(&small_k, curve) == 0 {
            return u32::MAX;
        }
        if fp2_is_zero(&small_k.x) != 0 {
            if special {
                let mut b24 = EcPoint::default();
                xisog_2_singular(&mut kps, &mut b24, a24);
                {
                    let bk = [big_k];
                    xeval_2_singular(core::slice::from_mut(&mut big_k), &bk, &kps);
                }
                xeval_2_singular_inplace(points, &kps);
                a24 = b24;
            } else {
                return u32::MAX;
            }
        } else {
            xisog_2(&mut kps, &mut a24, &small_k);
            {
                let bk = [big_k];
                xeval_2(core::slice::from_mut(&mut big_k), &bk, &kps);
            }
            xeval_2_inplace(points, &kps);
        }
    }
    a24_to_ac(curve, &a24);
    curve.is_a24_computed_and_normalized = false;
    0
}

/// Compute an isomorphism `from → to`. Returns 0xFFFFFFFF on degenerate output.
pub fn ec_isomorphism(isom: &mut EcIsom, from: &EcCurve, to: &EcCurve) -> u32 {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();
    let mut t4 = Fp2::default();

    fp2_mul(&mut t0, &from.a, &from.c);
    fp2_mul(&mut t1, &to.a, &to.c);

    fp2_mul(&mut t2, &t1, &to.c);
    fp2_add(&mut t3, &t2, &t2);
    fp2_dbl_ip(&mut t3);
    fp2_dbl_ip(&mut t3);
    fp2_add_ip(&mut t2, &t3);
    fp2_sqr(&mut t3, &to.a);
    fp2_mul_ip(&mut t3, &to.a);
    fp2_dbl_ip(&mut t3);
    fp2_sub(&mut isom.nx, &t3, &t2);
    fp2_mul(&mut t2, &t0, &from.a);
    fp2_sqr(&mut t3, &from.c);
    fp2_mul_ip(&mut t3, &from.c);
    fp2_add(&mut t4, &t3, &t3);
    fp2_add_ip(&mut t3, &t4);
    fp2_sub_ip(&mut t3, &t2);
    let nx = isom.nx;
    fp2_mul(&mut isom.nx, &nx, &t3);

    fp2_mul(&mut t2, &t0, &from.c);
    fp2_add(&mut t3, &t2, &t2);
    fp2_dbl_ip(&mut t3);
    fp2_dbl_ip(&mut t3);
    fp2_add_ip(&mut t2, &t3);
    fp2_sqr(&mut t3, &from.a);
    fp2_mul_ip(&mut t3, &from.a);
    fp2_dbl_ip(&mut t3);
    fp2_sub(&mut isom.d, &t3, &t2);
    fp2_mul(&mut t2, &t1, &to.a);
    fp2_sqr(&mut t3, &to.c);
    fp2_mul_ip(&mut t3, &to.c);
    fp2_add(&mut t4, &t3, &t3);
    fp2_add_ip(&mut t3, &t4);
    fp2_sub_ip(&mut t3, &t2);
    let d = isom.d;
    fp2_mul(&mut isom.d, &d, &t3);

    fp2_mul(&mut t0, &to.c, &from.a);
    fp2_mul_ip(&mut t0, &isom.nx);
    fp2_mul(&mut t1, &from.c, &to.a);
    fp2_mul_ip(&mut t1, &isom.d);
    fp2_sub(&mut isom.nz, &t0, &t1);
    fp2_mul(&mut t0, &from.c, &to.c);
    fp2_add(&mut t1, &t0, &t0);
    fp2_add_ip(&mut t0, &t1);
    let d = isom.d;
    fp2_mul(&mut isom.d, &d, &t0);
    let nx = isom.nx;
    fp2_mul(&mut isom.nx, &nx, &t0);

    fp2_is_zero(&isom.nx) | fp2_is_zero(&isom.d)
}

/// In-place evaluation of an isomorphism on a point.
pub fn ec_iso_eval(p: &mut EcPoint, isom: &EcIsom) {
    let mut tmp = Fp2::default();
    let px = p.x;
    fp2_mul(&mut p.x, &px, &isom.nx);
    fp2_mul(&mut tmp, &p.z, &isom.nz);
    let px = p.x;
    fp2_add(&mut p.x, &px, &tmp);
    let pz = p.z;
    fp2_mul(&mut p.z, &pz, &isom.d);
}
