// SPDX-License-Identifier: Apache-2.0
//! Hand-rolled rationals `Ibq = [Ibz; 2]` (num/denom), ported from `lll/rationals.c`.

use super::intbig::*;
use super::types::*;

pub type IbqVec4 = [Ibq; 4];
pub type IbqMat4x4 = [[Ibq; 4]; 4];

#[inline]
pub fn ibq_init() -> Ibq {
    [Ibz::default(), Ibz::from(1)]
}
#[inline]
pub fn ibq_vec_4_init() -> IbqVec4 {
    core::array::from_fn(|_| ibq_init())
}
#[inline]
pub fn ibq_mat_4x4_init() -> IbqMat4x4 {
    core::array::from_fn(|_| core::array::from_fn(|_| ibq_init()))
}

pub fn ibq_reduce(x: &mut Ibq) {
    let mut gcd = Ibz::default();
    let mut r = Ibz::default();
    ibz_gcd(&mut gcd, &x[0], &x[1]);
    let n = x[0].clone();
    ibz_div(&mut x[0], &mut r, &n, &gcd);
    debug_assert!(ibz_is_zero(&r) != 0);
    let d = x[1].clone();
    ibz_div(&mut x[1], &mut r, &d, &gcd);
    debug_assert!(ibz_is_zero(&r) != 0);
}

pub fn ibq_add(sum: &mut Ibq, a: &Ibq, b: &Ibq) {
    let mut add = Ibz::default();
    let mut prod = Ibz::default();
    ibz_mul(&mut add, &a[0], &b[1]);
    ibz_mul(&mut prod, &b[0], &a[1]);
    ibz_add(&mut sum[0], &add, &prod);
    ibz_mul(&mut sum[1], &a[1], &b[1]);
}

pub fn ibq_neg(neg: &mut Ibq, x: &Ibq) {
    ibz_copy(&mut neg[1], &x[1]);
    ibz_neg(&mut neg[0], &x[0]);
}

pub fn ibq_sub(diff: &mut Ibq, a: &Ibq, b: &Ibq) {
    let mut neg = ibq_init();
    ibq_neg(&mut neg, b);
    ibq_add(diff, a, &neg);
}

pub fn ibq_abs(abs: &mut Ibq, x: &Ibq) {
    let mut neg = ibq_init();
    ibq_neg(&mut neg, x);
    if ibq_cmp(x, &neg) < 0 {
        ibq_copy(abs, &neg);
    } else {
        ibq_copy(abs, x);
    }
}

pub fn ibq_mul(prod: &mut Ibq, a: &Ibq, b: &Ibq) {
    ibz_mul(&mut prod[0], &a[0], &b[0]);
    ibz_mul(&mut prod[1], &a[1], &b[1]);
}

pub fn ibq_inv(inv: &mut Ibq, x: &Ibq) -> i32 {
    let res = (ibq_is_zero(x) == 0) as i32;
    if res != 0 {
        ibz_copy(&mut inv[0], &x[0]);
        ibz_copy(&mut inv[1], &x[1]);
        inv.swap(0, 1);
    }
    res
}

pub fn ibq_cmp(a: &Ibq, b: &Ibq) -> i32 {
    let mut x = Ibz::default();
    let mut y = Ibz::default();
    ibz_mul(&mut x, &a[0], &b[1]);
    ibz_mul(&mut y, &b[0], &a[1]);
    // Negating both operands flips the comparison; track the sign instead.
    let mut sign = 1;
    if ibz_cmp(&a[1], ibz_const_zero()) > 0 {
        sign = -sign;
    }
    if ibz_cmp(&b[1], ibz_const_zero()) > 0 {
        sign = -sign;
    }
    sign * ibz_cmp(&x, &y)
}

pub fn ibq_is_zero(x: &Ibq) -> i32 {
    ibz_is_zero(&x[0])
}

pub fn ibq_is_one(x: &Ibq) -> i32 {
    (ibz_cmp(&x[0], &x[1]) == 0) as i32
}

pub fn ibq_set(q: &mut Ibq, a: &Ibz, b: &Ibz) -> i32 {
    ibz_copy(&mut q[0], a);
    ibz_copy(&mut q[1], b);
    (ibz_is_zero(b) == 0) as i32
}

pub fn ibq_copy(target: &mut Ibq, value: &Ibq) {
    ibz_copy(&mut target[0], &value[0]);
    ibz_copy(&mut target[1], &value[1]);
}

pub fn ibq_is_ibz(q: &Ibq) -> i32 {
    let mut r = Ibz::default();
    ibz_mod(&mut r, &q[0], &q[1]);
    ibz_is_zero(&r)
}

pub fn ibq_to_ibz(z: &mut Ibz, q: &Ibq) -> i32 {
    let mut r = Ibz::default();
    ibz_div(z, &mut r, &q[0], &q[1]);
    ibz_is_zero(&r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quaternion::intbig::ibz_from_i64;

    fn z(v: i64) -> Ibz {
        ibz_from_i64(v)
    }
    fn q(n: i64, d: i64) -> Ibq {
        [z(n), z(d)]
    }

    #[test]
    fn add_sub_mul() {
        let mut s = ibq_init();
        ibq_add(&mut s, &q(1, 2), &q(1, 3));
        ibq_reduce(&mut s);
        assert_eq!(s, [z(5), z(6)]);
        ibq_sub(&mut s, &q(1, 2), &q(1, 3));
        ibq_reduce(&mut s);
        assert_eq!(s, [z(1), z(6)]);
        ibq_mul(&mut s, &q(2, 3), &q(3, 5));
        ibq_reduce(&mut s);
        assert_eq!(s, [z(2), z(5)]);
    }

    #[test]
    fn cmp_handles_signs() {
        assert!(ibq_cmp(&q(1, 2), &q(1, 3)) > 0);
        assert!(ibq_cmp(&q(-1, 2), &q(1, 3)) < 0);
        assert!(ibq_cmp(&q(1, -2), &q(1, 3)) < 0);
        assert!(ibq_cmp(&q(-1, -2), &q(1, 3)) > 0);
        assert_eq!(ibq_cmp(&q(2, 4), &q(1, 2)), 0);
    }

    #[test]
    fn inv_abs_ibz() {
        let mut i = ibq_init();
        assert_eq!(ibq_inv(&mut i, &q(3, 7)), 1);
        assert_eq!(i, [z(7), z(3)]);
        assert_eq!(ibq_inv(&mut i, &q(0, 5)), 0);
        let mut a = ibq_init();
        ibq_abs(&mut a, &q(-3, 5));
        assert!(ibq_cmp(&a, &q(3, 5)) == 0);
        assert_eq!(ibq_is_ibz(&q(6, 3)), 1);
        assert_eq!(ibq_is_ibz(&q(5, 3)), 0);
        let mut zz = Ibz::default();
        assert_eq!(ibq_to_ibz(&mut zz, &q(6, 3)), 1);
        assert_eq!(zz, z(2));
    }
}
