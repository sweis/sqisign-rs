// SPDX-License-Identifier: Apache-2.0
//! Number-theoretic helpers, ported from `integers.c`.

use super::intbig::*;

/// Generate a random prime of `bitsize` bits, optionally ≡ 3 (mod 4).
pub fn ibz_generate_random_prime(
    p: &mut Ibz,
    is3mod4: bool,
    bitsize: u32,
    probability_test_iterations: i32,
) -> i32 {
    debug_assert!(bitsize != 0);
    let extra = is3mod4 as u32;
    let mut two_pow = Ibz::new();
    let mut two_powp = Ibz::new();
    ibz_pow(&mut two_pow, ibz_const_two(), bitsize - 1 - extra);
    ibz_pow(&mut two_powp, ibz_const_two(), bitsize - extra);

    let mut found = 0;
    let mut cnt: u64 = 0;
    while found == 0 {
        cnt += 1;
        if cnt % 100_000 == 0 {
            eprintln!(
                "Random prime generation is still running after {cnt} attempts, this is not normal! \
                 The expected number of attempts is {bitsize}"
            );
        }
        ibz_rand_interval(p, &two_pow, &two_powp);
        let t = p.clone();
        ibz_add(p, &t, &t);
        if is3mod4 {
            let t = p.clone();
            ibz_add(p, &t, &t);
            let t = p.clone();
            ibz_add(p, ibz_const_two(), &t);
        }
        let t = p.clone();
        ibz_add(p, ibz_const_one(), &t);
        found = ibz_probab_prime(p, probability_test_iterations);
    }
    found
}

/// Solve `x² + n·y² = p` for positive `x, y`. Assumes `p` prime.
/// Returns 1 on success, 0 if no solution exists.
pub fn ibz_cornacchia_prime(x: &mut Ibz, y: &mut Ibz, n: &Ibz, p: &Ibz) -> i32 {
    let mut r0 = Ibz::new();
    let mut r1 = Ibz::new();
    let mut r2 = Ibz::new();
    let mut a = Ibz::new();
    let mut prod = Ibz::new();

    // p = 2
    if ibz_cmp(p, ibz_const_two()) == 0 {
        if ibz_is_one(n) != 0 {
            ibz_set(x, 1);
            ibz_set(y, 1);
            return 1;
        }
        return 0;
    }
    // p = n
    if ibz_cmp(p, n) == 0 {
        ibz_set(x, 0);
        ibz_set(y, 1);
        return 1;
    }

    ibz_gcd(&mut r2, p, n);
    if ibz_is_one(&r2) == 0 {
        return 0;
    }

    // r2 ← √(-n) mod p
    ibz_neg(&mut r2, n);
    let nr2 = r2.clone();
    if ibz_sqrt_mod_p(&mut r2, &nr2, p) == 0 {
        return 0;
    }

    // Euclidean descent.
    ibz_copy(&mut prod, p);
    ibz_copy(&mut r1, p);
    ibz_copy(&mut r0, p);
    while ibz_cmp(&prod, p) >= 0 {
        ibz_div(&mut a, &mut r0, &r2, &r1);
        ibz_mul(&mut prod, &r0, &r0);
        ibz_copy(&mut r2, &r1);
        ibz_copy(&mut r1, &r0);
    }
    // Check (p - r0²) / n is a perfect square.
    ibz_sub(&mut a, p, &prod);
    let t = a.clone();
    ibz_div(&mut a, &mut r2, &t, n);
    if ibz_is_zero(&r2) == 0 {
        return 0;
    }
    if ibz_sqrt(y, &a) == 0 {
        return 0;
    }

    ibz_copy(x, &r0);
    ibz_mul(&mut a, y, y);
    let t = a.clone();
    ibz_mul(&mut a, &t, n);
    let t = prod.clone();
    ibz_add(&mut prod, &t, &a);
    (ibz_cmp(&prod, p) == 0) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;

    fn z(v: i64) -> Ibz {
        Integer::from(v)
    }

    #[test]
    fn cornacchia_small() {
        let (mut x, mut y) = (Ibz::new(), Ibz::new());
        // n=1: x² + y² = p
        for &(p, ex, ey) in &[(5, 1, 2), (13, 2, 3), (29, 2, 5), (97, 4, 9)] {
            assert_eq!(
                ibz_cornacchia_prime(&mut x, &mut y, &z(1), &z(p)),
                1,
                "p={p}"
            );
            let xx = x.clone() * x.clone() + y.clone() * y.clone();
            assert_eq!(xx, p);
            // Unordered pair.
            let (lo, hi) = if x < y {
                (x.clone(), y.clone())
            } else {
                (y.clone(), x.clone())
            };
            assert_eq!((lo, hi), (z(ex), z(ey)));
        }
        // n=1, p=2.
        assert_eq!(ibz_cornacchia_prime(&mut x, &mut y, &z(1), &z(2)), 1);
        assert_eq!((x.clone(), y.clone()), (z(1), z(1)));
        // n=1, p ≡ 3 mod 4: no solution.
        assert_eq!(ibz_cornacchia_prime(&mut x, &mut y, &z(1), &z(7)), 0);
        // n=3: x² + 3y² = 7 → (2,1).
        assert_eq!(ibz_cornacchia_prime(&mut x, &mut y, &z(3), &z(7)), 1);
        assert_eq!(x.clone() * x.clone() + z(3) * y.clone() * y.clone(), z(7));
        // p = n.
        assert_eq!(ibz_cornacchia_prime(&mut x, &mut y, &z(13), &z(13)), 1);
        assert_eq!((x, y), (z(0), z(1)));
    }

    #[test]
    fn cornacchia_large() {
        let p = Integer::from_str_radix("100000000000000000000000000000000000133", 10).unwrap();
        assert!(ibz_probab_prime(&p, 30) > 0);
        let (mut x, mut y) = (Ibz::new(), Ibz::new());
        // p ≡ 1 mod 4, so x²+y²=p solvable.
        assert_eq!(p.mod_u(4), 1);
        assert_eq!(ibz_cornacchia_prime(&mut x, &mut y, &z(1), &p), 1);
        assert_eq!(x.clone() * x.clone() + y.clone() * y.clone(), p);
    }
}
