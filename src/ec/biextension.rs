//! Pairings via cubical/biextension arithmetic, plus 2-adic Fp² discrete log.
//! Ported from `src/ec/ref/lvlx/biextension.c`.
//!
//! Used only on the signing path (id2iso). Not required for verification.
//!
//! Note from C: `cubical_add` is off by a factor x4 from the strict cubical
//! formula. This cancels in the Weil pairing and in the reduced Tate pairing
//! over Fp² (final exponentiation), but would be wrong for Tate over Fp.

use crate::gf::{
    fp2_add, fp2_batched_inv, fp2_frob, fp2_inv, fp2_is_equal, fp2_is_one, fp2_is_zero, fp2_mul,
    fp2_mul_ip, fp2_select, fp2_sqr, fp2_sqr_ip, fp2_sub, Fp2,
};
use crate::mp::{mp_add, multiple_mp_shiftl, Digit};
use crate::precomp::{NWORDS_ORDER, P_COFACTOR_FOR_2F, TORSION_EVEN_POWER};

use super::jac::jac_add;
#[cfg(any(debug_assertions, test))]
use super::{ec_biscalar_mul, ec_dbl_iter};
use super::{
    ec_curve_normalize_a24, ec_is_equal, ec_is_zero, lift_basis_normalized, test_basis_order_twof,
    test_point_order_twof, EcBasis, EcCurve, EcPoint, JacPoint,
};

#[derive(Clone, Copy, Default, Debug)]
pub struct PairingParams {
    pub e: u32,
    pub p: EcPoint,
    pub q: EcPoint,
    pub pq: EcPoint,
    pub ixp: Fp2,
    pub ixq: Fp2,
    pub a24: EcPoint,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct PairingDlogDiffPoints {
    pub pmr: EcPoint,
    pub pms: EcPoint,
    pub rmq: EcPoint,
    pub smq: EcPoint,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct PairingDlogParams {
    pub e: u32,
    pub pq: EcBasis,
    pub rs: EcBasis,
    pub diff: PairingDlogDiffPoints,
    pub ixp: Fp2,
    pub ixq: Fp2,
    pub ixr: Fp2,
    pub ixs: Fp2,
    pub a24: EcPoint,
}

// ---------------------------------------------------------------------------
// Cubical-torsor primitives
// ---------------------------------------------------------------------------

/// Cubical differential addition: like `xadd` but for `(1 : ixPQ)`-normalised
/// difference. Cost: 3M + 2S + 3a + 3s.
#[inline]
fn cubical_add(r: &mut EcPoint, p: &EcPoint, q: &EcPoint, ix_pq: &Fp2) {
    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut t2, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    let s0 = t0;
    fp2_mul(&mut t0, &s0, &t3);
    let s1 = t1;
    fp2_mul(&mut t1, &s1, &t2);
    fp2_add(&mut t2, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    fp2_sqr(&mut r.z, &t3);
    let s2 = t2;
    fp2_sqr(&mut t2, &s2);
    fp2_mul(&mut r.x, ix_pq, &t2);
}

/// Cubical combined double-and-add: given P, Q, ix(P-Q), compute P+Q and [2]Q.
/// Cost: 6M + 4S + 4a + 4s. Assumes A24.z == 1.
#[inline]
fn cubical_dbladd(
    ppq: &mut EcPoint,
    qq: &mut EcPoint,
    p: &EcPoint,
    q: &EcPoint,
    ix_pq: &Fp2,
    a24: &EcPoint,
) {
    debug_assert!(fp2_is_one(&a24.z) != 0);

    let mut t0 = Fp2::default();
    let mut t1 = Fp2::default();
    let mut t2 = Fp2::default();
    let mut t3 = Fp2::default();

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut ppq.x, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    fp2_sqr(&mut t2, &ppq.x);
    fp2_sqr(&mut qq.z, &t3);
    let s0 = t0;
    fp2_mul(&mut t0, &s0, &t3);
    let s1 = t1;
    fp2_mul(&mut t1, &s1, &ppq.x);
    fp2_add(&mut ppq.x, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    fp2_sqr(&mut ppq.z, &t3);
    let sx = ppq.x;
    fp2_sqr(&mut ppq.x, &sx);
    let sx = ppq.x;
    fp2_mul(&mut ppq.x, ix_pq, &sx);
    fp2_sub(&mut t3, &t2, &qq.z);
    fp2_mul(&mut qq.x, &t2, &qq.z);
    fp2_mul(&mut t0, &t3, &a24.x);
    let sz = qq.z;
    let s0 = t0;
    fp2_add(&mut t0, &s0, &sz);
    fp2_mul(&mut qq.z, &t0, &t3);
}

#[inline]
fn cubical_add_ip(r: &mut EcPoint, q: &EcPoint, ix_pq: &Fp2) {
    let p = *r;
    cubical_add(r, &p, q, ix_pq);
}

#[inline]
fn cubical_dbladd_ip(ppq: &mut EcPoint, qq: &mut EcPoint, ix_pq: &Fp2, a24: &EcPoint) {
    let p = *ppq;
    let q = *qq;
    cubical_dbladd(ppq, qq, &p, &q, ix_pq, a24);
}

/// Iterative biextension doubling: compute (PnQ, nQ) = (P + 2^e·Q-ish, 2^e·Q)
/// via e steps of `cubical_dbladd`.
fn biext_ladder_2e(
    e: u32,
    pnq: &mut EcPoint,
    nq: &mut EcPoint,
    pq: &EcPoint,
    q: &EcPoint,
    ixp: &Fp2,
    a24: &EcPoint,
) {
    *pnq = *pq;
    *nq = *q;
    for _ in 0..e {
        cubical_dbladd_ip(pnq, nq, ixp, a24);
    }
}

/// Monodromy ratio (X:Z) using `(1,0)` as cubical lift of 0_E.
#[inline]
fn point_ratio(r: &mut EcPoint, pnq: &EcPoint, nq: &EcPoint, p: &EcPoint) {
    debug_assert!(ec_is_zero(nq) != 0);
    debug_assert!(ec_is_equal(pnq, p) != 0);
    fp2_mul(&mut r.x, &nq.x, &p.x);
    r.z = pnq.x;
}

/// Cubical translation of `p` by a 2-torsion point `t`, in constant time.
fn translate(p: &mut EcPoint, t: &EcPoint) {
    let mut px_new = Fp2::default();
    let mut pz_new = Fp2::default();
    {
        let mut t0 = Fp2::default();
        let mut t1 = Fp2::default();
        fp2_mul(&mut t0, &t.x, &p.x);
        fp2_mul(&mut t1, &t.z, &p.z);
        fp2_sub(&mut px_new, &t0, &t1);
        fp2_mul(&mut t0, &t.z, &p.x);
        fp2_mul(&mut t1, &t.x, &p.z);
        fp2_sub(&mut pz_new, &t0, &t1);
    }
    let ta_is_zero = fp2_is_zero(&t.x);
    let snx = px_new;
    fp2_select(&mut px_new, &snx, &p.z, ta_is_zero);
    let snz = pz_new;
    fp2_select(&mut pz_new, &snz, &p.x, ta_is_zero);

    let tb_is_zero = fp2_is_zero(&t.z);
    let snx = px_new;
    fp2_select(&mut px_new, &snx, &p.x, tb_is_zero);
    let snz = pz_new;
    fp2_select(&mut pz_new, &snz, &p.z, tb_is_zero);

    p.x = px_new;
    p.z = pz_new;
}

#[inline]
fn translate_self(p: &mut EcPoint) {
    let t = *p;
    translate(p, &t);
}

/// Compute the biextension monodromy g_{P,Q}^{2^e} (level 1) via cubical
/// arithmetic of P + 2^e·Q. `swap_pq` chooses between (P, ixP) and (Q, ixQ).
fn monodromy_i(r: &mut EcPoint, pd: &PairingParams, swap_pq: bool) {
    let (p, q, ixp) = if swap_pq {
        (pd.q, pd.p, pd.ixq)
    } else {
        (pd.p, pd.q, pd.ixp)
    };
    let mut pnq = EcPoint::default();
    let mut nq = EcPoint::default();
    biext_ladder_2e(pd.e - 1, &mut pnq, &mut nq, &pd.pq, &q, &ixp, &pd.a24);
    translate(&mut pnq, &nq);
    translate_self(&mut nq);
    point_ratio(r, &pnq, &nq, &p);
}

/// Normalize P, Q to (X/Z : 1) and store PZ/PX, QZ/QX.
fn cubical_normalization(pd: &mut PairingParams, p: &EcPoint, q: &EcPoint) {
    let mut t = [p.x, p.z, q.x, q.z];
    fp2_batched_inv(&mut t);
    fp2_mul(&mut pd.ixp, &p.z, &t[0]);
    fp2_mul(&mut pd.ixq, &q.z, &t[2]);
    fp2_mul(&mut pd.p.x, &p.x, &t[1]);
    fp2_mul(&mut pd.q.x, &q.x, &t[3]);
    pd.p.z = Fp2::ONE;
    pd.q.z = Fp2::ONE;
}

/// Weil pairing on already-normalized data: r = e_{2^e}(P, Q).
fn weil_n(r: &mut Fp2, pd: &PairingParams) {
    let mut r0 = EcPoint::default();
    let mut r1 = EcPoint::default();
    monodromy_i(&mut r0, pd, true);
    monodromy_i(&mut r1, pd, false);
    fp2_mul(r, &r0.x, &r1.z);
    fp2_inv(r);
    fp2_mul_ip(r, &r0.z);
    fp2_mul_ip(r, &r1.x);
}

// ---------------------------------------------------------------------------
// Public pairings
// ---------------------------------------------------------------------------

/// Weil pairing e_{2^e}(P, Q). `pq` should be x(P+Q) (or x(P-Q); the sign
/// flips the result). Crashes (div-by-0) if P or Q is (0:1).
pub fn weil(r: &mut Fp2, e: u32, p: &EcPoint, q: &EcPoint, pq: &EcPoint, curve: &mut EcCurve) {
    let mut pd = PairingParams {
        e,
        ..Default::default()
    };
    cubical_normalization(&mut pd, p, q);
    pd.pq = *pq;
    ec_curve_normalize_a24(curve);
    pd.a24 = curve.a24;
    weil_n(r, &pd);
}

/// Raise an Fp² element to the cofactor `(p+1) / 2^f` (the odd part of the
/// group order divided by 2).
pub fn clear_cofac(r: &mut Fp2, a: &Fp2) {
    let mut exp: Digit = P_COFACTOR_FOR_2F[0];
    exp >>= 1;
    let x = *a;
    *r = *a;
    while exp > 0 {
        fp2_sqr_ip(r);
        if exp & 1 != 0 {
            fp2_mul_ip(r, &x);
        }
        exp >>= 1;
    }
}

/// Reduced 2^e-Tate pairing t_{2^e}(P, Q)^((p²-1)/2^e).
pub fn reduced_tate(
    r: &mut Fp2,
    e: u32,
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    curve: &mut EcCurve,
) {
    let e_full = TORSION_EVEN_POWER as u32;
    let e_diff = e_full - e;
    let mut pd = PairingParams {
        e,
        ..Default::default()
    };
    cubical_normalization(&mut pd, p, q);
    pd.pq = *pq;
    ec_curve_normalize_a24(curve);
    pd.a24 = curve.a24;

    let mut rr = EcPoint::default();
    monodromy_i(&mut rr, &pd, true);

    // reduced Tate = -(R.z/R.x)^((p^2-1)/2^f); split ^(p-1) into Frobenius * inverse.
    let mut frob = Fp2::default();
    let tmp = rr.x;
    fp2_frob(&mut frob, &rr.x);
    let s = rr.z;
    fp2_mul(&mut rr.x, &s, &frob);
    fp2_frob(&mut frob, &s);
    fp2_mul(&mut rr.z, &tmp, &frob);
    fp2_inv(&mut rr.x);
    fp2_mul(r, &rr.x, &rr.z);

    let s = *r;
    clear_cofac(r, &s);
    for _ in 0..e_diff {
        fp2_sqr_ip(r);
    }
}

// ---------------------------------------------------------------------------
// Fp² 2-adic discrete log
// ---------------------------------------------------------------------------

/// Upper bound on the recursion stack: ⌈log₂(TORSION_EVEN_POWER)⌉ + 1 ≤ 11 at lvl5.
const DLOG_DEPTH: usize = 16;

fn fp2_dlog_2e_rec(
    a: &mut [Digit; NWORDS_ORDER],
    len: i64,
    pows_f: &mut [Fp2],
    pows_g: &mut [Fp2],
    stacklen: usize,
) -> bool {
    if len == 0 {
        *a = [0; NWORDS_ORDER];
        return true;
    }
    if len == 1 {
        if fp2_is_one(&pows_f[stacklen - 1]) != 0 {
            *a = [0; NWORDS_ORDER];
            for i in 0..stacklen - 1 {
                fp2_sqr_ip(&mut pows_g[i]);
            }
            return true;
        }
        if fp2_is_equal(&pows_f[stacklen - 1], &pows_g[stacklen - 1]) != 0 {
            *a = [0; NWORDS_ORDER];
            a[0] = 1;
            for i in 0..stacklen - 1 {
                let g = pows_g[i];
                let f = pows_f[i];
                fp2_mul(&mut pows_f[i], &f, &g);
                fp2_sqr(&mut pows_g[i], &g);
            }
            return true;
        }
        return false;
    }
    // C uses (long)(len * 0.5); equivalent to integer division for len >= 0.
    let right = len / 2;
    let left = len - right;
    pows_f[stacklen] = pows_f[stacklen - 1];
    pows_g[stacklen] = pows_g[stacklen - 1];
    for _ in 0..left {
        let sf = pows_f[stacklen];
        fp2_sqr(&mut pows_f[stacklen], &sf);
        let sg = pows_g[stacklen];
        fp2_sqr(&mut pows_g[stacklen], &sg);
    }
    let mut dlp1 = [0; NWORDS_ORDER];
    let mut dlp2 = [0; NWORDS_ORDER];
    if !fp2_dlog_2e_rec(&mut dlp1, right, pows_f, pows_g, stacklen + 1) {
        return false;
    }
    if !fp2_dlog_2e_rec(&mut dlp2, left, pows_f, pows_g, stacklen) {
        return false;
    }
    multiple_mp_shiftl(&mut dlp2, right as u32, NWORDS_ORDER);
    mp_add(a, &dlp2, &dlp1, NWORDS_ORDER);
    true
}

/// Solve f = g^scal over the 2^e-roots of unity in Fp², given `g_inverse`.
fn fp2_dlog_2e(scal: &mut [Digit; NWORDS_ORDER], f: &Fp2, g_inverse: &Fp2, e: i32) -> bool {
    let mut log: usize = 0;
    let mut len = e as i64;
    while len > 1 {
        len >>= 1;
        log += 1;
    }
    log += 1;
    debug_assert!(log <= DLOG_DEPTH);
    let mut pows_f = [Fp2::default(); DLOG_DEPTH];
    let mut pows_g = [Fp2::default(); DLOG_DEPTH];
    pows_f[0] = *f;
    pows_g[0] = *g_inverse;
    *scal = [0; NWORDS_ORDER];
    let ok = fp2_dlog_2e_rec(scal, e as i64, &mut pows_f[..log], &mut pows_g[..log], 1);
    debug_assert!(ok);
    ok
}

// ---------------------------------------------------------------------------
// Batched pairings + dlog for change-of-basis recovery
// ---------------------------------------------------------------------------

fn cubical_normalization_dlog(pd: &mut PairingDlogParams, curve: &mut EcCurve) {
    let mut t = [
        pd.pq.p.x,
        pd.pq.p.z,
        pd.pq.q.x,
        pd.pq.q.z,
        pd.pq.pmq.x,
        pd.pq.pmq.z,
        pd.rs.p.x,
        pd.rs.p.z,
        pd.rs.q.x,
        pd.rs.q.z,
        curve.c,
    ];
    fp2_batched_inv(&mut t);

    let pz = pd.pq.p.z;
    fp2_mul(&mut pd.ixp, &pz, &t[0]);
    let px = pd.pq.p.x;
    fp2_mul(&mut pd.pq.p.x, &px, &t[1]);
    pd.pq.p.z = Fp2::ONE;

    let qz = pd.pq.q.z;
    fp2_mul(&mut pd.ixq, &qz, &t[2]);
    let qx = pd.pq.q.x;
    fp2_mul(&mut pd.pq.q.x, &qx, &t[3]);
    pd.pq.q.z = Fp2::ONE;

    let pmqx = pd.pq.pmq.x;
    fp2_mul(&mut pd.pq.pmq.x, &pmqx, &t[5]);
    pd.pq.pmq.z = Fp2::ONE;

    let rz = pd.rs.p.z;
    fp2_mul(&mut pd.ixr, &rz, &t[6]);
    let rx = pd.rs.p.x;
    fp2_mul(&mut pd.rs.p.x, &rx, &t[7]);
    pd.rs.p.z = Fp2::ONE;

    let sz = pd.rs.q.z;
    fp2_mul(&mut pd.ixs, &sz, &t[8]);
    let sx = pd.rs.q.x;
    fp2_mul(&mut pd.rs.q.x, &sx, &t[9]);
    pd.rs.q.z = Fp2::ONE;

    let a = curve.a;
    fp2_mul(&mut curve.a, &a, &t[10]);
    curve.c = Fp2::ONE;
}

fn compute_difference_points(pd: &mut PairingDlogParams, curve: &EcCurve) {
    let mut xy_p = JacPoint::default();
    let mut xy_q = JacPoint::default();
    let mut xy_r = JacPoint::default();
    let mut xy_s = JacPoint::default();

    let mut pq = pd.pq;
    let mut rs = pd.rs;
    lift_basis_normalized(&mut xy_p, &mut xy_q, &mut pq, curve);
    lift_basis_normalized(&mut xy_r, &mut xy_s, &mut rs, curve);
    pd.pq = pq;
    pd.rs = rs;

    let diff = |a: &JacPoint, b: &JacPoint| -> EcPoint {
        let mut t = JacPoint::default();
        jac_add(&mut t, &-*b, a, curve);
        t.into()
    };
    pd.diff.pmr = diff(&xy_p, &xy_r);
    pd.diff.pms = diff(&xy_p, &xy_s);
    pd.diff.rmq = diff(&xy_r, &xy_q);
    pd.diff.smq = diff(&xy_s, &xy_q);
}

#[allow(clippy::too_many_arguments)]
fn weil_dlog(
    r1: &mut [Digit; NWORDS_ORDER],
    r2: &mut [Digit; NWORDS_ORDER],
    s1: &mut [Digit; NWORDS_ORDER],
    s2: &mut [Digit; NWORDS_ORDER],
    pd: &PairingDlogParams,
) {
    let mut np = pd.pq.p;
    let mut nq = pd.pq.q;
    let mut nr = pd.rs.p;
    let mut ns = pd.rs.q;
    let mut npq = pd.pq.pmq;
    let mut pnq = pd.pq.pmq;
    let mut npr = pd.diff.pmr;
    let mut nps = pd.diff.pms;
    let mut pnr = pd.diff.pmr;
    let mut pns = pd.diff.pms;
    let mut nrq = pd.diff.rmq;
    let mut nsq = pd.diff.smq;
    let mut rnq = pd.diff.rmq;
    let mut snq = pd.diff.smq;

    for _ in 0..pd.e - 1 {
        cubical_add_ip(&mut npq, &np, &pd.ixq);
        cubical_add_ip(&mut npr, &np, &pd.ixr);
        cubical_dbladd_ip(&mut nps, &mut np, &pd.ixs, &pd.a24);

        cubical_add_ip(&mut pnq, &nq, &pd.ixp);
        cubical_add_ip(&mut rnq, &nq, &pd.ixr);
        cubical_dbladd_ip(&mut snq, &mut nq, &pd.ixs, &pd.a24);

        cubical_add_ip(&mut pnr, &nr, &pd.ixp);
        cubical_dbladd_ip(&mut nrq, &mut nr, &pd.ixq, &pd.a24);

        cubical_add_ip(&mut pns, &ns, &pd.ixp);
        cubical_dbladd_ip(&mut nsq, &mut ns, &pd.ixq, &pd.a24);
    }

    translate(&mut npq, &np);
    translate(&mut npr, &np);
    translate(&mut nps, &np);
    translate(&mut pnq, &nq);
    translate(&mut rnq, &nq);
    translate(&mut snq, &nq);
    translate(&mut pnr, &nr);
    translate(&mut nrq, &nr);
    translate(&mut pns, &ns);
    translate(&mut nsq, &ns);

    translate_self(&mut np);
    translate_self(&mut nq);
    translate_self(&mut nr);
    translate_self(&mut ns);

    let mut t0 = EcPoint::default();
    let mut t1 = EcPoint::default();
    let mut w1 = [Fp2::default(); 5];
    let mut w2 = [Fp2::default(); 5];

    // e(P,Q) = w0 — swap w1/w2 here so w1[0] is already 1/w0 after batch inv.
    point_ratio(&mut t0, &npq, &np, &pd.pq.q);
    point_ratio(&mut t1, &pnq, &nq, &pd.pq.p);
    fp2_mul(&mut w2[0], &t0.x, &t1.z);
    fp2_mul(&mut w1[0], &t1.x, &t0.z);

    // e(P,R) = w0^r2
    point_ratio(&mut t0, &npr, &np, &pd.rs.p);
    point_ratio(&mut t1, &pnr, &nr, &pd.pq.p);
    fp2_mul(&mut w1[1], &t0.x, &t1.z);
    fp2_mul(&mut w2[1], &t1.x, &t0.z);

    // e(R,Q) = w0^r1
    point_ratio(&mut t0, &nrq, &nr, &pd.pq.q);
    point_ratio(&mut t1, &rnq, &nq, &pd.rs.p);
    fp2_mul(&mut w1[2], &t0.x, &t1.z);
    fp2_mul(&mut w2[2], &t1.x, &t0.z);

    // e(P,S) = w0^s2
    point_ratio(&mut t0, &nps, &np, &pd.rs.q);
    point_ratio(&mut t1, &pns, &ns, &pd.pq.p);
    fp2_mul(&mut w1[3], &t0.x, &t1.z);
    fp2_mul(&mut w2[3], &t1.x, &t0.z);

    // e(S,Q) = w0^s1
    point_ratio(&mut t0, &nsq, &ns, &pd.pq.q);
    point_ratio(&mut t1, &snq, &nq, &pd.rs.q);
    fp2_mul(&mut w1[4], &t0.x, &t1.z);
    fp2_mul(&mut w2[4], &t1.x, &t0.z);

    fp2_batched_inv(&mut w1);
    for i in 0..5 {
        fp2_mul_ip(&mut w1[i], &w2[i]);
    }

    fp2_dlog_2e(r2, &w1[1], &w1[0], pd.e as i32);
    fp2_dlog_2e(r1, &w1[2], &w1[0], pd.e as i32);
    fp2_dlog_2e(s2, &w1[3], &w1[0], pd.e as i32);
    fp2_dlog_2e(s1, &w1[4], &w1[0], pd.e as i32);
}

/// Given bases ⟨P,Q⟩ and ⟨R,S⟩ (all of order 2^e), compute scalars with
/// R = [r1]P + [r2]Q, S = [s1]P + [s2]Q via batched Weil pairings + 2-adic dlog.
#[allow(clippy::too_many_arguments)]
pub fn ec_dlog_2_weil(
    r1: &mut [Digit; NWORDS_ORDER],
    r2: &mut [Digit; NWORDS_ORDER],
    s1: &mut [Digit; NWORDS_ORDER],
    s2: &mut [Digit; NWORDS_ORDER],
    pq: &EcBasis,
    rs: &EcBasis,
    curve: &mut EcCurve,
    e: i32,
) {
    debug_assert!(test_point_order_twof(&pq.q, curve, e));

    ec_curve_normalize_a24(curve);
    let mut pd = PairingDlogParams {
        e: e as u32,
        pq: *pq,
        rs: *rs,
        a24: curve.a24,
        ..Default::default()
    };
    cubical_normalization_dlog(&mut pd, curve);
    compute_difference_points(&mut pd, curve);
    weil_dlog(r1, r2, s1, s2, &pd);

    #[cfg(debug_assertions)]
    {
        let mut test;
        test = ec_biscalar_mul(r1, r2, e, pq, curve).unwrap();
        debug_assert!(ec_is_equal(&test, &rs.p) != 0);
        test = ec_biscalar_mul(s1, s2, e, pq, curve).unwrap();
        debug_assert!(ec_is_equal(&test, &rs.q) != 0);
    }
}

#[allow(clippy::too_many_arguments)]
fn tate_dlog_partial(
    r1: &mut [Digit; NWORDS_ORDER],
    r2: &mut [Digit; NWORDS_ORDER],
    s1: &mut [Digit; NWORDS_ORDER],
    s2: &mut [Digit; NWORDS_ORDER],
    pd: &PairingDlogParams,
) {
    let e_full = TORSION_EVEN_POWER as u32;
    let e_diff = e_full - pd.e;

    let mut np = pd.pq.p;
    let mut nq = pd.pq.q;
    let mut nr = pd.rs.p;
    let mut ns = pd.rs.q;
    let mut npq = pd.pq.pmq;
    let mut pnr = pd.diff.pmr;
    let mut pns = pd.diff.pms;
    let mut nrq = pd.diff.rmq;
    let mut nsq = pd.diff.smq;

    for _ in 0..e_full - 1 {
        cubical_dbladd_ip(&mut npq, &mut np, &pd.ixq, &pd.a24);
    }
    for _ in 0..pd.e - 1 {
        cubical_add_ip(&mut pnr, &nr, &pd.ixp);
        cubical_dbladd_ip(&mut nrq, &mut nr, &pd.ixq, &pd.a24);
        cubical_add_ip(&mut pns, &ns, &pd.ixp);
        cubical_dbladd_ip(&mut nsq, &mut ns, &pd.ixq, &pd.a24);
    }

    translate(&mut npq, &np);
    translate(&mut pnr, &nr);
    translate(&mut nrq, &nr);
    translate(&mut pns, &ns);
    translate(&mut nsq, &ns);

    translate_self(&mut np);
    translate_self(&mut nq);
    translate_self(&mut nr);
    translate_self(&mut ns);

    let mut t0 = EcPoint::default();
    let mut w1 = [Fp2::default(); 5];
    let mut w2 = [Fp2::default(); 5];

    point_ratio(&mut t0, &npq, &np, &pd.pq.q);
    w1[0] = t0.x;
    w2[0] = t0.z;

    point_ratio(&mut t0, &pnr, &nr, &pd.pq.p);
    w1[1] = t0.x;
    w2[1] = t0.z;

    point_ratio(&mut t0, &nrq, &nr, &pd.pq.q);
    w2[2] = t0.x;
    w1[2] = t0.z;

    point_ratio(&mut t0, &pns, &ns, &pd.pq.p);
    w1[3] = t0.x;
    w2[3] = t0.z;

    point_ratio(&mut t0, &nsq, &ns, &pd.pq.q);
    w2[4] = t0.x;
    w1[4] = t0.z;

    for i in 0..5 {
        let mut frob = Fp2::default();
        let tmp = w1[i];
        fp2_frob(&mut frob, &w1[i]);
        let s = w2[i];
        fp2_mul(&mut w1[i], &s, &frob);
        fp2_frob(&mut frob, &s);
        fp2_mul(&mut w2[i], &tmp, &frob);
    }

    fp2_batched_inv(&mut w2);
    for i in 0..5 {
        fp2_mul_ip(&mut w1[i], &w2[i]);
    }

    for i in 0..5 {
        let s = w1[i];
        clear_cofac(&mut w1[i], &s);
        for _ in 0..e_diff {
            fp2_sqr_ip(&mut w1[i]);
        }
    }

    fp2_dlog_2e(r2, &w1[1], &w1[0], pd.e as i32);
    fp2_dlog_2e(r1, &w1[2], &w1[0], pd.e as i32);
    fp2_dlog_2e(s2, &w1[3], &w1[0], pd.e as i32);
    fp2_dlog_2e(s1, &w1[4], &w1[0], pd.e as i32);
}

/// Tate-pairing variant of `ec_dlog_2_weil`. ⟨P,Q⟩ must be a basis for the
/// full 2^f-torsion; ⟨R,S⟩ may be of smaller order 2^e.
#[allow(clippy::too_many_arguments)]
pub fn ec_dlog_2_tate(
    r1: &mut [Digit; NWORDS_ORDER],
    r2: &mut [Digit; NWORDS_ORDER],
    s1: &mut [Digit; NWORDS_ORDER],
    s2: &mut [Digit; NWORDS_ORDER],
    pq: &EcBasis,
    rs: &EcBasis,
    curve: &mut EcCurve,
    e: i32,
) {
    debug_assert!(test_basis_order_twof(pq, curve, TORSION_EVEN_POWER as i32));

    ec_curve_normalize_a24(curve);
    let mut pd = PairingDlogParams {
        e: e as u32,
        pq: *pq,
        rs: *rs,
        a24: curve.a24,
        ..Default::default()
    };
    cubical_normalization_dlog(&mut pd, curve);
    compute_difference_points(&mut pd, curve);
    tate_dlog_partial(r1, r2, s1, s2, &pd);

    #[cfg(debug_assertions)]
    {
        let e_full = TORSION_EVEN_POWER as i32;
        let e_diff = e_full - e;
        let mut test;
        test = ec_biscalar_mul(r1, r2, e, pq, curve).unwrap();
        let t = test;
        ec_dbl_iter(&mut test, e_diff, &t, curve);
        debug_assert!(ec_is_equal(&test, &rs.p) != 0);
        test = ec_biscalar_mul(s1, s2, e, pq, curve).unwrap();
        let t = test;
        ec_dbl_iter(&mut test, e_diff, &t, curve);
        debug_assert!(ec_is_equal(&test, &rs.q) != 0);
    }
}

// ---------------------------------------------------------------------------
// Tests (ported from biextension-test.c)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ctrdrbg::DrbgState;
    use crate::ec::{ec_curve_to_basis_2f_to_hint, ec_is_zero, xadd, xdbl_a24};
    use crate::mp::mp_invert_matrix;

    fn fp2_exp_2e(r: &mut Fp2, e: u32, x: &Fp2) {
        *r = *x;
        for _ in 0..e {
            fp2_sqr_ip(r);
        }
    }

    fn setup_e6() -> (EcCurve, EcBasis, u32) {
        let mut curve = EcCurve::e0();
        curve.a = Fp2::from_small(6);
        curve.c = Fp2::ONE;
        ec_curve_normalize_a24(&mut curve);
        let e = TORSION_EVEN_POWER as u32;
        let mut basis = EcBasis::default();
        let _hint = ec_curve_to_basis_2f_to_hint(&mut basis, &mut curve, e as i32);
        (curve, basis, e)
    }

    #[test]
    fn weil_order_and_bilinearity() {
        let (mut curve, basis, e) = setup_e6();
        let p = basis.p;
        let q = basis.q;
        let pmq = basis.pmq;

        // [2^e]P = [2^e]Q = [2^e](P-Q) = 0
        for pt in [&p, &q, &pmq] {
            let mut tmp = EcPoint::default();
            ec_dbl_iter(&mut tmp, e as i32, pt, &mut curve);
            assert!(ec_is_zero(&tmp) != 0);
        }

        let mut ppq = EcPoint::default();
        xadd(&mut ppq, &p, &q, &pmq);

        let mut r1 = Fp2::default();
        weil(&mut r1, e, &p, &q, &ppq, &mut curve);

        // e(P,Q) has exact order 2^e in μ_{2^e}.
        let one = Fp2::ONE;
        let mut r2 = Fp2::default();
        fp2_exp_2e(&mut r2, e - 1, &r1);
        assert!(fp2_is_equal(&r2, &one) == 0);
        fp2_exp_2e(&mut r2, e, &r1);
        assert!(fp2_is_equal(&r2, &one) != 0);

        // Antisymmetry: e(P,Q,P+Q) = 1/e(P,Q,P-Q).
        let mut r2 = Fp2::default();
        weil(&mut r2, e, &p, &q, &pmq, &mut curve);
        fp2_inv(&mut r2);
        assert!(fp2_is_equal(&r1, &r2) != 0);

        // Bilinearity: e([2]P,Q) = e(P,[2]Q) = e(P,Q)^2.
        let a24 = curve.a24;
        let mut pp = EcPoint::default();
        let mut qq = EcPoint::default();
        let mut ppq2 = EcPoint::default();
        let mut pqq = EcPoint::default();
        xdbl_a24(&mut pp, &p, &a24, false);
        xdbl_a24(&mut qq, &q, &a24, false);
        xadd(&mut ppq2, &ppq, &p, &q);
        xadd(&mut pqq, &ppq, &q, &p);

        let mut r3 = Fp2::default();
        weil(&mut r2, e, &pp, &q, &ppq2, &mut curve);
        weil(&mut r3, e, &p, &qq, &pqq, &mut curve);
        assert!(fp2_is_equal(&r2, &r3) != 0);
        let mut rr1 = Fp2::default();
        fp2_sqr(&mut rr1, &r1);
        assert!(fp2_is_equal(&rr1, &r2) != 0);

        // Cubic case.
        let mut ppp = EcPoint::default();
        let mut qqq = EcPoint::default();
        let mut pppq = EcPoint::default();
        let mut pqqq = EcPoint::default();
        xadd(&mut ppp, &pp, &p, &p);
        xadd(&mut qqq, &qq, &q, &q);
        xadd(&mut pppq, &ppq2, &p, &ppq);
        xadd(&mut pqqq, &pqq, &q, &ppq);
        weil(&mut r2, e, &ppp, &q, &pppq, &mut curve);
        weil(&mut r3, e, &p, &qqq, &pqqq, &mut curve);
        assert!(fp2_is_equal(&r2, &r3) != 0);
        let mut rrr1 = Fp2::default();
        fp2_mul(&mut rrr1, &rr1, &r1);
        assert!(fp2_is_equal(&rrr1, &r2) != 0);
    }

    #[test]
    fn reduced_tate_order() {
        let (mut curve, basis, e) = setup_e6();
        let mut ppq = EcPoint::default();
        xadd(&mut ppq, &basis.p, &basis.q, &basis.pmq);
        let mut tp = Fp2::default();
        reduced_tate(&mut tp, e, &basis.p, &basis.q, &ppq, &mut curve);
        let one = Fp2::ONE;
        let mut r2 = Fp2::default();
        fp2_exp_2e(&mut r2, e - 1, &tp);
        assert!(fp2_is_equal(&r2, &one) == 0);
        fp2_exp_2e(&mut r2, e, &tp);
        assert!(fp2_is_equal(&r2, &one) != 0);
    }

    fn random_mixed_basis(
        curve: &mut EcCurve,
        bpq: &EcBasis,
        e: u32,
    ) -> (
        EcBasis,
        [Digit; NWORDS_ORDER],
        [Digit; NWORDS_ORDER],
        [Digit; NWORDS_ORDER],
        [Digit; NWORDS_ORDER],
    ) {
        let mut drbg = DrbgState::new(&[0u8; 48], None);
        let mut d1 = [0; NWORDS_ORDER];
        let mut d2 = [0; NWORDS_ORDER];
        let mut s1 = [0; NWORDS_ORDER];
        let mut s2 = [0; NWORDS_ORDER];
        for buf in [&mut d1, &mut d2, &mut s1, &mut s2] {
            let mut bytes = [0u8; (NWORDS_ORDER - 1) * 8];
            drbg.fill(&mut bytes);
            for (i, w) in buf[..NWORDS_ORDER - 1].iter_mut().enumerate() {
                *w = u64::from_le_bytes(bytes[i * 8..i * 8 + 8].try_into().unwrap());
            }
        }
        // s1, s2 odd; d1 even; d2 odd → r1 odd, r2 even → det odd.
        s1[0] |= 1;
        s2[0] |= 1;
        d1[0] &= !1;
        d2[0] |= 1;
        let mut r1 = [0; NWORDS_ORDER];
        let mut r2 = [0; NWORDS_ORDER];
        mp_add(&mut r1, &d1, &s1, NWORDS_ORDER);
        mp_add(&mut r2, &d2, &s2, NWORDS_ORDER);

        let mut brs = EcBasis::default();
        brs.p = ec_biscalar_mul(&r1, &r2, e as i32, bpq, curve).unwrap();
        brs.q = ec_biscalar_mul(&s1, &s2, e as i32, bpq, curve).unwrap();
        brs.pmq = ec_biscalar_mul(&d1, &d2, e as i32, bpq, curve).unwrap();
        (brs, r1, r2, s1, s2)
    }

    #[test]
    fn dlog_weil_and_tate() {
        let (mut curve, bpq, e) = setup_e6();
        let (brs, _r1, _r2, _s1, _s2) = random_mixed_basis(&mut curve, &bpq, e);

        let mut sr1 = [0; NWORDS_ORDER];
        let mut sr2 = [0; NWORDS_ORDER];
        let mut ss1 = [0; NWORDS_ORDER];
        let mut ss2 = [0; NWORDS_ORDER];

        // Weil-based dlog (debug_asserts inside check biscalar reconstruction).
        ec_dlog_2_weil(
            &mut sr1, &mut sr2, &mut ss1, &mut ss2, &bpq, &brs, &mut curve, e as i32,
        );
        let mut tmp = ec_biscalar_mul(&sr1, &sr2, e as i32, &bpq, &curve).unwrap();
        assert!(ec_is_equal(&tmp, &brs.p) != 0);
        tmp = ec_biscalar_mul(&ss1, &ss2, e as i32, &bpq, &curve).unwrap();
        assert!(ec_is_equal(&tmp, &brs.q) != 0);

        // Tate-based dlog at full order.
        ec_dlog_2_tate(
            &mut sr1, &mut sr2, &mut ss1, &mut ss2, &bpq, &brs, &mut curve, e as i32,
        );
        tmp = ec_biscalar_mul(&sr1, &sr2, e as i32, &bpq, &curve).unwrap();
        assert!(ec_is_equal(&tmp, &brs.p) != 0);
        tmp = ec_biscalar_mul(&ss1, &ss2, e as i32, &bpq, &curve).unwrap();
        assert!(ec_is_equal(&tmp, &brs.q) != 0);

        // Partial-torsion Tate.
        let e_full = TORSION_EVEN_POWER as i32;
        let e_partial: i32 = 126;
        let mut brs_p = brs;
        for pt in [&mut brs_p.p, &mut brs_p.q, &mut brs_p.pmq] {
            let s = *pt;
            ec_dbl_iter(pt, e_full - e_partial, &s, &mut curve);
        }
        ec_dlog_2_tate(
            &mut sr1, &mut sr2, &mut ss1, &mut ss2, &bpq, &brs_p, &mut curve, e_partial,
        );
        tmp = ec_biscalar_mul(&sr1, &sr2, e as i32, &bpq, &curve).unwrap();
        let s = tmp;
        ec_dbl_iter(&mut tmp, e_full - e_partial, &s, &mut curve);
        assert!(ec_is_equal(&tmp, &brs_p.p) != 0);
        tmp = ec_biscalar_mul(&ss1, &ss2, e as i32, &bpq, &curve).unwrap();
        let s = tmp;
        ec_dbl_iter(&mut tmp, e_full - e_partial, &s, &mut curve);
        assert!(ec_is_equal(&tmp, &brs_p.q) != 0);

        // Tate "to full basis": invert the recovered matrix.
        ec_dlog_2_tate(
            &mut sr1, &mut sr2, &mut ss1, &mut ss2, &bpq, &brs_p, &mut curve, e_partial,
        );
        mp_invert_matrix(
            &mut sr1,
            &mut sr2,
            &mut ss1,
            &mut ss2,
            e_partial,
            NWORDS_ORDER,
        );
        tmp = ec_biscalar_mul(&sr1, &sr2, e as i32, &brs_p, &curve).unwrap();
        let mut tmp2 = EcPoint::default();
        ec_dbl_iter(&mut tmp2, e_full - e_partial, &bpq.p, &mut curve);
        assert!(ec_is_equal(&tmp, &tmp2) != 0);
        tmp = ec_biscalar_mul(&ss1, &ss2, e as i32, &brs_p, &curve).unwrap();
        ec_dbl_iter(&mut tmp2, e_full - e_partial, &bpq.q, &mut curve);
        assert!(ec_is_equal(&tmp, &tmp2) != 0);
    }
}
