//! Microbenchmark: AVX-512 IFMA `fp_mul_ifma` (radix-2⁵², R=2²⁶⁰) vs the
//! production scalar BMI2/ADX `fp_mul` (4-limb, R=2²⁵⁶).
//!
//! Run: RUSTFLAGS="-C target-cpu=native" cargo bench --no-default-features \
//!        --features lvl1 --bench bench_fp_ifma
#![allow(unused, unsafe_code)]

use std::hint::black_box;
use std::time::Instant;

#[cfg(not(all(
    feature = "lvl1",
    target_arch = "x86_64",
    target_feature = "avx512ifma"
)))]
fn main() {
    println!("(avx512ifma not enabled — set RUSTFLAGS=\"-C target-cpu=native\" on Ice Lake/SPR)");
}

#[cfg(all(
    feature = "lvl1",
    target_arch = "x86_64",
    target_feature = "avx512ifma"
))]
fn main() {
    use sqisign_rs::gf::fp::FpInner;
    use sqisign_rs::gf::fp_lvl1_ifma::{fp_mul_ifma, fp_mul_ref52, FpIfma, MASK52};

    fn time_n<F: FnMut()>(label: &str, n: u64, mut f: F) -> f64 {
        f();
        let mut best = f64::MAX;
        for _ in 0..7 {
            let t = Instant::now();
            for _ in 0..n {
                f();
            }
            let per = t.elapsed().as_secs_f64() / n as f64;
            best = best.min(per);
        }
        println!("{label:24}  {:>8.2} ns  ({n} iters, min of 7)", best * 1e9);
        best
    }

    // Seed values (Montgomery-domain in their respective representations).
    let x4 = FpInner::from(0x1234_5678_9ABC_DEF0u64) * FpInner::from(0xCAFE_BABE_DEAD_BEEFu64);
    let mut x5 = FpIfma::ZERO;
    let mut y5 = FpIfma::ZERO;
    let seed = 0x1234_5678_9ABC_DEF0u64;
    let mut s = seed;
    for i in 0..5 {
        s = s.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        x5.0[i] = s & MASK52;
        s = s.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        y5.0[i] = s & MASK52;
    }

    println!("== Single-element fp_mul (latency, dependent chain) ==");
    let mut acc = x4;
    let t_asm = time_n("fp_mul (BMI2/ADX asm)", 10_000_000, || {
        acc.set_mul(&black_box(x4));
        black_box(&acc);
    });
    let mut acc = x4;
    let t_asm_sq = time_n("fp_sqr (BMI2/ADX asm)", 10_000_000, || {
        acc.set_square();
        black_box(&acc);
    });
    let mut acc = x5;
    let t_ifma = time_n("fp_mul_ifma (AVX512)", 10_000_000, || {
        acc = fp_mul_ifma(&acc, &black_box(y5));
        black_box(&acc);
    });
    let mut acc = x5;
    let t_ref52 = time_n("fp_mul_ref52 (scalar)", 10_000_000, || {
        acc = fp_mul_ref52(&acc, &black_box(y5));
        black_box(&acc);
    });
    println!("  ifma/asm     = {:.3}x", t_ifma / t_asm);
    println!("  ref52/asm    = {:.3}x", t_ref52 / t_asm);

    println!("\n== Throughput (8 independent chains) ==");
    let mut a4 = [x4; 8];
    let t_asm8 = time_n("fp_mul asm ×8 indep", 2_000_000, || {
        for i in 0..8 {
            a4[i].set_mul(&black_box(x4));
        }
        black_box(&a4);
    });
    let mut a5 = [x5; 8];
    let t_ifma8 = time_n("fp_mul_ifma ×8 indep", 2_000_000, || {
        for i in 0..8 {
            a5[i] = fp_mul_ifma(&a5[i], &black_box(y5));
        }
        black_box(&a5);
    });
    println!(
        "  per-elem: asm {:.2} ns, ifma {:.2} ns ({:.3}x)",
        t_asm8 / 8.0 * 1e9,
        t_ifma8 / 8.0 * 1e9,
        t_ifma8 / t_asm8
    );

    println!("\n== Strategy B: batched 8-way IFMA (SoA) ==");
    use sqisign_rs::gf::fp_lvl1_ifma::{fp_mul_ifma8, FpIfma8};
    let a8 = FpIfma8::from_singles(&[x5; 8]);
    let b8 = FpIfma8::from_singles(&[y5; 8]);
    let mut acc8 = a8;
    let t_b8 = time_n("fp_mul_ifma8 (8 muls)", 5_000_000, || {
        acc8 = fp_mul_ifma8(&acc8, &black_box(b8));
        black_box(&acc8);
    });
    println!(
        "  per-elem: {:.2} ns  vs asm-throughput {:.2} ns  ({:.3}x)",
        t_b8 / 8.0 * 1e9,
        t_asm8 / 8.0 * 1e9,
        t_b8 / t_asm8
    );
    println!(
        "  vs asm-latency 1elem {:.2} ns: {:.3}x",
        t_asm * 1e9,
        (t_b8 / 8.0) / t_asm
    );
}
