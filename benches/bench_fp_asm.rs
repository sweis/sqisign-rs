//! Microbenchmark: `fp_sqr_asm` vs `fp_mul_asm(x,x)` and full verify under
//! the gf-fast backend.
//! Run: RUSTFLAGS="-C target-cpu=native" cargo bench --no-default-features \
//!        --features lvl1,gf-fast,bench-internals --bench bench_fp_asm

use sqisign_rs::gf::fp::FpInner;
use sqisign_rs::nistapi::crypto_sign_open;
use std::hint::black_box;
use std::time::Instant;

const KAT0_PK: &str = "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B";
const KAT0_SM: &str = "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";

fn time_n<F: FnMut()>(label: &str, n: u64, mut f: F) -> f64 {
    f();
    let mut best = f64::MAX;
    for _ in 0..5 {
        let t = Instant::now();
        for _ in 0..n {
            f();
        }
        let per = t.elapsed().as_secs_f64() / n as f64;
        best = best.min(per);
    }
    println!("{label:20}  {:>8.2} ns  ({n} iters, min of 5)", best * 1e9);
    best
}

fn main() {
    // Warm a non-trivial element.
    let x = FpInner::from(0x1234_5678_9ABC_DEF0u64) * FpInner::from(0xCAFE_BABE_DEAD_BEEFu64);

    #[cfg(gf5_248_asm)]
    {
        println!("== fp microbench (asm) ==");
        let mut acc = x;
        let t_sqr = time_n("fp_sqr_asm", 10_000_000, || {
            acc.set_square();
            black_box(&acc);
        });
        let mut acc = x;
        let t_mul = time_n("fp_mul_asm(x,x)", 10_000_000, || {
            acc.set_square_via_mul();
            black_box(&acc);
        });
        println!("  ratio sqr/mul = {:.3}", t_sqr / t_mul);
    }

    #[cfg(not(gf5_248_asm))]
    println!("(asm not enabled — set RUSTFLAGS=\"-C target-cpu=native\")");

    println!("\n== verify (KAT-0, 200×) ==");
    let pk = hex::decode(KAT0_PK).unwrap();
    let sm = hex::decode(KAT0_SM).unwrap();
    time_n("verify", 200, || {
        let r = crypto_sign_open(&sm, &pk);
        assert!(r.is_ok());
    });
}
