// SPDX-License-Identifier: Apache-2.0
//! Wall-clock benchmarks for keygen/sign/verify on KAT vector 0.
//! Reports the minimum of N independent timings to reduce noise from
//! concurrent system load.
//! Run: cargo bench --features lvl1,sign

use sqisign_rs::nistapi::crypto_sign_open;
use std::time::{Duration, Instant};

const KAT0_PK: &str = "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B";
const KAT0_SM: &str = "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
const KAT0_SEED: &str = "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1";

/// Time `f` `n` times, report min/median/mean per call.
fn time_n<F: FnMut()>(label: &str, n: usize, mut f: F) {
    f(); // warm-up
    let mut samples: Vec<Duration> = (0..n)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed()
        })
        .collect();
    samples.sort();
    let min = samples[0];
    let med = samples[n / 2];
    let mean = samples.iter().sum::<Duration>() / n as u32;
    println!("{label:12}  min={min:>10.3?}  med={med:>10.3?}  mean={mean:>10.3?}  ({n} iters)");
}

fn main() {
    let pk = hex::decode(KAT0_PK).unwrap();
    let sm = hex::decode(KAT0_SM).unwrap();

    time_n("verify", 200, || {
        let r = crypto_sign_open(&sm, &pk);
        assert!(r.is_ok());
    });

    #[cfg(feature = "sign")]
    {
        use sqisign_rs::common::ctrdrbg;
        use sqisign_rs::nistapi::{crypto_sign, crypto_sign_keypair};
        let seed = hex::decode(KAT0_SEED).unwrap();

        time_n("keygen", 50, || {
            ctrdrbg::randombytes_init(&seed);
            let _ = crypto_sign_keypair().unwrap();
        });

        ctrdrbg::randombytes_init(&seed);
        let (_, sk) = crypto_sign_keypair().unwrap();
        let post_kg = ctrdrbg::snapshot();
        let msg = &sm[148..];
        time_n("sign", 50, || {
            ctrdrbg::restore(&post_kg);
            let _ = crypto_sign(msg, &sk).unwrap();
        });
    }
}
