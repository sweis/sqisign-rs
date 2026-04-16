//! Wall-clock benchmarks for keygen/sign/verify on KAT vector 0.
//! Reports the minimum of N independent timings to reduce noise from
//! concurrent system load. Reads KAT-0 from the level-appropriate .rsp file.

use sqisign_rs::nistapi::crypto_sign_open;
use sqisign_rs::params::{CRYPTO_ALGNAME, CRYPTO_BYTES, CRYPTO_SECRETKEYBYTES};
use std::time::{Duration, Instant};

fn kat0() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let path = format!(
        "{}/KAT/PQCsignKAT_{}_{}.rsp",
        env!("CARGO_MANIFEST_DIR"),
        CRYPTO_SECRETKEYBYTES,
        CRYPTO_ALGNAME
    );
    let s = std::fs::read_to_string(path).expect("KAT file");
    let field = |k: &str| {
        let v = s.lines().find_map(|l| l.strip_prefix(k)).unwrap().trim();
        hex::decode(v).unwrap()
    };
    (field("seed = "), field("pk = "), field("sm = "))
}

fn time_n<F: FnMut()>(label: &str, n: usize, mut f: F) {
    f();
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
    let (seed, pk, sm) = kat0();

    time_n("verify", 200, || {
        let r = crypto_sign_open(&sm, &pk);
        assert!(r.is_ok());
    });

    #[cfg(feature = "sign")]
    {
        use sqisign_rs::common::ctrdrbg;
        use sqisign_rs::nistapi::{crypto_sign, crypto_sign_keypair};

        time_n("keygen", 50, || {
            ctrdrbg::randombytes_init(&seed);
            let _ = crypto_sign_keypair().unwrap();
        });

        ctrdrbg::randombytes_init(&seed);
        let (_, sk) = crypto_sign_keypair().unwrap();
        let post_kg = ctrdrbg::snapshot();
        let msg = &sm[CRYPTO_BYTES..];
        time_n("sign", 50, || {
            ctrdrbg::restore(&post_kg);
            let _ = crypto_sign(msg, &sk).unwrap();
        });
    }

    let _ = seed;
}
