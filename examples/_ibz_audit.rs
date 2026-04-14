// Run a few KAT keypair+sign cycles and report the max ibz width seen.
// RUSTFLAGS='--cfg ibz_audit' cargo run --release --features lvl1,sign --example _ibz_audit
fn main() {
    use sqisign_rs::common::ctrdrbg;
    use sqisign_rs::nistapi::{crypto_sign, crypto_sign_keypair};

    let seeds: &[&str] = &[
        "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1",
    ];
    let msg = b"audit";
    for s in seeds {
        let mut seed = [0u8; 48];
        for i in 0..48 {
            seed[i] = u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap();
        }
        ctrdrbg::randombytes_init(&seed);
        let (_, sk) = crypto_sign_keypair().unwrap();
        let _ = crypto_sign(msg, &sk).unwrap();
    }
    #[cfg(ibz_audit)]
    sqisign_rs::quaternion::ibz_audit_report();
}
