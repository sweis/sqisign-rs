// SPDX-License-Identifier: Apache-2.0
//! Known Answer Tests against NIST KAT response files.

use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone)]
struct KatVector {
    count: usize,
    seed: Vec<u8>,
    msg: Vec<u8>,
    pk: Vec<u8>,
    sk: Vec<u8>,
    sm: Vec<u8>,
}

fn parse_kat_file(path: &str) -> Vec<KatVector> {
    let content = fs::read_to_string(path).expect("KAT file");
    let mut vectors = Vec::new();
    let mut cur: Option<KatVector> = None;

    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("count = ") {
            if let Some(v) = cur.take() {
                vectors.push(v);
            }
            cur = Some(KatVector {
                count: rest.parse().unwrap(),
                seed: vec![],
                msg: vec![],
                pk: vec![],
                sk: vec![],
                sm: vec![],
            });
        } else if let Some(v) = cur.as_mut() {
            if let Some(rest) = line.strip_prefix("seed = ") {
                v.seed = hex::decode(rest).unwrap();
            } else if let Some(rest) = line.strip_prefix("msg = ") {
                v.msg = hex::decode(rest).unwrap();
            } else if let Some(rest) = line.strip_prefix("pk = ") {
                v.pk = hex::decode(rest).unwrap();
            } else if let Some(rest) = line.strip_prefix("sk = ") {
                v.sk = hex::decode(rest).unwrap();
            } else if let Some(rest) = line.strip_prefix("sm = ") {
                v.sm = hex::decode(rest).unwrap();
            }
        }
    }
    if let Some(v) = cur.take() {
        vectors.push(v);
    }
    vectors
}

fn kat_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("KAT");
    #[cfg(all(feature = "lvl1", not(feature = "lvl3"), not(feature = "lvl5")))]
    p.push("PQCsignKAT_353_SQIsign_lvl1.rsp");
    #[cfg(feature = "lvl3")]
    p.push("PQCsignKAT_529_SQIsign_lvl3.rsp");
    #[cfg(feature = "lvl5")]
    p.push("PQCsignKAT_701_SQIsign_lvl5.rsp");
    p
}

#[test]
fn kat_file_parses() {
    let vectors = parse_kat_file(kat_path().to_str().unwrap());
    assert_eq!(vectors.len(), 100);
    assert_eq!(vectors[0].seed.len(), 48);
    assert_eq!(
        vectors[0].pk.len(),
        sqisign_rs::params::CRYPTO_PUBLICKEYBYTES
    );
    assert_eq!(
        vectors[0].sk.len(),
        sqisign_rs::params::CRYPTO_SECRETKEYBYTES
    );
}

#[test]
fn kat_verify() {
    use sqisign_rs::nistapi;
    let vectors = parse_kat_file(kat_path().to_str().unwrap());
    for v in &vectors {
        let m = nistapi::crypto_sign_open(&v.sm, &v.pk).expect("verify");
        assert_eq!(m, v.msg, "count {}", v.count);
    }
}

#[test]
fn kat_verify_reject_tampered() {
    use sqisign_rs::nistapi;
    let vectors = parse_kat_file(kat_path().to_str().unwrap());
    for v in vectors.iter().take(5) {
        // Flip one signature byte; verification must fail.
        let mut sm = v.sm.clone();
        sm[10] ^= 0x01;
        assert!(
            nistapi::crypto_sign_open(&sm, &v.pk).is_err(),
            "count {}",
            v.count
        );
        // Flip one message byte; verification must fail.
        let mut sm = v.sm.clone();
        let last = sm.len() - 1;
        sm[last] ^= 0x01;
        assert!(
            nistapi::crypto_sign_open(&sm, &v.pk).is_err(),
            "count {}",
            v.count
        );
    }
}

#[test]
#[cfg(feature = "sign")]
fn kat_sign() {
    use sqisign_rs::{common::ctrdrbg, nistapi};
    let vectors = parse_kat_file(kat_path().to_str().unwrap());
    let limit = std::env::var("KAT_SIGN_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(vectors.len());
    for v in vectors.iter().take(limit) {
        ctrdrbg::randombytes_init(&v.seed);
        let (pk, sk) = nistapi::crypto_sign_keypair().expect("keypair");
        assert_eq!(pk, v.pk.as_slice(), "pk count {}", v.count);
        assert_eq!(sk, v.sk.as_slice(), "sk count {}", v.count);
        let sm = nistapi::crypto_sign(&v.msg, &sk).expect("sign");
        assert_eq!(sm, v.sm, "sm count {}", v.count);
    }
    #[cfg(ibz_audit)]
    sqisign_rs::quaternion::ibz_audit_report();
}
