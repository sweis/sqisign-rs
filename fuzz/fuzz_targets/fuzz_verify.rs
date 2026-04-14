#![no_main]
//! Fuzz `crypto_sign_open` with arbitrary (pk || sm). Must never panic.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::nistapi::crypto_sign_open;
use sqisign_rs::params::{CRYPTO_BYTES, CRYPTO_PUBLICKEYBYTES};

fuzz_target!(|data: &[u8]| {
    if data.len() < CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES {
        return;
    }
    let pk = &data[..CRYPTO_PUBLICKEYBYTES];
    let sm = &data[CRYPTO_PUBLICKEYBYTES..];
    let _ = crypto_sign_open(sm, pk);
});
