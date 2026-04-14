#![no_main]
//! Fuzz public key byte decode. Must never panic.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::precomp::PUBLICKEY_BYTES;
use sqisign_rs::verification::{public_key_from_bytes, PublicKey};

fuzz_target!(|data: &[u8]| {
    if data.len() < PUBLICKEY_BYTES {
        return;
    }
    let mut pk = PublicKey::default();
    public_key_from_bytes(&mut pk, &data[..PUBLICKEY_BYTES]);
});
