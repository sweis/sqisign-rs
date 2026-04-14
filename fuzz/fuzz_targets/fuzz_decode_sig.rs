#![no_main]
//! Fuzz signature byte decode/encode roundtrip. Must never panic.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::precomp::SIGNATURE_BYTES;
use sqisign_rs::verification::{signature_from_bytes, signature_to_bytes, Signature};

fuzz_target!(|data: &[u8]| {
    if data.len() < SIGNATURE_BYTES {
        return;
    }
    let enc = &data[..SIGNATURE_BYTES];
    let mut sig = Signature::default();
    signature_from_bytes(&mut sig, enc);
    let mut out = [0u8; SIGNATURE_BYTES];
    signature_to_bytes(&mut out, &sig);
});
