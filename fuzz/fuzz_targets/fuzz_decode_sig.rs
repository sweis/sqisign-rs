#![no_main]
//! Fuzz signature byte decode/encode roundtrip. Must never panic.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::verification::Signature;

fuzz_target!(|data: &[u8]| {
    if let Ok(sig) = Signature::try_from(data) {
        let _ = sig.to_bytes();
    }
});
