#![no_main]
//! Fuzz public key byte decode. Must never panic.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::verification::PublicKey;

fuzz_target!(|data: &[u8]| {
    let _ = PublicKey::try_from(data);
});
