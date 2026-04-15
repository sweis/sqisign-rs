#![no_main]
//! Fuzz `Fp2::try_decode`: if input is accepted as canonical, encode must roundtrip.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::gf::{Fp2, FP2_ENCODED_BYTES};

fuzz_target!(|data: &[u8]| {
    if data.len() < FP2_ENCODED_BYTES {
        return;
    }
    let src = &data[..FP2_ENCODED_BYTES];
    if let Some(x) = Fp2::try_decode(src) {
        assert_eq!(&x.encode()[..], src, "canonical input must roundtrip");
    }
});
