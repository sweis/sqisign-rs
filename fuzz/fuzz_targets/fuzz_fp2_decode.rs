#![no_main]
//! Fuzz `fp2_decode`: if input is accepted as canonical, encode must roundtrip.

use libfuzzer_sys::fuzz_target;
use sqisign_rs::gf::{fp2_decode, fp2_encode, Fp2, FP2_ENCODED_BYTES};

fuzz_target!(|data: &[u8]| {
    if data.len() < FP2_ENCODED_BYTES {
        return;
    }
    let src = &data[..FP2_ENCODED_BYTES];
    let mut x = Fp2::default();
    let ok = fp2_decode(&mut x, src);
    if ok == 0xFFFF_FFFF {
        let mut out = [0u8; FP2_ENCODED_BYTES];
        fp2_encode(&mut out, &x);
        assert_eq!(&out[..], src, "canonical input must roundtrip");
    }
});
