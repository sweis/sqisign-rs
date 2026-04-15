// SPDX-License-Identifier: Apache-2.0
//! Test-only helpers shared across modules.

use crate::gf::{fp_decode_reduce, Fp, Fp2, FP_ENCODED_BYTES};

/// xorshift64* PRNG for property tests (not cryptographic).
pub struct Prng(pub u64);

impl Prng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    pub fn fill(&mut self, buf: &mut [u8]) {
        for chunk in buf.chunks_mut(8) {
            let x = self.next().to_le_bytes();
            chunk.copy_from_slice(&x[..chunk.len()]);
        }
    }
    pub fn fp(&mut self) -> Fp {
        let mut buf = [0u8; FP_ENCODED_BYTES];
        self.fill(&mut buf);
        let mut a = Fp::default();
        fp_decode_reduce(&mut a, &buf);
        a
    }
    pub fn fp2(&mut self) -> Fp2 {
        Fp2 {
            re: self.fp(),
            im: self.fp(),
        }
    }
}

#[track_caller]
#[allow(dead_code)] // golden-vector tests are lvl1-only
pub fn assert_hex(buf: &[u8], expected: &str) {
    assert_eq!(hex::encode(buf), expected);
}
