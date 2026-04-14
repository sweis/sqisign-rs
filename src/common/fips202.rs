// SPDX-License-Identifier: Apache-2.0
//! SHAKE256 wrappers matching the C `fips202.h` incremental API.
//!
//! The C code only uses SHAKE256 (one-shot and `shake256_inc_*`), so only
//! those are provided. Backed by the `sha3` crate.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::{Shake256, Shake256Reader};

pub const SHAKE256_RATE: usize = 136;

/// One-shot SHAKE256: hash `input` and squeeze `output.len()` bytes.
pub fn shake256(output: &mut [u8], input: &[u8]) {
    let mut hasher = Shake256::default();
    hasher.update(input);
    hasher.finalize_xof().read(output);
}

/// Incremental SHAKE256 context mirroring the C `shake256incctx` lifecycle:
/// `init` → `absorb`* → `finalize` → `squeeze`*.
///
/// Clone is supported (matches `shake256_inc_ctx_clone`). Dropping the value
/// releases it (matches `shake256_inc_ctx_release`).
#[derive(Clone)]
pub struct Shake256Inc {
    state: State,
}

#[derive(Clone)]
enum State {
    Absorbing(Shake256),
    Squeezing(Shake256Reader),
}

impl Default for Shake256Inc {
    fn default() -> Self {
        Self::new()
    }
}

impl Shake256Inc {
    /// `shake256_inc_init`
    pub fn new() -> Self {
        Self {
            state: State::Absorbing(Shake256::default()),
        }
    }

    /// `shake256_inc_absorb` — may be called repeatedly before `finalize`.
    pub fn absorb(&mut self, input: &[u8]) {
        match &mut self.state {
            State::Absorbing(h) => h.update(input),
            State::Squeezing(_) => panic!("absorb called after finalize"),
        }
    }

    /// `shake256_inc_finalize` — switch from absorb to squeeze phase.
    pub fn finalize(&mut self) {
        let old = core::mem::replace(&mut self.state, State::Absorbing(Shake256::default()));
        match old {
            State::Absorbing(h) => self.state = State::Squeezing(h.finalize_xof()),
            State::Squeezing(r) => self.state = State::Squeezing(r),
        }
    }

    /// `shake256_inc_squeeze` — may be called repeatedly after `finalize`.
    pub fn squeeze(&mut self, output: &mut [u8]) {
        match &mut self.state {
            State::Squeezing(r) => r.read(output),
            State::Absorbing(_) => panic!("squeeze called before finalize"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shake256_nist_vector() {
        // NIST CAVP: SHAKE256("") first 32 bytes
        let mut out = [0u8; 32];
        shake256(&mut out, b"");
        assert_eq!(
            hex::encode(out),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
    }

    #[test]
    fn shake256_incremental_matches_oneshot() {
        let mut one = [0u8; 64];
        shake256(&mut one, b"hello world");

        let mut ctx = Shake256Inc::new();
        ctx.absorb(b"hello ");
        ctx.absorb(b"world");
        ctx.finalize();
        let mut a = [0u8; 20];
        let mut b = [0u8; 44];
        ctx.squeeze(&mut a);
        ctx.squeeze(&mut b);

        assert_eq!(&one[..20], &a[..]);
        assert_eq!(&one[20..], &b[..]);
    }

    #[test]
    fn shake256_inc_clone() {
        let mut ctx = Shake256Inc::new();
        ctx.absorb(b"abc");
        ctx.finalize();
        let mut ctx2 = ctx.clone();
        let mut o1 = [0u8; 16];
        let mut o2 = [0u8; 16];
        ctx.squeeze(&mut o1);
        ctx2.squeeze(&mut o2);
        assert_eq!(o1, o2);
    }
}
