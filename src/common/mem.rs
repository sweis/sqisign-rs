// SPDX-License-Identifier: Apache-2.0
//! Secure memory clearing. Thin wrapper over `zeroize`.

use zeroize::Zeroize;

/// Zero `mem` in a way the optimizer will not elide.
/// Mirrors `sqisign_secure_clear` in the C code.
pub fn secure_clear(mem: &mut [u8]) {
    mem.zeroize();
}
