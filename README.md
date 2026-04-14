# sqisign-rs

A Rust implementation of [SQIsign](https://sqisign.org/), a post-quantum
digital signature scheme based on isogenies of supersingular elliptic curves.

Ported from the [C reference implementation](https://github.com/SQIsign/the-sqisign)
and validated against the official NIST KAT test vectors.

> **Warning:** SQIsign is a NIST PQC Round 2 candidate and this port has not
> been independently audited. Do not use in production.

## Status

| Level | Prime | Verify | Sign / Keygen | KAT |
|-------|-------|--------|---------------|-----|
| lvl1  | 5·2²⁴⁸ − 1 | ✅ | ✅ | 100/100 sign + 100/100 verify |
| lvl3  | 65·2³⁷⁶ − 1 | — | — | pending |
| lvl5  | 27·2⁵⁰⁰ − 1 | — | — | pending |

See [PORTING.md](PORTING.md) for module status, C-vs-spec discrepancies found
during the port, and design decisions.

## Building

Requires Rust 1.75+ (stable). The signing path additionally requires GMP
(via the `rug` crate) for arbitrary-precision integer arithmetic.

```sh
# Full build (keygen + sign + verify)
cargo build --release --features lvl1,sign

# Verification only (no GMP dependency)
cargo build --release --no-default-features --features lvl1
```

## Testing

```sh
# All unit tests + KAT integration tests (~8s release)
cargo test --release --features lvl1,sign

# Run only the first N sign vectors (each takes ~80ms)
KAT_SIGN_LIMIT=5 cargo test --features lvl1,sign kat_sign
```

## Usage

```rust
use sqisign_rs::common::ctrdrbg;
use sqisign_rs::nistapi::{crypto_sign, crypto_sign_keypair, crypto_sign_open};
use sqisign_rs::params::CRYPTO_BYTES;

// Seed the deterministic RNG (in production, draw 48 bytes from the OS RNG).
ctrdrbg::randombytes_init(&[42u8; 48]);

// Generate a keypair.
let (pk, sk) = crypto_sign_keypair().unwrap();

// Sign. Output is signature || message.
let msg = b"Per aspera ad astra";
let sm = crypto_sign(msg, &sk).unwrap();

// Verify a valid signature: returns the recovered message.
let recovered = crypto_sign_open(&sm, &pk).unwrap();
assert_eq!(recovered, msg);

// A tampered signature is rejected.
let mut bad = sm.clone();
bad[0] ^= 0x01;
assert!(crypto_sign_open(&bad, &pk).is_err());
```

A runnable version is in [`examples/sign_verify.rs`](examples/sign_verify.rs):

```sh
cargo run --release --features lvl1,sign --example sign_verify
```

## Feature flags

- `lvl1` / `lvl3` / `lvl5` — select the security level (mutually exclusive;
  exactly one must be enabled).
- `sign` — enable key generation and signing. Pulls in `rug` (GMP). Without
  this flag the crate builds verification only, matching the C reference's
  `ENABLE_SIGN=OFF` mode.

## Fuzzing

Fuzz targets live in `fuzz/` and require `cargo install cargo-fuzz` (nightly
toolchain). Targets cover the verification path only (no GMP dependency).

```sh
cargo fuzz run fuzz_verify
```

Available targets: `fuzz_verify`, `fuzz_decode_sig`, `fuzz_decode_pk`,
`fuzz_fp2_decode`. A small KAT-derived seed corpus is checked in under
`fuzz/corpus/fuzz_verify/`; regenerate with `python3 fuzz/seed_corpus.py`.

## License

Apache-2.0. See `LICENSE`.
