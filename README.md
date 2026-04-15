# sqisign-rs

[![CI](https://github.com/sweis/sqisign-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sweis/sqisign-rs/actions/workflows/ci.yml)

A Rust implementation of [SQIsign](https://sqisign.org/), a post-quantum
digital signature scheme based on isogenies of supersingular elliptic curves.

Ported from the [C reference implementation](https://github.com/SQIsign/the-sqisign)
and validated against the official NIST KAT test vectors.

> **Warning:** SQIsign is a NIST PQC Round 2 candidate and `sqisign-rs` is AI generated. It
> has not been independently audited. Do not use in production.

## Status

Performance is at parity with the C `ref` build on identical workloads. See [PORTING.md](PORTING.md) 
for C-vs-spec discrepancies found during the port, security hardening notes, and mutation-testing results.

## Building

Requires Rust 1.75+ (stable). Enable the pre-commit hook (rustfmt + clippy):

```sh
git config core.hooksPath .githooks
```

The signing path uses arbitrary-precision
integers; the default backend is [`malachite`](https://crates.io/crates/malachite-nz)
(pure Rust, no C dependency).

At lvl1 the default GF(p) backend is a vendored saturated 4-limb Montgomery
implementation with hand-written x86_64 BMI2/ADX inline asm for `mul`/`sqr`;
build with `RUSTFLAGS="-C target-cpu=native"` to enable the asm path.

```sh
# Full build (keygen + sign + verify)
RUSTFLAGS="-C target-cpu=native" cargo build --release --features lvl1,sign

# Opt-in GMP backend via rug (requires libgmp)
cargo build --release --features lvl1,gmp

# Verification only (no bigint dependency at all)
cargo build --release --no-default-features --features lvl1
```

## Testing

```sh
# Unit tests + KAT + Wycheproof-style edge cases (~11s release at lvl1)
cargo test --release --features lvl1,sign

# Run only the first N sign vectors
KAT_SIGN_LIMIT=5 cargo test --features lvl1,sign kat_sign

# Mutation testing
cargo mutants --no-default-features --features lvl1
```

The Wycheproof-style suite (`tests/wycheproof.rs`) generates ~1.8K-3.4K edge
cases per level (bit-flip sweep, non-canonical encodings, structure boundaries,
malleability). Export JSON vectors with `WYCHEPROOF_EXPORT=1`.

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
- `sign` — enable key generation and signing using the pure-Rust `malachite`
  bigint backend. Without this flag the crate builds verification only,
  matching the C reference's `ENABLE_SIGN=OFF` mode.
- `gmp` — use GMP (via `rug`) instead of `malachite` for the bigint backend.
  Slightly faster (~5% on sign); requires libgmp.
- `cryptobigint` — use [`crypto-bigint`](https://crates.io/crates/crypto-bigint)
  as a fixed-precision bigint backend. Much slower (~400× on sign) since every
  integer is sized to the worst-case HNF intermediate; provided for users who
  want a minimal, audited dependency set.
- `gf-portable` — opt out of the vendored 4-limb Montgomery GF backend at lvl1
  and use the in-tree unsaturated modarith reference instead. lvl3/lvl5 always
  use modarith.
- `bench-internals` — exposes DRBG snapshot/restore for deterministic benchmark
  replays. Do not enable in production.

## Benchmarking

```sh
tools/record_bench.sh "my-change"
```

appends a `[commit, timestamp, verify, keygen, sign, note]` row to
[`BENCH_HISTORY.csv`](BENCH_HISTORY.csv). For a controlled C-vs-Rust comparison
see `tools/perf/README.md`.

## Fuzzing

Fuzz targets live in `fuzz/` and require `cargo install cargo-fuzz` (nightly
toolchain). Targets cover the verification path only (no GMP dependency).

```sh
cargo fuzz run fuzz_verify
```

Available targets: `fuzz_verify`, `fuzz_decode_sig`, `fuzz_decode_pk`,
`fuzz_fp2_decode`. A small KAT-derived seed corpus is checked in under
`fuzz/corpus/fuzz_verify/`; regenerate with `python3 fuzz/seed_corpus.py`.
