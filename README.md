# sqisign-rs

[![CI](https://github.com/sweis/sqisign-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sweis/sqisign-rs/actions/workflows/ci.yml)

A Rust implementation of [SQIsign](https://sqisign.org/), a post-quantum
digital signature scheme based on isogenies of supersingular elliptic curves.

Ported from the [C reference implementation](https://github.com/SQIsign/the-sqisign)
and validated against the official NIST KAT test vectors.

> **Warning:** SQIsign is a NIST PQC Round 2 candidate and `sqisign-rs` is AI generated. It
> has not been independently audited. Do not use in production.

## Status

See [PORTING.md](PORTING.md) for C-vs-spec discrepancies found during the
port, security hardening notes, and mutation-testing results.

## Backends

There are two pluggable backend axes; both default to the fastest option
available on the build target.

**GF(p)/GF(p²) field arithmetic** (verification hot path; lvl1 only):

- *default* — vendored saturated 4-limb Montgomery with x86_64 inline asm
  (`MULX`/`ADCX`/`ADOX` for `fp_mul`/`fp_sqr` and a fused sum/diff-of-products
  kernel for `Fp2::mul`). On CPUs with AVX-512 IFMA, the theta-isogeny inner
  loop additionally runs in batched 8-way radix-2⁵² SoA. Requires
  `RUSTFLAGS="-C target-cpu=native"` (or at minimum
  `-C target-feature=+adx,+bmi2`); without it the build script falls back to
  a portable intrinsics path with a warning.
- `gf-portable` — opt-out to the in-tree unsaturated 5-limb modarith
  reference (the readable Rust-only implementation). lvl3/lvl5 always use
  this backend.

**Arbitrary-precision integers** (signing only; verification has no bigint
dependency):

- *default* (`sign` feature) — [`malachite`](https://crates.io/crates/malachite-nz),
  pure Rust.
- `gmp` — GMP via [`rug`](https://crates.io/crates/rug). ~5% faster signing;
  requires `libgmp`.
- `cryptobigint` — [`crypto-bigint`](https://crates.io/crates/crypto-bigint)
  fixed-precision. ~400× slower (every integer is sized to the worst-case
  HNF intermediate); provided for users who want a minimal dependency set.

## Performance

lvl1 on Intel Xeon Platinum 8488C (Sapphire Rapids), KAT-0, release, min-of-N:

| GF backend | bigint | verify | keygen | sign |
|---|---|---|---|---|
| asm + AVX-512 IFMA | gmp | **1.46 ms** | **13.8 ms** | **33.8 ms** |
| asm + AVX-512 IFMA | malachite | 1.44 ms | 15.9 ms | 39.3 ms |
| asm (BMI2/ADX) | gmp | 1.64 ms | 16.0 ms | 36.5 ms |
| asm (BMI2/ADX) | malachite | 1.54 ms | 16.5 ms | 39.7 ms |
| `gf-portable` | gmp | 3.20 ms | 23.0 ms | 50.6 ms |
| `gf-portable` | malachite | 3.16 ms | 22.8 ms | 53.6 ms |
| C reference (`broadwell`) | (gmp) | 1.56 ms | — | — |
| C reference (`ref`) | (gmp) | 2.50 ms | 16.9 ms | 40.0 ms |

`asm + AVX-512 IFMA` with malachite is the default; build with
`RUSTFLAGS="-C target-cpu=native"`. The GF backend determines verify time;
the bigint backend determines the keygen/sign delta. See
[`BENCH_HISTORY.csv`](BENCH_HISTORY.csv) and `tools/perf/README.md` for the
controlled instruction-count comparison.

## Building

Requires Rust 1.75+ (stable). Enable the pre-commit hook (rustfmt + clippy):

```sh
git config core.hooksPath .githooks
```

```sh
# Full build (keygen + sign + verify)
RUSTFLAGS="-C target-cpu=native" cargo build --release --features lvl1,sign

# Opt-in GMP bigint backend (requires libgmp)
RUSTFLAGS="-C target-cpu=native" cargo build --release --features lvl1,gmp

# Pure-Rust reference (no asm, no GMP)
cargo build --release --no-default-features --features lvl1,gf-portable,sign

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

- `lvl1` / `lvl3` / `lvl5` — security level (mutually exclusive).
- `sign` — enable keygen and signing (default bigint: malachite). Without
  this flag the crate is verification-only.
- `gmp` / `cryptobigint` / `gf-portable` — see [Backends](#backends).
- `bench-internals` — exposes DRBG snapshot/restore for deterministic
  benchmark replays. Do not enable in production.

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
