# Controlled C-vs-Rust performance harnesses

These run identical workloads (KAT-0 seed, deterministic) so that
`perf stat -e instructions:u` gives directly comparable counts.

## C side

Requires the C reference built at `../the-sqisign/build-sysgmp` with
`-DSQISIGN_BUILD_TYPE=ref -DGMP_LIBRARY=SYSTEM -DCMAKE_BUILD_TYPE=Release`.

```sh
B=/root/src/personal-hacking/the-sqisign/build-sysgmp
C=/root/src/personal-hacking/the-sqisign

# verify ×200
gcc -O3 -DNDEBUG -DSQISIGN_VARIANT=lvl1 -DRADIX_64 -DSQISIGN_BUILD_TYPE_REF \
  -DSQISIGN_GF_IMPL_REF -I$C/include -I$C/src/nistapi/lvl1 \
  c_verify_bench.c -o /tmp/c_verify_bench \
  -Wl,--start-group $B/src/libsqisign_lvl1_nistapi.a $B/src/libsqisign_lvl1.a \
  $B/src/verification/ref/lvl1/libsqisign_verification_lvl1.a \
  $B/src/hd/ref/lvl1/libsqisign_hd_lvl1.a $B/src/ec/ref/lvl1/libsqisign_ec_lvl1.a \
  $B/src/gf/ref/lvl1/libsqisign_gf_lvl1.a $B/src/mp/ref/generic/libsqisign_mp_generic.a \
  $B/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a \
  $B/src/common/generic/libsqisign_common_sys.a \
  $B/src/signature/ref/lvl1/libsqisign_signature_lvl1.a \
  $B/src/id2iso/ref/lvl1/libsqisign_id2iso_lvl1.a \
  $B/src/quaternion/ref/generic/libsqisign_quaternion_generic.a \
  -Wl,--end-group -lgmp -lm
perf stat -e instructions:u,cycles:u /tmp/c_verify_bench

# keygen / sign ×20
gcc -O3 -DNDEBUG -DSQISIGN_VARIANT=lvl1 -DRADIX_64 -DSQISIGN_BUILD_TYPE_REF \
  -DSQISIGN_GF_IMPL_REF -DENABLE_SIGN -I$C/include -I$C/src/nistapi/lvl1 \
  -I$C/src/common/ref/include c_kgsign_bench.c -o /tmp/c_kgsign_bench \
  -Wl,--start-group $B/src/libsqisign_lvl1_test_nistapi.a $B/src/libsqisign_lvl1_test.a \
  $B/src/signature/ref/lvl1/libsqisign_signature_lvl1.a \
  $B/src/id2iso/ref/lvl1/libsqisign_id2iso_lvl1.a \
  $B/src/quaternion/ref/generic/libsqisign_quaternion_generic.a \
  $B/src/verification/ref/lvl1/libsqisign_verification_lvl1.a \
  $B/src/hd/ref/lvl1/libsqisign_hd_lvl1.a $B/src/ec/ref/lvl1/libsqisign_ec_lvl1.a \
  $B/src/gf/ref/lvl1/libsqisign_gf_lvl1.a $B/src/mp/ref/generic/libsqisign_mp_generic.a \
  $B/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a \
  $B/src/common/generic/libsqisign_common_test.a \
  -Wl,--end-group -lgmp -lm
perf stat -e instructions:u,cycles:u /tmp/c_kgsign_bench kg
perf stat -e instructions:u,cycles:u /tmp/c_kgsign_bench kgsign
```

## Rust side

These are not registered as `[[example]]` (they're moved out of `examples/` to
keep `cargo build` clean). Run them ad hoc:

```sh
mkdir -p examples && cp tools/perf/rust_*.rs examples/
# add to Cargo.toml [[example]] entries with required-features = ["sign"] for kg_sign
cargo build --release --features lvl1,sign --example rust_kg_sign
perf stat -e instructions:u,cycles:u ./target/release/examples/rust_kg_sign kg
perf stat -e instructions:u,cycles:u ./target/release/examples/rust_kg_sign kgsign
```

Per-op derivation: `keygen = kg/20`, `sign = (kgsign − kg)/20`, `verify = verify200/200`.

## Caveat

The C `ref` build uses bitsliced software AES (`br_aes_ct64_*`); the Rust `aes`
crate auto-detects AES-NI. This dominates keygen/sign DRBG cost — see PORTING.md.
