# SQIsign C→Rust Porting Status

Tracking the port from <https://github.com/SQIsign/the-sqisign> (reference impl).

## Module status (lvl1)

| Module | C LOC | Status | Notes |
|--------|-------|--------|-------|
| common | 4557 | **done** | SHAKE256 via `sha3`, AES-CTR-DRBG via `aes`; 6 tests |
| mp | 445 | **done** | fixed-precision limb arithmetic; 13 tests |
| gf | 15182 | **done** | fp_p5248_64 (5-limb unsaturated) + fp2; 27 tests, C-cross-checked |
| precomp (verify) | ~12K | **done** | constants for verification path |
| ec | 4766 | **done** | x-only/Jacobian/isog/basis (biextension stubbed); 10 tests |
| hd | 2042 | **done** | theta isogenies; 4 tests |
| verification | 740 | **done** | **all 100 KAT verify vectors pass** (0.55s) + tamper-reject |
| nistapi | 193 | **done** | crypto_sign_open works |
| quaternion (core) | ~2400 | **done** | intbig/dim/algebra/integers via `rug`; 30 tests |
| quaternion (adv) | ~2200 | **done** | HNF/LLL/lattice/ideal/normeq; clean-room DPE; 26 tests |
| precomp (sign) | ~263K | **done** | scripted port via `tools/gen_precomp_sign.py`; 8 tests |
| ec (biextension) | 770 | **done** | Weil/Tate pairings + dlog; 3 tests |
| id2iso | 2594 | **done** | id2iso + dim2id2iso (find_uv, clapotis); 2 tests |
| signature | 1370 | **done** | **all 100 KAT sign vectors pass** (8s release) |

**Status: full lvl1 keygen/sign/verify parity with C reference.** 131 unit tests + 4 KAT integration tests (200 vectors total exercised).

## Test vectors

- KAT files copied to `KAT/` from C repo (100 vectors per level).
- `tests/kat.rs` parses .rsp files; verify/sign tests are `#[ignore]` until impl lands.

## C oddities flagged for spec cross-check

Each item analysed against the SQIsign Round-2 spec (2025-07-07) and call-site usage.

| Item | Verdict | Notes |
|---|---|---|
| `mp_mul2` omits `a[1]*b[0]` | **Dead code** | Never called anywhere. Keep faithful port. |
| `mp_neg` no carry propagation | **Correct in context** | Only used by `mp_inv_2e`/`mp_invert_matrix`, both mod-2^e where carry into the modulus bit is irrelevant. Keep faithful. |
| `mp_shiftr` returns `x[0]&1` | **Correct in context** | Only callers (ec.c:438-439) pass `shift=1`. Return value wrong for `shift>1` but unreached. Keep faithful. |
| `ec_is_equal` precedence | **Doc bug, harmless** | `~a & ~b * c` evaluates as `~a & (~b * c)`; with all-ones masks `(-1)*(-1)=1` (u32), so returns 1 not 0xFFFFFFFF. Doc-comment claims 0xFFFFFFFF. All callers test as boolean. Keep faithful; consider fixing doc upstream. |
| `xMUL` calls `xDBLADD(...,true)` unconditionally | **Latent bug, unreached** | If `is_A24_computed_and_normalized==false`, local A24 has z=4C and the `true` flag wrongly skips the z-multiply at xDBLADD:325. All current callers (`basis.c` via `ec_curve_normalize_A24`; `sign.c` after curve normalization) ensure normalization first. Keep faithful; worth an upstream `assert`. |
| `jac_from_ws` skips x-copy when A=0 | **Latent bug, unreached** | Would leave Q->x uninitialized for distinct Q,P. All callers pass `Q==P` (same pointer). Keep faithful. |
| `ibz_rand_interval_bits` subtracts `m` | **Real bug, test-only** | intbig.c:569 `mpz_sub_ui(*rand,*rand,m)` shifts the interval to `[-2^m-m, 2^m-m]`. Only callers are in `test/`. Keep faithful (irrelevant to KAT). |
| `ibz_div` doc | **No issue** | Header (intbig.h:105) actually says "rounded towards zero" — matches `mpz_tdiv_qr`. Earlier port note was wrong. |
| `ibz_mat_4x4_is_hnf` `linestart` never updated | **Real bug, debug-only** | hnf.c:37-46: `linestart` stays -1 so `linestart < i` is vacuous. Only used in NDEBUG asserts. Keep faithful. |
| `quat_sampling_random_ideal_O0_given_norm` clobbers `found` | **Real bug vs spec** | normeq.c:319 resets `found=0` to enter the rerandomize loop, discarding the represent_integer success flag. Spec Alg 3.10 step 10 says "might raise an exception"; C silently proceeds with an invalid `gen`. Unreachable with valid SQIsign parameters (represent_integer succeeds w.h.p.). Keep faithful for KAT; worth upstream issue. |
| `normeq.c:188` C `%` on signed | **Real bug vs spec; Rust must match C** | Spec Alg 3.12 step 15: `(x−t ≡ 2 mod 4)`. C uses `(ibz_get(x)-ibz_get(t)) % 4 == 2`; for negative differences C gives `-2`, failing the check. Rejects valid samples → extra retry iterations → different DRBG consumption. **The Rust port currently uses `rem_euclid(4)` which is spec-correct but breaks KAT determinism. Must revert to C's truncating `%` for KAT parity.** |

## Spec discrepancies found (verify.c vs Algorithm 4.9)

- **Extra rejection of e'_rsp = 1**: verify.c:258 rejects `pow_dim2_deg_resp == 1`; spec only requires `≥ 0` (step 7). C comment: "dim2 isogeny embeds odd-degree dim1 isogeny, never length 2." Honest signers never produce e'_rsp=1, so this is defense-in-depth. Spec should arguably include it.
- **Matrix storage is transposed vs spec notation**: C stores `mat[i][j]` = spec `M[j][i]`. Self-consistent (encode/decode/apply all use the same convention); not a bug. Spec step 16 indices `M[0][0],M[0][1]` correspond to C's `mat[0][0],mat[1][0]` — match.
- **Basis-scaling exponent**: spec text reads "2^{f−e'_rsp+2}"; C computes `f − e'_rsp − 2`. The PDF rendering drops parentheses; spec intent is `f−(e'_rsp+2)`. Match.
- **Supersingularity check (spec steps 3-4)** is implicit in C (per spec §4.5 prose: "byproduct of a successful computation of an isogeny"). C checks `ec_curve_verify_A` (A∉{0,±2}) for both curves, then relies on isogeny computation to fail for non-supersingular input. Matches spec's stated implementation strategy.

## Bugs found in our own port during development

- `gf/fp.rs`: `NRES_C[4]` had a hex typo (12 c's vs 11) breaking nres/decode.
- `quaternion/dpe.rs`: clean-room `set_z` extracted 64 bits with rounding instead of 53-bit truncation (mismatching `mpz_get_d_2exp`).
- `quaternion/dpe.rs`: clean-room `round` used ties-to-even; C `round()` is half-away-from-zero. This single tie-case at LLL step 1102 (`u=4.5`) desynchronized the entire KAT-0 keygen.
- `quaternion/normeq.rs`: used spec-correct `rem_euclid(4)` instead of C's truncating `%`, breaking DRBG sync.

## Security hardening (Rust improvements over C)

| Area | Finding | Status |
|---|---|---|
| **DoS via `backtracking`** | A malicious signature with `backtracking ≥ 128` makes `check_canonical_basis_change_matrix` compute a negative shift; `as u32` wraps to ~4e9 and `multiple_mp_shiftl` burns ~18ms (release) or panics on `>> 64` overflow at exactly 128 (debug). **Same bug in C (UB shift).** | **Fixed**: early `return false` when `shift ≤ 0`. Output unchanged (signature would be rejected later anyway). Regression test `reject_adversarial_backtracking`. Worth upstream issue. |
| Constant-time | `select_ct`/`swap_ct`/`fp_select`/`fp_cswap` are mask-xor based, no branches. `xDBLMUL` ladder is CT. The `kbits=1` branch in `ec_biscalar_mul` is data-dependent but only reachable with `f=1`, which never occurs on the verify path (always `f ≥ HD_EXTRA_TORSION=2`). Signing is variable-time by design (uses GMP). | OK |
| Panics on adversarial input (verify) | Audited `unwrap`/`expect`/`panic!`/`unreachable!` reachable from `crypto_sign_open`: `ec/mod.rs:791` and `hd/theta_isogenies.rs:51` are genuinely unreachable (`x & 1` ∈ {0,1}, `x & 3` ∈ {0..3}). `try_into().unwrap()` calls are on fixed-size stack-array slices (infallible). `fips202` panics are programmer-error only. **Exception**: `find_na_x_coord` `debug_assert!` at basis.rs:245 fires if attacker supplies `hint` with `hint_a` bit not matching `is_square(A)` — debug-mode crash only; release proceeds and rejects later. | Debug-only; matches C `assert`. |
| Loop bounds | `find_na_x_coord`/`find_nqr_factor` loops are probabilistically bounded (~2 iterations expected, density ½). All other verify-path loops bounded by `TORSION_EVEN_POWER`/`SQISIGN_RESPONSE_LENGTH` constants. | OK |
| `unsafe` | Single block: `intbig.rs` calls `gmp::mpz_get_si` directly for exact `ibz_get` semantics (signing-only, KAT-determinism). `#![warn(unsafe_code)]` lint enabled. | Justified |
| Secret zeroization | `DrbgState` derives `ZeroizeOnDrop`. `SecretKey` impls `Drop` zeroizing field-element/basis components and overwriting `rug::Integer` values with 0. **Limitation**: GMP may have realloc'd intermediate limb buffers that are not scrubbed; `rug` provides no zeroize hook. The serialized sk byte buffer should also be zeroized by callers. | Best-effort |
| Unused deps | `subtle` was declared but never used (CT primitives are hand-masked). | **Removed** |

## Performance

Controlled comparison on **identical workload** (KAT-0 seed, deterministic), measured
with `perf stat -e instructions:u` (run-to-run stable to ±1k; cycle counts are
noisy under concurrent load and shown for reference only). C is the `ref` build
linked against system GMP 6.2.1. Harnesses: `tools/perf/`.

| Op (KAT-0) | C ref instr | Rust instr | Rust/C | C cycles | Rust cycles |
|---|---|---|---|---|---|
| verify | 49.3 M | 52.1 M | **1.06** | 20.5 M | 23.7 M¹ |
| keygen | 682.1 M | 339.8 M | **0.50** | 287.3 M | 161.9 M |
| sign | 1440.9 M | 784.5 M | **0.54** | 551.7 M | 339.6 M |

¹ Cycle counts are from single runs under ~100-load-average (cargo-mutants -j16
running concurrently); treat as ±15%.

**Verify** (no DRBG, no GMP): Rust executes 6% more instructions; cycles roughly
at parity. The remaining gap is the Fp2-temp pattern that LLVM mostly but not
fully eliminates.

**Keygen/sign**: Rust ~2× faster — but this is **AES-NI vs bitsliced software
AES**, not algorithm-level. `perf record` shows the C `ref` build spends 52% of
keygen in `br_aes_ct64_*` (BearSSL constant-time bitsliced AES, since `ref`
deliberately avoids CPU-specific code), while the Rust `aes` crate auto-detects
AES-NI and the DRBG cost is ~0%. C's `broadwell` build type would also use
AES-NI. Excluding AES, the GMP-heavy portion is comparable (both use fat-binary
GMP with runtime CPU dispatch).

The earlier table that showed Rust 1.3-1.5× **slower** was wrong — it compared
Rust on KAT-0 against C `benchmark_lvl1` medians over random instances. Keygen
and sign have high variance (C's own benchmark: keygen stddev 70.6 Mcyc on a
117.7 Mcyc median), and KAT-0 happens to be ~2.4× harder than median.

**Wall-clock** (`benches/bench.rs`, release, lvl1, KAT-0, light load):
verify 3.2 ms · keygen 24.9 ms · sign 58.5 ms · 100× full-KAT (kg+sign+verify) 8.0 s.

Optimizations applied: `fp2_*_ip` in-place ops (288 sites) and DRBG key-schedule
caching (~1% combined improvement; LTO already eliminated most aliasing temps).
`#[inline(always)]` on the modarith backend was tried and reverted (I-cache bloat).

Open opportunities:
- Port the C `broadwell` AVX2 `fp_*` backend (largest remaining win for verify).

## Mutation testing (verify path)

`cargo mutants` configured in `.cargo/mutants.toml` — scoped to `mp/gf/ec/hd/verification` with `--lib` only. Sign-only/dead/debug code is excluded by name there with rationale comments.

| Run | Mutants | Missed | Caught | Notes |
|---|---|---|---|---|
| Initial (`--lib`) | 1681 | 279 | 1265 | 75% |
| + targeted unit tests | 1606 | 142 | 1317 | 82% |
| + diverse-KAT lib test, equiv-mut excludes | 1461 | 103 | 1235 | 93% (49 unviable, 74 timeout = effectively caught) |

Tests added (15): `mp::{parity, digit_nonzero_ct_top_bit_only, multiple_shiftl_exact_radix}`, `ec::{validity_predicates, torsion_predicates, biscalar_mul_kbits1, dbl_iter_normalize_threshold, jac_predicates_and_init, singular_isogeny_consistency, basis_from_bad_hint_rejected, xdbladd_nonnormalized_a24}`, `verification::{reject_malformed_inputs, hash_to_challenge_stable, verify_diverse_kats}`.

**Real bug found via mutation analysis**: `ec_is_basis_four_torsion` used Rust bitwise `!` on a u32 where C uses logical `!`; on a degenerate basis with P=Q it returned `0xFFFFFFFE` (truthy) instead of 0. Unreachable on the verify path (only called when `pow_dim2_deg_resp == 0`, never in lvl1 KATs) but a real divergence from C. Fixed.

Residual survivors fall into:
- **Equivalent mutations** (~20): bit-disjoint `^`↔`|` in funnel shifts/carry-OR; extra Newton iterations in `mp_inv_2e`; `hash_to_challenge` mask is always all-ones (`2·SECURITY_BITS` is a multiple of 64 at every level).
- **Error-path validity checks in `hd/theta_isogenies.rs`** (~40): `verify_two_torsion`, `gluing_compute`, `theta_isogeny_compute`, `splitting_compute` zero-coordinate guards. Valid signatures never trigger these; tampering tends to fail upstream first. Would need a fuzzer-found invariant-violating sig to exercise.
- **Strategy/allocation tuning** in `theta_chain_compute_impl`/`ec_eval_even_strategy` (~15): mutations change the doubling/isogeny tradeoff but produce the same correct result via a different (slower) path.

## Design decisions

- Feature flags `lvl1`/`lvl3`/`lvl5` select the prime; only one active per build.
- Feature flag `sign` gates the GMP/quaternion dependency (verify-only builds are GMP-free, matching C's `ENABLE_SIGN=OFF`).
- Field element internal repr matches C's auto-generated Montgomery form for bit-exact KAT compatibility.
