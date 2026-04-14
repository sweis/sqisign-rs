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

## Design decisions

- Feature flags `lvl1`/`lvl3`/`lvl5` select the prime; only one active per build.
- Feature flag `sign` gates the GMP/quaternion dependency (verify-only builds are GMP-free, matching C's `ENABLE_SIGN=OFF`).
- Field element internal repr matches C's auto-generated Montgomery form for bit-exact KAT compatibility.
