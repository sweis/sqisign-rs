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
| quaternion (adv) | ~2200 | in progress | HNF/LLL/lattice/ideal/normeq; clean-room DPE |
| precomp (sign) | ~263K | pending | quaternion_data, endomorphism_action |
| ec (biextension) | 770 | pending | pairings, sign-only |
| id2iso | 2594 | pending | needs quaternion, hd |
| signature | 1370 | stub | needs id2iso |

## Test vectors

- KAT files copied to `KAT/` from C repo (100 vectors per level).
- `tests/kat.rs` parses .rsp files; verify/sign tests are `#[ignore]` until impl lands.

## C oddities flagged for spec cross-check

These are ported faithfully but look suspicious or underspecified:

- `mp/mp_mul2` omits the `a[1]*b[0]` partial product.
- `mp/mp_neg` does not propagate the +1 carry past limb 0.
- `mp/mp_shiftr` returns only `x[0]&1` regardless of shift amount.
- `ec/ec_is_equal` uses C precedence `~l_zero & ~r_zero * lr_equal` (mul binds tighter than and).
- `ec/xMUL` calls `xDBLADD(...,true)` even when A24 was computed with z≠1.
- `ec/jac_from_ws` does not copy `Q->x` when A=0 — dead path in current callers but a latent bug.
- `quaternion/ibz_rand_interval_bits` subtracts the bit-count `m` from the result — likely a C bug.
- `quaternion/ibz_div` doc-comment says floor but impl is truncating (`mpz_tdiv_qr`).

## Spec discrepancies found

(spec PDF cross-check pending)

## Design decisions

- Feature flags `lvl1`/`lvl3`/`lvl5` select the prime; only one active per build.
- Feature flag `sign` gates the GMP/quaternion dependency (verify-only builds are GMP-free, matching C's `ENABLE_SIGN=OFF`).
- Field element internal repr matches C's auto-generated Montgomery form for bit-exact KAT compatibility.
