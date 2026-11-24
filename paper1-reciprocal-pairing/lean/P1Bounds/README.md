# P1Bounds — Lean 4 formalization of two unitary-matrix inequalities

Machine-checked (Lean 4 + Mathlib, no `sorry`, no extra axioms, no `native_decide`)
proofs of the lower bounds used in `papers/paper1_mu_note`:
the "Lower bound" part of Theorem `thm:closed` (n = 3) and Proposition `prop:lower4` (n = 4).

All code is in `P1Bounds/Basic.lean`, namespace `P1Bounds`. Indices are 0-based (`Fin 3`, `Fin 4`);
`U : Matrix.unitaryGroup (Fin n) ℂ`, `U.val` is the underlying matrix, `‖·‖` is the complex modulus.

## What is proved and how it maps to the paper

| Lean theorem | Statement | Paper |
|---|---|---|
| `cofactor_identity_three` | `i ≠ j → i ≠ k → j ≠ k →` `U i i * U j j - U i j * U j i = det U * conj (U k k)` | eq. `eq:cofactor` |
| `pair_lower_bound_three` | same hypotheses: `\|‖U i i‖ * ‖U j j‖ - ‖U k k‖\| ≤ ‖U i j * U j i‖` | Thm `thm:closed`, "More precisely…" inequality |
| `sum_lower_bound_three` | `\|ab - c\| + \|ac - b\| + \|bc - a\| ≤ ‖U01 U10‖ + ‖U02 U20‖ + ‖U12 U21‖`, with `a,b,c = ‖U00‖,‖U11‖,‖U22‖` | Thm `thm:closed`, lower bound `Σ_{i<j} ≥ …` |
| `jacobi_two_by_two_four` | `i,j,k,l` pairwise distinct `→ U i i * U j j - U i j * U j i = det U * conj (U k k * U l l - U k l * U l k)` | eq. `eq:jacobi4` |
| `pair_lower_bound_four` | same hypotheses: `\|‖U i i‖ * ‖U j j‖ - ‖U k k‖ * ‖U l l‖\| ≤ ‖U i j * U j i‖ + ‖U k l * U l k‖` | Prop `prop:lower4`, first display |
| `sum_lower_bound_four` | `\|a0 a1 - a2 a3\| + \|a0 a2 - a1 a3\| + \|a0 a3 - a1 a2\| ≤` sum of the six `‖U i j * U j i‖`, `i < j`, written out; `a_i = ‖U i i‖` | Prop `prop:lower4`, `μ ≥ λ` (pointwise in `U`) |
| `sum_lower_bound_four'` | the same with the right-hand side as a `Finset` sum over `{p : Fin 4 × Fin 4 \| p.1 < p.2}` | idem |

Supporting lemmas (over an arbitrary commutative ring `R`, hypothesis `M * N = 1`):

- `adjugate_eq_det_smul_of_mul_eq_one`: `adjugate M = det M • N` (any finite index type).
- `cofactor_three_of_mul_eq_one`: 3×3 principal 2×2 minor `= det M * N k k`.
- `jacobi_block_four`: Jacobi's complementary-minor identity for the partition `{0,1} ⊔ {2,3}`,
  by the standard determinant trick `M * [[1, N₁₂],[0, N₂₂]] = [[M₁₁, 0],[M₂₁, 1]]`.
- `jacobi_four_of_mul_eq_one`: the same for any pairwise distinct `i j k l`, by reindexing with the
  permutation `![i, j, k, l]`.
- `norm_det_unitary`: `‖det U‖ = 1`; `abs_norm_sub_norm_le_of_unit`: reverse triangle inequality.

The proofs follow the paper's arguments step by step (`adj U = det U · U*`, `|det U| = 1`,
reverse triangle inequality, sum over the three partitions).

## What is NOT formalized

- The definition of `μ(d)`, the fibre `F_d` and the feasible set; the theorems are pointwise
  statements about an arbitrary unitary `U` (the paper's bounds on `F_d` follow by substituting
  `‖U i i‖ = a_i`).
- Attainment for n = 3 (the Euler–Rodrigues/quaternion rotation), hence the equality
  `μ(d) = |ab-c| + |ac-b| + |bc-a|`.
- The zero set of `μ`/`λ` (n = 3 and n = 4) and the feasibility (Thompson/Horn) conditions.
- Anything for n ≥ 5.

## Versions

- Lean toolchain: `leanprover/lean4:v4.37.0-rc1`
  (Lean 4.37.0-rc1, commit `0b77e12952ad39077d80d8bd6b66d601f225c321`)
- Mathlib: tag `v4.37.0-rc1`, commit `31ce06eaadf4a200c0e1407c5f10953cf606b8f8`
  (pinned in `lake-manifest.json`)

## How to reproduce (Docker; nothing is installed on the host)

From this folder (PowerShell or bash; adjust the host path):

```sh
docker run -d --name p1lean --cpus 2 --memory 6g ubuntu:24.04 sleep infinity
docker exec p1lean bash -c "apt-get update && apt-get install -y curl git ca-certificates"
docker exec p1lean bash -c "curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none"
docker exec p1lean mkdir -p /root/P1Bounds
docker cp . p1lean:/root/P1Bounds
# fetch prebuilt Mathlib (do NOT compile Mathlib from source), then build
docker exec p1lean bash -c "cd /root/P1Bounds && ~/.elan/bin/lake exe cache get && ~/.elan/bin/lake build"
docker exec p1lean bash -c "cd /root/P1Bounds && grep -n 'sorry\|axiom\|native_decide' *.lean P1Bounds/*.lean"
# clean up
docker rm -f p1lean
```

The `.lake` directory (Mathlib, ~6 GB with the toolchain) stays inside the container; only the
sources, `lakefile.toml`, `lake-manifest.json` and `lean-toolchain` live here.

## Verification output (2026-11-23)

`lake build` (the repeated notice "imports `Mathlib` … consider adding `module`" is filtered out;
it is the only warning):

```
⚠ [18291/18295] Built P1Bounds.Basic (19s)
info: P1Bounds/Basic.lean:262:0: 'P1Bounds.cofactor_identity_three' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:263:0: 'P1Bounds.pair_lower_bound_three' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:264:0: 'P1Bounds.sum_lower_bound_three' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:265:0: 'P1Bounds.jacobi_two_by_two_four' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:266:0: 'P1Bounds.pair_lower_bound_four' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:267:0: 'P1Bounds.sum_lower_bound_four' depends on axioms: [propext, Classical.choice, Quot.sound]
info: P1Bounds/Basic.lean:268:0: 'P1Bounds.sum_lower_bound_four'' depends on axioms: [propext, Classical.choice, Quot.sound]
⚠ [18292/18295] Ran P1Bounds.Basic:irArts
⚠ [18293/18295] Built P1Bounds (5.9s)
⚠ [18294/18295] Ran P1Bounds:irArts
Build completed successfully (18295 jobs).
```

`grep -n "sorry\|axiom\|native_decide" *.lean P1Bounds/*.lean` — the only hits are the seven
`#print axioms` audit commands themselves (the word "axioms"); there is no `sorry`, no `axiom`
declaration and no `native_decide`:

```
P1Bounds/Basic.lean:262:#print axioms P1Bounds.cofactor_identity_three
P1Bounds/Basic.lean:263:#print axioms P1Bounds.pair_lower_bound_three
P1Bounds/Basic.lean:264:#print axioms P1Bounds.sum_lower_bound_three
P1Bounds/Basic.lean:265:#print axioms P1Bounds.jacobi_two_by_two_four
P1Bounds/Basic.lean:266:#print axioms P1Bounds.pair_lower_bound_four
P1Bounds/Basic.lean:267:#print axioms P1Bounds.sum_lower_bound_four
P1Bounds/Basic.lean:268:#print axioms P1Bounds.sum_lower_bound_four'
```

One deviation from the default `lake new … math` template: the line
`weak.linter.mathlibStandardSet = true` was removed from `lakefile.toml` (Mathlib's style linters
only; no effect on soundness).
