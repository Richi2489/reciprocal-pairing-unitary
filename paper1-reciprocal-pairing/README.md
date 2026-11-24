# Paper 1 — Reciprocal Pairing in Unitary Matrices with Prescribed Diagonal Moduli

## Version 2 of the paper (November 2026)

The manuscript was rewritten (15 pages). Its statements are checked by
`scripts/paper1_v2_check.py` (run it from this folder; exit code 0 = all pass):
the one-invariant lemma and the closed form in dimension three, the
factorization identity and the maximum of M on the fibre, the inequality for
symmetric unitary 3x3 matrices, the n=4 lower bound, zero family and circulant
bound, the Paley/doubly-regular-tournament zeros and the count of doubly regular
tournaments of order 7, the splitting lemma, the order-11 matrix and its
permutation form, and the merging theorem with the decay bounds for reciprocal
primes.
`scripts/paper1_const4_search.py` is the second optimizer used as evidence for
the circulant conjecture. `lean/P1Bounds/` is the Lean formalization.
The other scripts and the reports in `docs/` belong to the previous, longer
version; they remain valid, but their docstrings use the old numbering of
theorems, tables and problems.


Replication package. Contents:

- `qgraphs/unistochastic.py` — the exact one-variable reduction at n=3
  (the unistochastic sub-interval via the chain-link condition, and the
  explicit construction/verification of a unitary at a given point).
  Its `barrer` is the original numerical sweep of that reduction: a fine
  grid plus local refinement, accurate to grid resolution, returning
  `None` on infeasible diagonals. It is kept as the historical tool the
  early tables came from. The exact procedures are the closed form of
  Theorem 1 (`scripts/paper1_closed_form_check.py`) and the quartic root
  solver for the boundary of K (`scripts/paper1_feasible_set_check.py`).
- `scripts/run_mu_scaling.py` — the SLSQP multi-start scaling study
  (n = 3..10; stable optimizer range n <= 6, values feasible to tolerance,
  not certified upper bounds), with the decision rules printed before
  computation.
- `scripts/paper1_closed_form_check.py` — Theorem 1, the n=3 closed form
  mu(d) = |ab-c| + |ac-b| + |bc-a| (a, b, c = square roots of the diagonal):
  builds the quaternion rotation that attains it, checks it against the
  quartic solver and Haar-random unitaries, and checks the Hoelder-1/2
  continuity bound (Corollary 1).
- `scripts/paper1_feasible_set_check.py` — checks Theorem 2 and Table 1
  (n=3 minima, components of the feasible set K, the exact merger at c=1/9,
  the Hellinger identity) and the doubly-stochastic vs. unitary feasibility
  condition against direct construction.
- `scripts/paper1_thompson_groundtruth.py` — ground truth for the unitary
  feasibility condition: for each d it tries to construct U(3) with
  |U_ii|^2 = d directly and compares with
  sqrt(d_i) + sqrt(d_j) <= 1 + sqrt(d_k) (no grids).
- `qgraphs/mu_general.py` — the SLSQP engine the scaling study imports.
- `scripts/paper1_extra_checks.py` — Problem O4: exhaustive enumeration of
  oriented graphs for n <= 6 and the Paley construction at n = 7.
- `scripts/paper1_n4_check.py` — the n=4 section: Jacobi identity and the
  lower bound on Haar-random U(4), the zero family U_theta, the value at
  the von Mangoldt diagonal, the circulant bound (Propositions 4-6, Corollary 2).
- `lean/P1Bounds/` — Lean 4 + Mathlib formalization of the n=3 and n=4 lower
  bounds (no `sorry`; see its README for the build and `#print axioms`).
- `data/mu_scaling_results.csv` — output table used in the paper.
  CAUTION: `python -m scripts.run_mu_scaling` with no argument reruns the
  whole study and rewrites this file and `docs/mu_scaling_report.md` as it
  goes; `--report-only` rebuilds report and figures from the CSV, and any
  other argument is refused.
- `docs/mu_scaling_prereg.md` — the frozen decision rules of the scaling
  study; `docs/mu_ruta2_budget.md` (budget bound) and
  `docs/next_step_prereg.md` are the notes it refers to.
- `docs/mu_scaling_report.md`, `docs/mu_d_n3_verification.md` — the
  auto-generated reports the paper's tables are read from. The internal
  reports, plans and some script comments are in Spanish.

Reproduce: `python -m scripts.run_mu_scaling` (scaling study) and
`python scripts/paper1_feasible_set_check.py` (Theorem 2, Table 1) and
`python scripts/paper1_closed_form_check.py` (Theorem 1, Corollary 1) and
`python scripts/paper1_n4_check.py` (the n=4 section) from this folder's root
(Python 3.11+, `pip install -r ../requirements.txt`).

Exact values quoted in the paper (Table 1, the constant-family closed form
mu(c,c,c) = 3(sqrt(c)-c), and the corner asymptotics) can be re-derived with
`qgraphs/unistochastic.py`; see the paper's Section 3 for the algebra.
