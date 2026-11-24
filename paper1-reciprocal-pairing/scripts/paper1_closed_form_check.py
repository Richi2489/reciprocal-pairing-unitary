"""Paper 1: numerical companion to the n=3 closed form (Theorem 1, Corollary 1).

Theorem 1.  For feasible d = (d1,d2,d3), with a=sqrt(d1), b=sqrt(d2), c=sqrt(d3),

    mu(d) = |ab - c| + |ac - b| + |bc - a|,

attained by the quaternion rotation R with
    q0^2 = (1+a+b+c)/4, q1^2 = (1+a-b-c)/4, q2^2 = (1-a+b-c)/4, q3^2 = (1-a-b+c)/4.

The proof in the paper is complete without this script.  The script checks, to
machine precision:

  1. the construction: R R^T = I, det R = 1, diag R = (a,b,c), the cofactor
     identity R_ij R_ji = R_ii R_jj - R_kk, and objective(R) = closed form, on
     random feasible d, on the Thompson boundary, and at special points
     (0,0,0), (1,1,1), (1,0,0) and permutations, the footnote diagonal
     (0.04699, 0.59058, 0), and the rows of Table 1;
  2. feasibility: the four q_i^2 are >= 0 exactly when the Thompson condition
     (2) holds (so infeasible diagonals such as (1/2,1/2,0) and (0.8,0.5,0.2)
     get a negative q_i^2);
  3. agreement with the independent quartic solver of
     scripts/paper1_feasible_set_check.py (minimum of W over the outermost
     points of K) on interior d;
  4. the lower bound on Haar-random unitaries (never violated) and the complex
     cofactor identity U_ii U_jj - U_ij U_ji = det(U) conj(U_kk);
  5. the Hoelder estimate of Corollary 1 on random feasible pairs, and the
     sharpness witness mu(t,t,t) = 3(sqrt t - t).

Run:  python scripts/paper1_closed_form_check.py      (exit code 0 = all pass)
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paper1_feasible_set_check as quartic  # noqa: E402

PAIRS = ((0, 1, 2), (0, 2, 1), (1, 2, 0))  # (i, j, k)
TOL = 1e-12


def closed(d):
    a, b, c = np.sqrt(np.asarray(d, dtype=float))
    return abs(a * b - c) + abs(a * c - b) + abs(b * c - a)


def q_squares(d):
    a, b, c = np.sqrt(np.asarray(d, dtype=float))
    return np.array([1 + a + b + c, 1 + a - b - c, 1 - a + b - c, 1 - a - b + c]) / 4.0


def thompson(d):
    r = np.sqrt(np.asarray(d, dtype=float))
    return r.sum() - 2.0 * r.min() <= 1.0 + 1e-15


def rotation(d):
    q0, q1, q2, q3 = np.sqrt(np.maximum(q_squares(d), 0.0))
    return np.array([
        [1 - 2 * (q2 * q2 + q3 * q3), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
        [2 * (q1 * q2 + q0 * q3), 1 - 2 * (q1 * q1 + q3 * q3), 2 * (q2 * q3 - q0 * q1)],
        [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 1 - 2 * (q1 * q1 + q2 * q2)],
    ])


def objective(U):
    return sum(abs(U[i, j] * U[j, i]) for i, j, _ in PAIRS)


def haar_unitary(rng, n=3):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    ph = np.diag(R) / np.abs(np.diag(R))
    return Q * ph  # Mezzadri's phase fix -> Haar measure


def feasible_sample(rng, m):
    out = []
    while len(out) < m:
        d = rng.random(3)
        if thompson(d):
            out.append(d)
    return np.array(out)


def boundary_sample(rng, m):
    """Diagonals on the Thompson boundary a+b = 1+c (and permutations)."""
    out = []
    while len(out) < m:
        a, b = rng.random(2)
        c = a + b - 1.0
        if c < 0:
            continue
        r = np.array([a, b, c])
        perm = rng.permutation(3)
        out.append(r[perm] ** 2)
    return np.array(out)


def main():
    rng = np.random.default_rng(20260924)
    fails = []

    special = [(0, 0, 0), (1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1),
               (0.04699, 0.59058, 0.0), (0.5, 1 / 3, 0.2), (0.4, 0.4, 0.4),
               (0.5, 0.45, 0.4), (0.6, 0.4, 0.4), (0.01, 0.01, 0.01),
               (1 / 9, 1 / 9, 1 / 9), (0.25, 0.25, 0.25)]
    ds = np.vstack([np.array(special, dtype=float), feasible_sample(rng, 5000),
                    boundary_sample(rng, 1000),
                    np.array([[t, t, t] for t in np.linspace(0, 1, 201)])])

    # 1. construction ------------------------------------------------------
    e_orth = e_det = e_diag = e_obj = e_cof = 0.0
    for d in ds:
        if not thompson(d):
            fails.append(("sample not feasible", d))
            continue
        R = rotation(d)
        r = np.sqrt(d)
        e_orth = max(e_orth, np.abs(R @ R.T - np.eye(3)).max())
        e_det = max(e_det, abs(np.linalg.det(R) - 1.0))
        e_diag = max(e_diag, np.abs(np.diag(R) - r).max())
        e_obj = max(e_obj, abs(objective(R) - closed(d)))
        for i, j, k in PAIRS:
            e_cof = max(e_cof, abs(R[i, j] * R[j, i] - (R[i, i] * R[j, j] - R[k, k])))
    print(f"[1] construction on {len(ds)} feasible diagonals "
          f"({len(special)} special, 1000 on the Thompson boundary, 201 constant)")
    print(f"    max |RR^T-I| = {e_orth:.2e}   max |det R - 1| = {e_det:.2e}   "
          f"max |diag R - sqrt d| = {e_diag:.2e}")
    print(f"    max |R_ijR_ji - (R_iiR_jj - R_kk)| = {e_cof:.2e}   "
          f"max |objective(R) - closed form| = {e_obj:.2e}")
    for name, v in (("orth", e_orth), ("det", e_det), ("diag", e_diag),
                    ("objective", e_obj), ("cofactor", e_cof)):
        if v > TOL:
            fails.append(("construction " + name, v))

    # 2. feasibility <=> q_i^2 >= 0 ---------------------------------------
    bad = 0
    for d in rng.random((20000, 3)):
        if thompson(d) != bool((q_squares(d) >= -1e-15).all()):
            bad += 1
    for d in [(0.5, 0.5, 0.0), (0.8, 0.5, 0.2), (1.0, 1.0, 0.0)]:
        qs = q_squares(d)
        print(f"    infeasible d = {d}: min q_i^2 = {qs.min():+.6f}  (negative)")
        if qs.min() >= 0:
            fails.append(("infeasible d has q^2 >= 0", d))
    print(f"[2] Thompson condition vs (all q_i^2 >= 0) on 20000 random d: "
          f"{bad} disagreements")
    if bad:
        fails.append(("thompson vs q", bad))

    # 3. quartic solver (outermost points of K) on interior d --------------
    e_q = 0.0
    nq = 0
    for d in feasible_sample(rng, 3000):
        if d.min() < 1e-3 or d.max() > 1 - 1e-3:
            continue
        _, comps = quartic.components_of_K(tuple(d))
        if not comps:
            fails.append(("quartic solver found K empty", d))
            continue
        xm, xp = comps[0][0], comps[-1][1]
        w = quartic.W_of_x(tuple(d), np.array([xm, xp]))
        e_q = max(e_q, float(np.abs(w - closed(d)).max()))
        nq += 1
        if nq >= 1000:
            break
    print(f"[3] quartic solver vs closed form on {nq} interior diagonals: "
          f"max |W(x_pm) - closed| = {e_q:.2e}")
    if e_q > 1e-9:
        fails.append(("quartic", e_q))

    # 4. Haar lower bound ---------------------------------------------------
    gap = np.inf
    e_ccof = 0.0
    N = 20000
    for _ in range(N):
        U = haar_unitary(rng)
        dU = np.abs(np.diag(U)) ** 2
        gap = min(gap, objective(U) - closed(dU))
        det = np.linalg.det(U)
        for i, j, k in PAIRS:
            e_ccof = max(e_ccof, abs(U[i, i] * U[j, j] - U[i, j] * U[j, i]
                                     - det * np.conj(U[k, k])))
    print(f"[4] {N} Haar unitaries: min (objective - closed form) = {gap:+.3e}; "
          f"max cofactor-identity residual = {e_ccof:.2e}")
    if gap < -TOL:
        fails.append(("Haar lower bound violated", gap))
    if e_ccof > TOL:
        fails.append(("complex cofactor identity", e_ccof))

    # 5. Hoelder estimate of Corollary 1 -------------------------------------
    worst1 = worst2 = -np.inf
    A = feasible_sample(rng, 4000)
    B = feasible_sample(rng, 4000)
    near = np.clip(A + rng.normal(scale=1e-3, size=A.shape), 0, 1)
    pairs = list(zip(A, B)) + [(x, y) for x, y in zip(A, near) if thompson(y)]
    for d, e in pairs:
        lhs = abs(closed(d) - closed(e))
        mid = 3 * np.abs(np.sqrt(d) - np.sqrt(e)).sum()
        rhs = 9 * np.sqrt(np.abs(d - e).max())
        worst1 = max(worst1, lhs - mid)
        worst2 = max(worst2, mid - rhs)
    print(f"[5] Corollary 1 on {len(pairs)} feasible pairs: "
          f"max(|dmu| - 3 sum|dsqrt|) = {worst1:+.2e}, "
          f"max(3 sum|dsqrt| - 9 ||dd||^1/2) = {worst2:+.2e}")
    if worst1 > TOL or worst2 > TOL:
        fails.append(("Hoelder", worst1, worst2))
    for t in (1e-2, 1e-4, 1e-6, 1e-8):
        ratio = closed((t, t, t)) / np.sqrt(t)
        print(f"    mu(t,t,t)/sqrt(t) at t={t:.0e}: {ratio:.6f}  (-> 3: exponent 1/2 sharp)")

    # reference values --------------------------------------------------------
    print("[ref] closed-form values:")
    for d in [(0.5, 1 / 3, 0.2), (0.4, 0.4, 0.4), (0.5, 0.45, 0.4), (0.6, 0.4, 0.4),
              (0.01, 0.01, 0.01), (0.04699, 0.59058, 0.0), (1, 0, 0), (0, 0, 0),
              (1, 1, 1)]:
        v = closed(d)
        den = 3 - sum(d)
        rho = 2 * v / den if den > 0 else float("nan")
        print(f"    d = {tuple(round(x, 6) for x in d)}: mu = {v:.10f}   rho = {rho:.10f}")
    checks = {(0.5, 1 / 3, 0.2): 0.7489956996, (0.01, 0.01, 0.01): 0.27,
              (0.04699, 0.59058, 0.0): 1.151851, (0.4, 0.4, 0.4): 0.6973665961}
    for d, v in checks.items():
        if abs(closed(d) - v) > 1e-6:
            fails.append(("reference value", d, closed(d), v))

    if fails:
        print("FAIL:", fails)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
