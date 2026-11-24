"""Paper 1: numerical companion to Section "The case n=4".

Proposition (lower bound, n=4).  For U in U(4) and each partition
{i,j} | {k,l} of {1,2,3,4},

    |U_ij U_ji| + |U_kl U_lk|  >=  | |U_ii U_jj| - |U_kk U_ll| |,

which follows from Jacobi's complementary-minor identity

    U_ii U_jj - U_ij U_ji = det(U) * conj(U_kk U_ll - U_kl U_lk).

Proposition (zero set, n=4).  mu(d) = 0 exactly at (1,1,1,1) and on the
permutations of (0,0,0,s), 0 <= s <= 1.

The proofs in the paper are complete without this script.  It checks:

  1. the identity and the three inequalities on Haar-random unitaries;
  2. the family U_theta: orthogonality, diagonal (0,0,0,cos theta), zero
     objective;
  3. the value of the bound at the von Mangoldt diagonal (1/2,1/3,1/5,1/7),
     its feasibility (Thompson) and the bound on rho;
  4. the reduction of the bound to the n=3 closed form on the face d_4 = 1;
  5. the algebra of the zero set of the bound on a grid.
  6. the circulant family of the constant-diagonal upper bound nu(t).

Run:  python scripts/paper1_n4_check.py      (exit code 0 = all pass)
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

PARTS = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]


def haar(n, rng):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    ph = np.diag(r) / np.abs(np.diag(r))
    return q * ph


def objective(u):
    n = u.shape[0]
    return sum(abs(u[i, j] * u[j, i]) for i in range(n) for j in range(i + 1, n))


def lam(a):
    return sum(abs(a[i] * a[j] - a[k] * a[l]) for i, j, k, l in PARTS)


def main():
    rng = np.random.default_rng(2026)
    ok = True

    # 1. identity and inequalities
    err_id, slack_pair, slack_sum = 0.0, np.inf, np.inf
    for _ in range(20000):
        u = haar(4, rng)
        det = np.linalg.det(u)
        a = np.abs(np.diag(u))
        for i, j, k, l in PARTS:
            lhs = u[i, i] * u[j, j] - u[i, j] * u[j, i]
            rhs = det * np.conj(u[k, k] * u[l, l] - u[k, l] * u[l, k])
            err_id = max(err_id, abs(lhs - rhs))
            pair = abs(u[i, j] * u[j, i]) + abs(u[k, l] * u[l, k])
            slack_pair = min(slack_pair, pair - abs(a[i] * a[j] - a[k] * a[l]))
        slack_sum = min(slack_sum, objective(u) - lam(a))
    print(f"1. Jacobi identity: max error {err_id:.2e}")
    print(f"   pair inequality: min slack {slack_pair:.3e}; sum: min slack {slack_sum:.3e}")
    ok &= err_id < 1e-12 and slack_pair > -1e-12 and slack_sum > -1e-12

    # 2. the family U_theta
    worst = 0.0
    for th in np.linspace(0.0, np.pi / 2, 181):
        c, s = np.cos(th), np.sin(th)
        u = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [c, 0, 0, s], [-s, 0, 0, c]], float)
        worst = max(worst, np.abs(u @ u.T - np.eye(4)).max(),
                    np.abs(np.diag(u) - np.array([0, 0, 0, c])).max(), objective(u))
    print(f"2. U_theta: max defect (orthogonality, diagonal, objective) {worst:.2e}")
    ok &= worst < 1e-14

    # 3. von Mangoldt diagonal
    d = np.array([1 / 2, 1 / 3, 1 / 5, 1 / 7])
    a = np.sqrt(d)
    val = lam(a)
    explicit = (1 / np.sqrt(6) - 1 / np.sqrt(35)) + (1 / np.sqrt(10) - 1 / np.sqrt(21)) \
        + (1 / np.sqrt(14) - 1 / np.sqrt(15))
    thompson = a.sum() - 2 * a.min()
    rho = 2 * val / (4 - d.sum())
    print(f"3. von Mangoldt n=4: lambda = {val:.10f} (explicit {explicit:.10f}), "
          f"Thompson {thompson:.4f} <= 2, rho >= {rho:.4f}")
    ok &= abs(val - explicit) < 1e-14 and thompson <= 2 and 0.3462 < val < 0.3463 and rho > 0.245

    # 4. face d_4 = 1
    worst = 0.0
    for _ in range(2000):
        a3 = rng.uniform(size=3)
        closed = abs(a3[0] * a3[1] - a3[2]) + abs(a3[0] * a3[2] - a3[1]) + abs(a3[1] * a3[2] - a3[0])
        worst = max(worst, abs(lam(np.append(a3, 1.0)) - closed))
    print(f"4. face d_4 = 1: max |lambda - closed form| {worst:.2e}")
    ok &= worst < 1e-14

    # 5. zero set of lambda on a grid: constant or at most one non-zero entry
    grid = np.linspace(0, 1, 11)
    bad = 0
    for a in itertools.product(grid, repeat=4):
        a = np.array(a)
        if lam(a) < 1e-12:
            constant = np.ptp(a) < 1e-12
            one = (a > 1e-12).sum() <= 1
            bad += not (constant or one)
    print(f"5. zeros of lambda on an 11^4 grid outside the stated set: {bad}")
    ok &= bad == 0

    # 6. circulant upper bound on constant diagonals
    def nu(t):
        s = np.sqrt(t)
        return 2 * t if t <= 0.25 else 2 * (1 - s) * (3 * s - 1)

    worst = 0.0
    for t in np.linspace(0.0, 1.0, 201):
        s = np.sqrt(t)
        c_, x2 = (0.0, 2 * t) if t <= 0.25 else (2 * s - 1, s)
        e = np.exp(1j * np.arccos(c_))
        x, y = np.sqrt(x2), np.sqrt(1 - x2)
        u = [x / 2 * (1 + e), 1j * y / 2 * (1 - 1j * e), x / 2 * (1 - e), 1j * y / 2 * (1 + 1j * e)]
        m = np.array([[u[(k - j) % 4] for k in range(4)] for j in range(4)])
        worst = max(worst, np.abs(m @ m.conj().T - np.eye(4)).max(),
                    np.abs(np.abs(np.diag(m)) ** 2 - t).max(), abs(objective(m) - nu(t)))
    print(f"6. circulant family: max defect (unitarity, diagonal, objective - nu) {worst:.2e}; "
          f"nu(4/9) = {nu(4 / 9):.6f}")
    ok &= worst < 1e-13 and abs(nu(4 / 9) - 2 / 3) < 1e-14

    print("ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
