"""Paper 1: multi-start search of mu(t,t,t,t) over the whole fibre in U(4).

Second implementation, independent of qgraphs/mu_general.py, used as
numerical evidence for the conjecture mu(t,t,t,t) = nu(t). SLSQP over
U = exp(A), A skew-Hermitian, with a smoothed objective and hard diagonal
constraints (accepted if the diagonal error is below 1e-8). Outputs are
best values found, not certified minima.

Run:  python scripts/paper1_const4_search.py const        (16 values of t; slow)
      python scripts/paper1_const4_search.py show d1 d2 d3 d4
"""
import sys
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

N = 4
IU = np.triu_indices(N, 1)


def U_of(x):
    A = np.zeros((N, N), complex)
    k = 0
    for i in range(N):
        A[i, i] = 1j * x[k]; k += 1
    for i in range(N):
        for j in range(i + 1, N):
            z = x[k] + 1j * x[k + 1]; k += 2
            A[i, j] = z; A[j, i] = -np.conj(z)
    return expm(A)


def pairs(U):
    return np.abs(U[IU] * U.T[IU])


def mu_opt(d, starts=40, seed=0, eps=1e-10):
    d = np.asarray(d, float)
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    cons = {"type": "eq", "fun": lambda x: np.abs(np.diag(U_of(x))) ** 2 - d}
    for _ in range(starts):
        x0 = rng.normal(scale=1.5, size=N * N)
        for e in (1e-4, 1e-7, eps):
            f = lambda x: np.sum(np.sqrt(np.abs(U_of(x)[IU] * U_of(x).T[IU]) ** 2 + e * e))
            r = minimize(f, x0, method="SLSQP", constraints=[cons],
                         options={"maxiter": 400, "ftol": 1e-13})
            x0 = r.x
        U = U_of(x0)
        err = np.abs(np.abs(np.diag(U)) ** 2 - d).max()
        val = pairs(U).sum()
        if err < 1e-8 and val < best[0]:
            best = (val, U)
    return best


def lam(d):
    a = np.sqrt(d)
    return abs(a[0]*a[1]-a[2]*a[3]) + abs(a[0]*a[2]-a[1]*a[3]) + abs(a[0]*a[3]-a[1]*a[2])


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True, linewidth=150)
    mode = sys.argv[1]
    if mode == "const":
        for t in [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
            v, U = mu_opt([t] * 4, starts=25)
            print(f"t={t:.2f} mu={v:.6f} mu/t={v/t:.4f} budget={2-2*t:.3f} pairs={np.sort(pairs(U))}", flush=True)
    if mode == "show":
        d = [float(x) for x in sys.argv[2:6]]
        v, U = mu_opt(d, starts=40)
        print("mu", v, "lambda", lam(d)); print("|U|\n", np.abs(U)); print("pairs", pairs(U)); print("det", np.linalg.det(U))
