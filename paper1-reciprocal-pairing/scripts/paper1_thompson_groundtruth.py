"""Verdad de campo: para cada d, se intenta CONSTRUIR U(3) con |U_ii|^2 = d.
Se contrasta con la condicion sqrt(d_i)+sqrt(d_j) <= 1+sqrt(d_k). Sin mallas."""
import numpy as np, sys
from scipy.optimize import minimize
from scipy.linalg import expm
sys.stdout.reconfigure(encoding="utf-8")
rng = np.random.default_rng(2026)


def realizable(d, tries=25, tol=1e-12):
    for _ in range(tries):
        v0 = rng.normal(size=18) * 0.8

        def obj(v):
            A = (v[:9] + 1j*v[9:]).reshape(3, 3)
            U = expm(1j*(A + A.conj().T))
            return float(np.sum((np.abs(np.diag(U))**2 - d)**2))

        r = minimize(obj, v0, method="Nelder-Mead",
                     options=dict(maxiter=8000, maxfev=8000, fatol=1e-18, xatol=1e-13))
        if r.fun < tol:
            return True, r.fun
    return False, r.fun


def thompson3(d):
    r = np.sqrt(np.asarray(d, float))
    return all(r[i] + r[j] <= 1 + r[k] + 1e-9 for i, j, k in [(0,1,2),(0,2,1),(1,2,0)])


print("d                        Thompson  construible  residuo    veredicto")
bad = 0
N = 45
for t in range(N):
    u = rng.random()
    if u < 0.3:
        d = rng.random(3)
    elif u < 0.6:                              # cerca de la frontera de Thompson
        d = np.clip(rng.random(3) * 1.3 - 0.1, 0, 1)
    elif u < 0.8:                              # con un cero exacto
        d = rng.random(3); d[rng.integers(3)] = 0.0
    else:                                      # con un uno exacto
        d = rng.random(3); d[rng.integers(3)] = 1.0
    th = thompson3(d)
    ok, res = realizable(d)
    match = (th == ok)
    if not match:
        bad += 1
    print(f"{np.round(d,4)}  {str(th):8} {str(ok):11}  {res:.1e}  "
          f"{'ok' if match else '<<< DISCREPANCIA'}")
print()
print(f"{N} diagonales, discrepancias contra construccion directa: {bad}")
