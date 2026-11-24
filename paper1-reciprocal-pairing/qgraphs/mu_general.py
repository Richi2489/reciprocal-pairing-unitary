"""
mu_general.py — μ(d) para n ≥ 4 por optimización directa sobre U(n)
==================================================================

En n=3 la reducción unistocástica hace μ(d) un mínimo exacto de 1 variable
(`unistochastic.barrer`). Para n≥4 el conjunto unistocástico con diagonal fija no
tiene descripción por desigualdades simples, así que minimizamos M(U) directamente
sobre la fibra {U∈U(n): |U_ii|²=d_i}.

    μ(d) = mín_{U∈U(n)} M(U),   M1(U) = Σ_{i<j} |U_ij·U_ji|,   sujeto a |U_ii|² = d_i.

Parametrización: U = expm(A), A anti-hermitiana (n² parámetros reales). exp es
suprayectiva sobre U(n), así que todo unitario es alcanzable y el multi-start cubre
el paisaje. Restricciones de diagonal duras vía SLSQP.

d NO se normaliza: d_i = |U_ii|² ∈ [0,1] independientes (para von Mangoldt d_i=1/p_i).

Ver `docs/mu_scaling_prereg.md` para las reglas pre-registradas.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

# ---------------------------------------------------------------------------
# Familias de diagonales (SIN normalizar)
# ---------------------------------------------------------------------------


def _primos(n):
    """Los primeros n primos."""
    ps, c = [], 2
    while len(ps) < n:
        if all(c % p for p in ps):
            ps.append(c)
        c += 1
    return ps


def von_mangoldt_diag(n, power=1.0):
    """d_i = 1/p_i^power, p_i el i-ésimo primo. Diagonal objetivo (crudo)."""
    return np.array([1.0 / p ** power for p in _primos(n)], dtype=float)


def constant_diag(n, c=0.4):
    """d_i = c constante. Σd_i = c·n (crece con n: ojo confounder)."""
    return np.full(n, float(c))


def constant_diag_fixed_sum(n, s0=2.0):
    """d_i = s0/n constante, así Σd_i = s0 FIJO. Control C2 del confounder."""
    return np.full(n, float(s0) / n)


def linear_diag(n, dmax=0.5):
    """d_i ∝ (n-i+1), escalada a max d_i = dmax (rango fijo, indep. de n)."""
    v = np.arange(n, 0, -1, dtype=float)
    return v / v.max() * dmax


def random_diag(n, seed=0):
    """d_i ~ Uniforme(0,1) i.i.d. SIN normalizar."""
    return np.random.default_rng(seed).uniform(0.0, 1.0, size=n)


# ---------------------------------------------------------------------------
# Factibilidad (prefiltro necesario de doblemente-estocástica)
# ---------------------------------------------------------------------------


def diagonal_feasible_ds(d, tol=1e-12):
    """Condición NECESARIA (no suficiente) para que d sea la diagonal de una
    doblemente estocástica: para todo i, (1-d_i) ≤ Σ_{j≠i}(1-d_j), i.e.
    Σd ≤ 2 d_i + n - 2. Reproduce la infactibilidad de (0.8,0.5,0.2) en n=3.
    Unistocasticidad (n≥4) es más estricta; el optimizador es el árbitro final."""
    d = np.asarray(d, float)
    n = len(d)
    if np.any(d < -tol) or np.any(d > 1 + tol):
        return False
    S = d.sum()
    return bool(np.all(S <= 2.0 * d + (n - 2) + tol))


# ---------------------------------------------------------------------------
# Parametrización U = expm(A), A anti-hermitiana
# ---------------------------------------------------------------------------


def _n_params(n):
    return n * n


def antiherm_from_params(x, n):
    """Construye A anti-hermitiana (A†=-A) desde n² reales:
    n en la diagonal imaginaria + n(n-1)/2 pares (re,im) en el triángulo superior."""
    A = np.zeros((n, n), dtype=complex)
    A[np.diag_indices(n)] = 1j * x[:n]
    k = n
    iu = np.triu_indices(n, 1)
    for a, b in zip(*iu):
        re, im = x[k], x[k + 1]
        A[a, b] = re + 1j * im
        A[b, a] = -re + 1j * im  # anti-hermitiana: A_ba = -conj(A_ab)
        k += 2
    return A


def unitary_from_params(x, n):
    return expm(antiherm_from_params(x, n))


def M1(U):
    """Masa emparejada cruda Σ_{i<j} |U_ij·U_ji|."""
    n = U.shape[0]
    iu = np.triu_indices(n, 1)
    return float(np.sum(np.abs(U[iu] * U[(iu[1], iu[0])])))


def M2(U):
    """Σ_{i<j} |U_ij·U_ji|² (sonda de robustez, NO es μ)."""
    n = U.shape[0]
    iu = np.triu_indices(n, 1)
    return float(np.sum(np.abs(U[iu] * U[(iu[1], iu[0])]) ** 2))


_OBJ = {"M1": M1, "M2": M2}


# ---------------------------------------------------------------------------
# Minimización con restricción dura de diagonal
# ---------------------------------------------------------------------------


def _constraint_err_max(x, n, d):
    U = unitary_from_params(x, n)
    return float(np.max(np.abs(np.abs(np.diag(U)) ** 2 - d)))


def min_mu_unitary(d, objective="M1", num_starts=200, seed=0,
                   tol_feas=1e-6, maxiter=200, check_feasible=True):
    """Minimiza M(U) sobre {U∈U(n): |U_ii|²=d_i} con SLSQP multi-start.

    Devuelve dict con:
      mu           : mejor M factible (o None si infactible)
      U            : minimizador
      err          : error de restricción del mejor (max_i||U_ii|²-d_i|)
      feasible     : bool (algún arranque logró err<tol_feas)
      infeasible_apriori : bool (falla el prefiltro DS necesario)
      mins         : lista de M de TODOS los arranques factibles (para percentiles)
      errs         : errores de restricción por arranque
      n_starts_feasible : cuántos arranques quedaron factibles
    """
    d = np.asarray(d, float)
    n = len(d)
    obj = _OBJ[objective]
    out = {"mu": None, "U": None, "err": None, "feasible": False,
           "infeasible_apriori": False, "mins": [], "errs": [],
           "n_starts_feasible": 0, "n": n, "sum_d": float(d.sum())}

    if check_feasible and not diagonal_feasible_ds(d):
        out["infeasible_apriori"] = True
        return out

    rng = np.random.default_rng(seed)
    npar = _n_params(n)
    fun = lambda x: obj(unitary_from_params(x, n))
    cons = [{"type": "eq",
             "fun": lambda x: np.abs(np.diag(unitary_from_params(x, n))) ** 2 - d}]

    best_m, best_U, best_err = np.inf, None, np.inf
    for _ in range(num_starts):
        x0 = rng.standard_normal(npar) * 1.5
        try:
            res = minimize(fun, x0, method="SLSQP", constraints=cons,
                           options={"maxiter": maxiter, "ftol": 1e-10})
        except Exception:
            continue
        err = _constraint_err_max(res.x, n, d)
        out["errs"].append(err)
        if err < tol_feas:
            U = unitary_from_params(res.x, n)
            m = obj(U)
            out["mins"].append(m)
            out["n_starts_feasible"] += 1
            if m < best_m:
                best_m, best_U, best_err = m, U, err

    if out["n_starts_feasible"] > 0:
        out.update(mu=best_m, U=best_U, err=best_err, feasible=True)
    return out
