#!/usr/bin/env python3
"""Verificacion numerica de los resultados de la Seccion "Two further results"
del Paper 1 (anadidos 2026-09-09/10).

(A) Ceros interiores de mu en dimension n = 4m+3: si T es un torneo doblemente
    regular (T T^T = (m+1) I + m J), entonces U = (T - m I)/(m+1) es ortogonal
    real, |U_ii|^2 = (m/(m+1))^2 y U_ij U_ji = 0 para todo i != j, de modo que
    mu(c,...,c) = 0 con c = (m/(m+1))^2.  Torneos de Paley (n primo = 3 mod 4).
    Sumas directas: d = (1/4) 1_n tiene mu = 0 para n = 14, 21.

(B) Cota inferior (Lema): con sigma = sum_i (1 - sqrt d_i), Delta = n - sum d_i,
        M(U) >= 2 sigma - sigma^2 - Delta/2 >= sigma (1 - sigma),
    y rho >= 1 - sigma.  Identidad usada: Delta + sum (1-sqrt d_i)^2 = 2 sigma.
    Se prueba contra los valores exactos n=3 de la familia constante y contra
    miles de unitarias aleatorias (cualquier U de la fibra acota mu(d) por
    arriba, asi que M(U) >= cota(d(U)) es condicion necesaria).

(C) Lema combinatorio: no hay ceros INTERIORES (todos los d_i en (0,1)) para
    n <= 6.  Condicion necesaria sobre el grafo orientado i->j iff U_ij != 0:
    cada arco i->j tiene un vecino de salida comun y uno de entrada comun
    (ortogonalidad de filas y de columnas), luego grados >= 2.  Enumeracion
    exhaustiva de los 3^C(n,2) grafos orientados para n = 4, 5, 6: ninguno
    cumple.  El torneo de Paley n=7 si cumple.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(20260909)


def paley_tournament(p: int) -> np.ndarray:
    assert p % 4 == 3
    qr = {(x * x) % p for x in range(1, p)}
    T = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i, j] = 1.0
    return T


def pairing(U: np.ndarray) -> float:
    A = np.abs(U * U.T)
    return float(np.sum(np.triu(A, 1)))


def sigma_delta(d: np.ndarray) -> tuple[float, float]:
    return float(np.sum(1 - np.sqrt(d))), float(len(d) - d.sum())


def bound_full(d: np.ndarray) -> float:
    s, D = sigma_delta(d)
    return 2 * s - s * s - D / 2


def bound_sigma(d: np.ndarray) -> float:
    s, _ = sigma_delta(d)
    return s * (1 - s)


# ----------------------------------------------------------------------------
def check_A() -> None:
    print("(A) ceros interiores via torneos de Paley")
    blocks = {}
    for p in (3, 7, 11, 19, 23, 31, 43):
        m = (p - 3) // 4
        T = paley_tournament(p)
        J = np.ones((p, p))
        assert np.allclose(T + T.T, J - np.eye(p))
        assert np.allclose(T @ T.T, (m + 1) * np.eye(p) + m * J), "no doblemente regular"
        U = (T - m * np.eye(p)) / (m + 1)
        orth = np.max(np.abs(U @ U.T - np.eye(p)))
        d = np.diag(U) ** 2
        c = (m / (m + 1)) ** 2
        s, D = sigma_delta(d)
        print(f"  n={p:2d} m={m:2d}  ||UU^T-I||={orth:.1e}  |U_ii|^2={d[0]:.6f} (c={c:.6f})"
              f"  M(U)={pairing(U):.1e}  Delta={D:.3f}  sigma={s:.3f}")
        assert orth < 1e-12 and np.allclose(d, c) and pairing(U) < 1e-14
        blocks[p] = U
    # sumas directas: d = 1/4 uniforme en n = 14, 21
    for k in (2, 3):
        U = np.kron(np.eye(k), blocks[7])
        d = np.diag(U) ** 2
        assert np.allclose(d, 0.25) and pairing(U) < 1e-14
        assert np.max(np.abs(U @ U.T - np.eye(7 * k))) < 1e-12
        print(f"  suma directa {k}x(n=7): n={7 * k}, d=1/4 uniforme, M(U)={pairing(U):.1e}")


# ----------------------------------------------------------------------------
def check_B_identity_and_exact_n3() -> None:
    print("(B) identidad Delta + sum(1-sqrt d)^2 = 2 sigma, y cota vs exacto n=3")
    for _ in range(1000):
        d = rng.uniform(0, 1, size=rng.integers(2, 12))
        s, D = sigma_delta(d)
        assert abs(D + np.sum((1 - np.sqrt(d)) ** 2) - 2 * s) < 1e-12
    print("  identidad verificada en 1000 vectores aleatorios")
    for c in np.linspace(0.05, 0.999, 200):
        d = np.full(3, c)
        mu = 3 * (np.sqrt(c) - c)
        assert bound_full(d) <= mu + 1e-12 and bound_sigma(d) <= bound_full(d) + 1e-12
    for c in (0.9, 0.99, 0.999):
        d = np.full(3, c)
        mu = 3 * (np.sqrt(c) - c)
        rho = 2 * np.sqrt(c) / (1 + np.sqrt(c))
        s, D = sigma_delta(d)
        print(f"  c={c}: mu={mu:.6f}  cota={bound_full(d):.6f}  sigma(1-sigma)={bound_sigma(d):.6f}"
              f"  rho={rho:.4f}  1-sigma={1 - s:.4f}")


def check_B_random() -> None:
    print("(B) cota inferior vs unitarias aleatorias U = exp(eps A)")
    total = active = viol = 0
    for n in range(2, 13):
        for eps in (0.02, 0.05, 0.1, 0.2, 0.4, 0.8):
            for _ in range(300):
                A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
                A = (A - A.conj().T) / 2
                U = expm(eps * A)
                d = np.abs(np.diag(U)) ** 2
                lb = bound_full(d)
                M = pairing(U)
                total += 1
                active += lb > 0
                viol += M < lb - 1e-10
                assert bound_sigma(d) <= lb + 1e-12
    print(f"  {total} unitarias, cota activa (>0) en {active}, violaciones: {viol}")
    assert viol == 0


# ----------------------------------------------------------------------------
def oriented_graphs_passing(n: int) -> int:
    """Cuenta grafos orientados en n vertices donde todo vertice tiene grado de
    salida y de entrada >= 2 y todo arco i->j tiene un vecino de salida comun y
    un vecino de entrada comun (condicion necesaria para un cero interior)."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    P = len(pairs)
    total = 3 ** P
    POP = np.array([bin(v).count("1") for v in range(1 << n)], dtype=np.int8)
    bits = [1 << k for k in range(n)]
    count = 0
    chunk = 300000
    for start in range(0, total, chunk):
        end = min(start + chunk, total)
        idx = np.arange(start, end, dtype=np.int64)
        digits = np.empty((end - start, P), dtype=np.int8)
        t = idx.copy()
        for p in range(P):
            digits[:, p] = (t % 3).astype(np.int8)
            t //= 3
        out = np.zeros((end - start, n), dtype=np.int64)
        inn = np.zeros((end - start, n), dtype=np.int64)
        for p, (i, j) in enumerate(pairs):
            s = digits[:, p]
            a = s == 0
            b = s == 1
            out[:, i] |= a * bits[j]
            inn[:, j] |= a * bits[i]
            out[:, j] |= b * bits[i]
            inn[:, i] |= b * bits[j]
        ok = (POP[out].min(axis=1) >= 2) & (POP[inn].min(axis=1) >= 2)
        if not ok.any():
            continue
        d2, o2, i2 = digits[ok], out[ok], inn[ok]
        good = np.ones(d2.shape[0], dtype=bool)
        for p, (i, j) in enumerate(pairs):
            s = d2[:, p]
            for a, b, st in ((i, j, 0), (j, i, 1)):
                sel = s == st
                if not sel.any():
                    continue
                c1 = (o2[:, a] & o2[:, b]) != 0
                c2 = (i2[:, a] & i2[:, b]) != 0
                good &= ~sel | (c1 & c2)
        count += int(good.sum())
    return count


def check_C() -> None:
    print("(C) lema combinatorio: ningun grafo orientado con n<=6 cumple la condicion necesaria")
    for n in (4, 5, 6):
        k = oriented_graphs_passing(n)
        print(f"  n={n}: grafos orientados que cumplen = {k}  (de 3^{n * (n - 1) // 2})")
        assert k == 0
    T7 = paley_tournament(7).astype(int)
    ok = all(((T7[i] & T7[j]).any() and (T7[:, i] & T7[:, j]).any())
             for i in range(7) for j in range(7) if T7[i, j])
    print(f"  n=7 (Paley) cumple la condicion: {ok}")
    assert ok


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    check_A()
    check_B_identity_and_exact_n3()
    check_B_random()
    check_C()
    print("OK")
