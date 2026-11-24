"""
unistochastic.py — μ(d) exacta para n=3 vía la reducción unistocástica
======================================================================

La masa emparejada M(U) = Σ_{i<j} |U_ij·U_ji| = Σ_{i<j} √(B_ij·B_ji) depende SÓLO de
B = |U|∘² (moduli al cuadrado), no de las fases. Cuando U recorre las unitarias con
|U_ii|² = d_i fijo, B recorre exactamente las matrices UNISTOCÁSTICAS con diagonal d.

Para n=3, una doblemente estocástica con diagonal fija tiene UN grado de libertad x:

    B = [[d1,  x,           (1-d1)-x        ],
         [(1-d1-d2+d3)-x,  d2,   (d1-d3)+x  ],
         [(d2-d3)+x,       (1-d2)-x,   d3   ]]

(B_12=x=a, B_13=b, B_21=c, B_23=e, B_31=f, B_32=g; filas y columnas suman 1 por
construcción). Unistochasticidad de una 3×3 doblemente estocástica ⟺ los tres
términos de ortogonalidad de columnas (1,2) cierran un TRIÁNGULO:

    t1=√(B11·B12), t2=√(B21·B22), t3=√(B31·B32)  con  2·max(t) ≤ t1+t2+t3.

Por tanto (reducción EXACTA, no aproximada):

    μ(d) = mín { M(x) : x factible (B≥0) y B(x) unistocástica },   M(x)=√(ac)+√(bf)+√(eg).

Nota de normalización: M es la masa emparejada CRUDA (pares no ordenados, |·|), la de
`docs/cota_inferior_demostracion.md` — NO el W ponderado de FL-2.
"""

from __future__ import annotations

import numpy as np


def matriz_B(d, x):
    """La doblemente estocástica 3×3 con diagonal d y parámetro x."""
    d1, d2, d3 = d
    a = x
    b = (1.0 - d1) - x
    c = (1.0 - d1 - d2 + d3) - x
    e = (d1 - d3) + x
    fentry = (d2 - d3) + x
    g = (1.0 - d2) - x
    return np.array([[d1, a, b], [c, d2, e], [fentry, g, d3]], dtype=float)


def intervalo_factible(d):
    """[x_lo, x_hi] donde todas las entradas de B son ≥ 0."""
    d1, d2, d3 = d
    x_lo = max(0.0, d3 - d1, d3 - d2)
    x_hi = min(1.0 - d1, 1.0 - d1 - d2 + d3, 1.0 - d2)
    return x_lo, x_hi


def masa_emparejada(d, x):
    """M(x) = √(B12 B21) + √(B13 B31) + √(B23 B32) (masa emparejada cruda)."""
    B = matriz_B(d, x)
    prod = np.array([B[0, 1] * B[1, 0], B[0, 2] * B[2, 0], B[1, 2] * B[2, 1]])
    return float(np.sum(np.sqrt(np.clip(prod, 0.0, None))))


def es_unistochastica(d, x, tol=1e-12):
    """Triángulo de ortogonalidad de columnas (1,2): condición unistocástica n=3."""
    B = matriz_B(d, x)
    if np.any(B < -tol):
        return False
    t = np.sqrt(np.clip(np.array(
        [B[0, 0] * B[0, 1], B[1, 0] * B[1, 1], B[2, 0] * B[2, 1]]), 0.0, None))
    return bool(2.0 * t.max() <= t.sum() + tol)


def realizar(d, x, tol=1e-8):
    """Construye una U∈U(3) explícita con |U_ii|²=d_i y B=matriz_B(d,x), o None si
    B no es unistocástica. Verificación: fila 1 real √B_1j; fila 2 con fases que
    cierran el triángulo de ortogonalidad; fila 3 = complemento ortonormal (sus
    moduli salen √B_3j GRATIS por estocasticidad de columnas). Devuelve None si el
    triángulo no cierra o la verificación falla."""
    B = matriz_B(d, x)
    if np.any(B < -1e-12):
        return None
    B = np.clip(B, 0.0, None)
    r1 = np.sqrt(B[0])
    m = np.sqrt(B[0] * B[1])                        # m_j = √(B_1j B_2j)
    m0, m1, m2 = m
    if m0 < 1e-14 or m1 < 1e-14:
        return None                                 # borde: usar es_unistochastica
    cpsi = (m2 * m2 - m0 * m0 - m1 * m1) / (2.0 * m0 * m1)
    if abs(cpsi) > 1.0 + 1e-9:
        return None                                 # triángulo no cierra ⇒ no unist.
    a2 = float(np.arccos(np.clip(cpsi, -1.0, 1.0)))
    # r2_j = √B_2j · e^{iα_j}, con Σ_j m_j e^{iα_j}=0, α_1=0, α_2=a2:
    z = np.array([m0, m1 * np.exp(1j * a2), 0.0], dtype=complex)
    z[2] = -z[0] - z[1]
    fase = np.array([1.0 + 0j, z[1] / m1, z[2] / m2 if m2 > 1e-14 else 1.0])
    r2 = np.sqrt(B[1]) * fase
    # fila 3 = complemento ortonormal (Hermitiano) de r1, r2:
    Mrows = np.vstack([r1.astype(complex), r2])
    # base del núcleo por la izquierda: vector ortogonal a conj(r1), conj(r2)
    _, _, Vh = np.linalg.svd(Mrows.conj())
    r3 = Vh[-1].conj()
    U = np.vstack([r1.astype(complex), r2, r3])
    if np.max(np.abs(U.conj().T @ U - np.eye(3))) > tol:
        return None
    if np.max(np.abs(np.abs(np.diag(U)) ** 2 - np.array(d))) > tol:
        return None
    return U


def barrer(d, n=400_001, refinar=True):
    """μ_ds (mín sobre el intervalo factible) y μ_uni (mín sobre factible ∩ unist.),
    con el argmin de cada uno. Barrido fino + refinamiento local opcional."""
    x_lo, x_hi = intervalo_factible(d)
    if x_hi < x_lo - 1e-12:
        return None                                 # diagonal no factible
    x_hi = max(x_hi, x_lo)      # intervalo de un solo punto (p. ej. d=(1,0,0)): B determinada
    xs = np.linspace(x_lo, x_hi, n)
    B12 = xs
    d1, d2, d3 = d
    B21 = (1 - d1 - d2 + d3) - xs
    B13 = (1 - d1) - xs
    B31 = (d2 - d3) + xs
    B23 = (d1 - d3) + xs
    B32 = (1 - d2) - xs
    M = (np.sqrt(np.clip(B12 * B21, 0, None)) + np.sqrt(np.clip(B13 * B31, 0, None))
         + np.sqrt(np.clip(B23 * B32, 0, None)))
    # máscara unistocástica vectorizada (triángulo columnas 1,2)
    t1 = np.sqrt(np.clip(d1 * B12, 0, None))
    t2 = np.sqrt(np.clip(B21 * d2, 0, None))
    t3 = np.sqrt(np.clip(B31 * B32, 0, None))
    tmax = np.maximum(np.maximum(t1, t2), t3)
    uni = 2.0 * tmax <= (t1 + t2 + t3) + 1e-12
    i_ds = int(np.argmin(M))
    if not uni.any():
        return {"x_lo": x_lo, "x_hi": x_hi, "mu_ds": float(M[i_ds]),
                "x_ds": float(xs[i_ds]), "mu_uni": None, "x_uni": None,
                "frac_uni": 0.0}
    Mu = np.where(uni, M, np.inf)
    i_uni = int(np.argmin(Mu))
    res = {"x_lo": x_lo, "x_hi": x_hi, "mu_ds": float(M[i_ds]),
           "x_ds": float(xs[i_ds]), "mu_uni": float(Mu[i_uni]),
           "x_uni": float(xs[i_uni]), "frac_uni": float(uni.mean())}
    if refinar and np.isfinite(Mu[i_uni]):
        lo = xs[max(i_uni - 2, 0)]
        hi = xs[min(i_uni + 2, n - 1)]
        xr = np.linspace(lo, hi, 4001)
        for x in xr:
            if es_unistochastica(d, x):
                m = masa_emparejada(d, x)
                if m < res["mu_uni"]:
                    res["mu_uni"], res["x_uni"] = float(m), float(x)
    return res
