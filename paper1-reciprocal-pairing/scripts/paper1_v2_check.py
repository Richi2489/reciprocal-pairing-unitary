"""Paper 1 (version 2): numerical and symbolic companion.

The proofs in the paper are complete without this script, except where the
paper says a polynomial identity is "verified by expansion" (checked here
symbolically) and for the count of doubly regular tournaments of order 7.

Sections:
  1. Lemma "one invariant" and Theorem "closed form" (n=3).
  2. Identity (factor) of Proposition "unistochastic matrices with a given
     diagonal" and the two factorizations used in Theorem "the maximum of M
     on the fibre" (symbolic).
  3. That theorem: both cases of the maximum; Corollary "symmetric unitary
     matrices".
  4. n=4: Jacobi identity, Proposition "lower bound", U_theta, circulant
     family, value at (1/2,1/3,1/5,1/7).
  5. Doubly regular tournaments: Paley orders 7, 11, 19, 23; count of
     labelled doubly regular tournaments of order 7.
  (6. removed: superseded by section 10.)
  7. Lemma "splitting a diagonal entry" and the zeros it gives for n = 8, 9.
  8. The orthogonal matrix of order 11 (exact integer arithmetic).
  9. Order 11: the matrix as a row/column permutation of T - 2I; counts of
     pattern automorphisms for p = 11 and p = 7.
 10. Theorem "merging" and Corollary "decay": best partitions for small n,
     an explicit matrix for n = 6, the two decay bounds.

Run:  python scripts/paper1_v2_check.py      (exit code 0 = all pass)
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
import sympy as sp

rng = np.random.default_rng(2026)
OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"[{'ok' if cond else 'FAIL'}] {name} {detail}")


def haar(n):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def M(u):
    n = u.shape[0]
    return sum(abs(u[i, j] * u[j, i]) for i in range(n) for j in range(i + 1, n))


def rotation(al, be, ga):
    """Euler-Rodrigues rotation with diagonal (al, be, ga) (signed), or None."""
    q2 = np.array([1 + al + be + ga, 1 + al - be - ga, 1 - al + be - ga, 1 - al - be + ga]) / 4
    if q2.min() < -1e-13:
        return None
    q0, q1, q2_, q3 = np.sqrt(np.clip(q2, 0, None))
    return np.array([
        [1 - 2 * (q2_**2 + q3**2), 2 * (q1 * q2_ - q0 * q3), 2 * (q1 * q3 + q0 * q2_)],
        [2 * (q1 * q2_ + q0 * q3), 1 - 2 * (q1**2 + q3**2), 2 * (q2_ * q3 - q0 * q1)],
        [2 * (q1 * q3 - q0 * q2_), 2 * (q2_ * q3 + q0 * q1), 1 - 2 * (q1**2 + q2_**2)]])


def feasible3(a, b, c):
    return a + b <= 1 + c + 1e-15 and a + c <= 1 + b + 1e-15 and b + c <= 1 + a + 1e-15


# ---------------------------------------------------------------- 1
w = 0.0
for _ in range(20000):
    U = haar(3)
    a = np.abs(np.diag(U)); det = np.linalg.det(U)
    tau = (U[0, 0] * U[1, 1] * U[2, 2] * np.conj(det)).real
    for i, j, k in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
        w = max(w, abs(U[i, i] * U[j, j] - U[i, j] * U[j, i] - det * np.conj(U[k, k])),
                abs(abs(U[i, j] * U[j, i]) ** 2 - (a[i] ** 2 * a[j] ** 2 + a[k] ** 2 - 2 * tau)))
    assert abs(tau) <= a.prod() + 1e-12
check("1. Lemma (one invariant): cofactor identity and pair formula, 20000 Haar U(3)", w < 1e-12, f"max defect {w:.1e}")

w = 0.0; cnt = 0
while cnt < 5000:
    a, b, c = rng.uniform(0, 1, 3)
    if not feasible3(a, b, c):
        continue
    cnt += 1
    R = rotation(a, b, c)
    closed = abs(a * b - c) + abs(a * c - b) + abs(b * c - a)
    w = max(w, np.abs(R @ R.T - np.eye(3)).max(), abs(np.linalg.det(R) - 1),
            np.abs(np.diag(R) - [a, b, c]).max(), abs(M(R) - closed))
check("1. Theorem 1: rotation attains the closed form, 5000 feasible diagonals", w < 1e-12, f"max defect {w:.1e}")

# ---------------------------------------------------------------- 2
a, b, c, x = sp.symbols("a b c x", positive=True)
d1, d2, d3 = a**2, b**2, c**2
S = 1 - d1 - d2 + d3
B = sp.Matrix([[d1, x, 1 - d1 - x], [S - x, d2, d1 - d3 + x], [d2 - d3 + x, 1 - d2 - x, d3]])
L2 = [B[0, k] * B[1, k] for k in range(3)]
G = 2 * (L2[0] * L2[1] + L2[0] * L2[2] + L2[1] * L2[2]) - (L2[0]**2 + L2[1]**2 + L2[2]**2)
Pr = x * (S - x)
check("2. identity (factor): G(B(x)) = ((ab+c)^2 - x(S-x)) (x(S-x) - (ab-c)^2)",
      sp.expand(G - ((a * b + c)**2 - Pr) * (Pr - (a * b - c)**2)) == 0)
check("2. B(x) is doubly stochastic in form",
      all(sp.expand(sum(B.row(i)) - 1) == 0 for i in range(3)) and all(sp.expand(sum(B.col(i)) - 1) == 0 for i in range(3)))
f1 = sp.Rational(1, 4) * (1 + b + c - a) * (1 + a + c - b) * (1 + a + b - c) * (a + b + c - 1)
f2 = sp.Rational(1, 4) * (1 + a - b - c) * (1 + b - a - c) * (1 + c - a - b) * (1 + a + b + c)
check("2. factorizations of (ab+c)^2 - S^2/4 and S^2/4 - (ab-c)^2",
      sp.expand((a * b + c)**2 - S**2 / 4 - f1) == 0 and sp.expand(S**2 / 4 - (a * b - c)**2 - f2) == 0)

# ---------------------------------------------------------------- 3
w = 0.0; cnt = 0
while cnt < 3000:
    a_, b_, c_ = rng.uniform(0, 1, 3)
    if not feasible3(a_, b_, c_) or a_ + b_ + c_ > 1:
        continue
    cnt += 1
    R = rotation(-a_, -b_, -c_)
    target = a_ * b_ + a_ * c_ + b_ * c_ + a_ + b_ + c_
    budget = (3 - a_**2 - b_**2 - c_**2) / 2
    s_ = a_ + b_ + c_
    w = max(w, np.abs(R @ R.T - np.eye(3)).max(), np.abs(np.abs(np.diag(R)) - [a_, b_, c_]).max(),
            abs(M(R) - target), abs(budget - target - (1 - s_) * (3 + s_) / 2))
check("3. Theorem 2, case a+b+c <= 1: rotation with diagonal (-a,-b,-c), 3000 diagonals", w < 1e-12, f"max defect {w:.1e}")

# case a+b+c >= 1: the symmetric doubly stochastic matrix is unistochastic; build U explicitly
w = 0.0; cnt = 0
while cnt < 3000:
    a_, b_, c_ = rng.uniform(0, 1, 3)
    if not feasible3(a_, b_, c_) or a_ + b_ + c_ < 1:
        continue
    cnt += 1
    d = np.array([a_**2, b_**2, c_**2])
    Bs = np.diag(d).astype(float)
    for i, j, k in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
        Bs[i, j] = Bs[j, i] = (1 - d[i] - d[j] + d[k]) / 2
    assert Bs.min() > -1e-13
    r = np.sqrt(np.clip(Bs, 0, None))
    L = r[0] * r[1]                      # row pair (1,2): triangle with sides L
    Lc = r[0] * r[2]                     # row pair (1,3)
    tri = lambda v: v[0] <= v[1] + v[2] + 1e-9 and v[1] <= v[0] + v[2] + 1e-9 and v[2] <= v[0] + v[1] + 1e-9
    assert tri(L) and tri(Lc)
    # phases: row 1 real; rows 2, 3 from the two triangles (law of cosines)
    def angles(v):
        # phases t1, t2 with v0 + v1 e^{i t1} + v2 e^{i t2} = 0
        if v[1] < 1e-14 or v[2] < 1e-14:
            return (np.pi, np.pi)
        cos1 = (v[2]**2 - v[0]**2 - v[1]**2) / (2 * v[0] * v[1]) if v[0] > 1e-14 else 0.0
        t1 = np.arccos(np.clip(cos1, -1, 1))
        z = -(v[0] + v[1] * np.exp(1j * t1))
        return (t1, np.angle(z))
    best = np.inf
    for s2 in (1, -1):
        for s3 in (1, -1):
            t = angles(L); u = angles(Lc)
            U = r.astype(complex)
            U[1, 1] *= np.exp(-1j * s2 * t[0]); U[1, 2] *= np.exp(-1j * s2 * t[1])
            U[2, 1] *= np.exp(-1j * s3 * u[0]); U[2, 2] *= np.exp(-1j * s3 * u[1])
            best = min(best, np.abs(U @ U.conj().T - np.eye(3)).max())
    w = max(w, best)
check("3. Theorem 2, case a+b+c >= 1: a unitary with symmetric moduli exists, 3000 diagonals", w < 1e-7, f"max unitarity defect {w:.1e}")

mn = np.inf
for _ in range(20000):
    V = haar(3); Sm = V @ V.T
    mn = min(mn, np.abs(np.diag(Sm)).sum())
K = (2 / 3) * np.ones((3, 3)) - np.eye(3)
check("3. Corollary: |S11|+|S22|+|S33| >= 1 on 20000 symmetric unitaries; equality at (2/3)J - I",
      mn >= 1 - 1e-9 and np.allclose(K @ K.T, np.eye(3)) and abs(np.abs(np.diag(K)).sum() - 1) < 1e-15, f"min {mn:.6f}")

# ---------------------------------------------------------------- 4
PARTS = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]
lam = lambda a: sum(abs(a[i] * a[j] - a[k] * a[l]) for i, j, k, l in PARTS)
e1, sl = 0.0, np.inf
for _ in range(20000):
    U = haar(4); det = np.linalg.det(U); a = np.abs(np.diag(U))
    for i, j, k, l in PARTS:
        e1 = max(e1, abs(U[i, i] * U[j, j] - U[i, j] * U[j, i] - det * np.conj(U[k, k] * U[l, l] - U[k, l] * U[l, k])))
        sl = min(sl, abs(U[i, j] * U[j, i]) + abs(U[k, l] * U[l, k]) - abs(a[i] * a[j] - a[k] * a[l]))
check("4. Jacobi identity and Proposition (lower bound), 20000 Haar U(4)", e1 < 1e-12 and sl > -1e-12, f"identity {e1:.1e}, min slack {sl:.1e}")
v = lam(np.sqrt([1 / 2, 1 / 3, 1 / 5, 1 / 7]))
check("4. lambda(1/2,1/3,1/5,1/7) = 0.3462...", 0.3462 < v < 0.3463, f"{v:.10f}")
w = 0.0
for th in np.linspace(0, np.pi / 2, 181):
    cth, sth = np.cos(th), np.sin(th)
    U = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [cth, 0, 0, sth], [-sth, 0, 0, cth]], float)
    w = max(w, np.abs(U @ U.T - np.eye(4)).max(), M(U), np.abs(np.diag(U) - [0, 0, 0, cth]).max())
check("4. U_theta: orthogonal, diagonal (0,0,0,cos), M = 0", w < 1e-14)
nu = lambda t: 2 * t if t <= 0.25 else 2 * (1 - np.sqrt(t)) * (3 * np.sqrt(t) - 1)
w = 0.0
for t in np.linspace(0, 1, 201):
    s_ = np.sqrt(t)
    cph, x2 = (0.0, 2 * t) if t <= 0.25 else (2 * s_ - 1, s_)
    e = np.exp(1j * np.arccos(cph)); xx, yy = np.sqrt(x2), np.sqrt(1 - x2)
    u = [xx / 2 * (1 + e), 1j * yy / 2 * (1 - 1j * e), xx / 2 * (1 - e), 1j * yy / 2 * (1 + 1j * e)]
    U = np.array([[u[(k - j) % 4] for k in range(4)] for j in range(4)])
    w = max(w, np.abs(U @ U.conj().T - np.eye(4)).max(), np.abs(np.abs(np.diag(U))**2 - t).max(), abs(M(U) - nu(t)))
check("4. circulant family: unitary, diagonal t, M = nu(t), 201 values", w < 1e-13 and abs(nu(4 / 9) - 2 / 3) < 1e-14)

# ---------------------------------------------------------------- 5
for q in (7, 11, 19, 23):
    k = (q - 3) // 4
    qr = {(z * z) % q for z in range(1, q)}
    T = np.array([[1.0 if ((j - i) % q) in qr else 0.0 for j in range(q)] for i in range(q)])
    U = (T - k * np.eye(q)) / (k + 1)
    check(f"5. Paley q={q}: (T-kI)/(k+1) orthogonal, M=0, diagonal (k/(k+1))^2",
          np.allclose(U @ U.T, np.eye(q)) and M(U) == 0 and abs(U[0, 0]**2 - (k / (k + 1))**2) < 1e-15)

# labelled doubly regular tournaments of order 7: regular, every pair with exactly one common out-neighbour
pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
count = 0
def extend(idx, out):
    """out[i] = bitmask of out-neighbours; prune on out-degree."""
    global count
    if idx == len(pairs):
        if all(bin(o).count("1") == 3 for o in out) and \
           all(bin(out[i] & out[j]).count("1") == 1 for i, j in pairs):
            count += 1
        return
    i, j = pairs[idx]
    for (u, v) in ((i, j), (j, i)):
        if bin(out[u]).count("1") < 3:
            out[u] |= 1 << v
            extend(idx + 1, out)
            out[u] &= ~(1 << v)
extend(0, [0] * 7)
check("5. labelled doubly regular tournaments of order 7 = 240 = 7!/21", count == 240, f"count {count}")

# ---------------------------------------------------------------- 7
def split(U, v, s_):
    n = U.shape[0]; dv = abs(U[v, v]) ** 2; t_ = dv / s_
    T0 = np.array([[np.sqrt(s_), np.sqrt((1 - s_) * (1 - t_))], [0, np.sqrt(t_)]], complex)
    A, S1, Bh = np.linalg.svd(T0); A2, S2, B2h = np.linalg.svd(np.diag([U[v, v], 1.0]))
    assert np.allclose(S1, S2)
    L0 = A @ A2.conj().T; R0 = B2h.conj().T @ Bh
    W = np.zeros((n + 1, n + 1), complex); W[:n, :n] = U; W[n, n] = 1
    L = np.eye(n + 1, dtype=complex); R = np.eye(n + 1, dtype=complex)
    idx = [v, n]
    for p_ in range(2):
        for q_ in range(2):
            L[idx[p_], idx[q_]] = L0[p_, q_]; R[idx[p_], idx[q_]] = R0[p_, q_]
    return L @ W @ R, t_

w = 0.0; bad = 0
for _ in range(3000):
    n = int(rng.integers(3, 7)); U = haar(n); v = int(rng.integers(n)); dv = abs(U[v, v]) ** 2
    s_ = rng.uniform(dv, 1); Up, t_ = split(U, v, s_)
    dn = np.abs(np.diag(Up)) ** 2
    rest = [i for i in range(n) if i != v]
    w = max(w, np.abs(Up @ Up.conj().T - np.eye(n + 1)).max(), abs(dn[v] - s_), abs(dn[n] - t_),
            np.abs(dn[rest] - np.abs(np.diag(U))[rest] ** 2).max())
    bad += M(Up) > M(U) + 1e-10
check("7. splitting lemma on 3000 Haar unitaries: U' unitary, diagonal as stated, M(U') <= M(U)", w < 1e-12 and bad == 0, f"max defect {w:.1e}")
qr7 = {(z * z) % 7 for z in range(1, 7)}
T7 = np.array([[1.0 if ((j - i) % 7) in qr7 else 0.0 for j in range(7)] for i in range(7)])
U8, _ = split(((T7 - np.eye(7)) / 2).astype(complex), 0, 0.5)
U9, _ = split(U8, 3, 0.6)
d8 = np.abs(np.diag(U8)) ** 2
check("7. zeros with all d_i in (0,1) for n = 8 (diagonal 1/2, 1/4 x6, 1/2) and n = 9",
      M(U8) < 1e-14 and M(U9) < 1e-14 and np.allclose(d8, [0.5] + [0.25] * 6 + [0.5])
      and (np.abs(np.diag(U9)) ** 2).min() > 0 and (np.abs(np.diag(U9)) ** 2).max() < 1)

# ---------------------------------------------------------------- 8
E3 = np.array([[int(z) for z in r.split()] for r in """1 1 0 1 1 2 0 0 0 1 0
0 1 1 0 1 -1 -1 0 0 0 -2
1 0 1 1 0 -1 2 1 0 0 0
0 1 0 1 -2 0 -1 1 1 0 0
0 0 -2 0 1 -1 0 1 1 1 0
0 0 0 1 0 -1 -1 0 -2 1 1
1 0 0 0 1 0 -1 1 0 -2 1
1 1 0 0 0 -1 0 -2 1 0 1
1 -2 1 0 0 0 -1 0 1 1 0
0 1 1 -2 0 0 0 1 0 1 1
2 0 -1 -1 -1 0 0 0 -1 0 -1""".splitlines()])
qr11 = {(z * z) % 11 for z in range(1, 11)}
pat = np.array([[1 if (i == j or ((j - i) % 11) in qr11) else 0 for j in range(11)] for i in range(11)])
dg = sorted(abs(int(z)) for z in np.diag(E3))
check("8. order 11: (3E)(3E)^T = 9I exactly, support = T + I, diagonal moduli 1/3 (x10) and 2/3 (x1)",
      (E3 @ E3.T == 9 * np.eye(11, dtype=int)).all() and ((E3 != 0).astype(int) == pat).all()
      and dg == [1] * 10 + [2])

# ---------------------------------------------------------------- 9
# Remark after the order-11 proposition: 3E comes from T - 2I by two different
# permutations and sign changes; automorphism counts of the pattern.
T11 = pat - np.eye(11, dtype=int)
B11 = T11 - 2 * np.eye(11, dtype=int)
r_ = (0, 2, 6, 5, 7, 8, 1, 10, 3, 9, 4)
c_ = (4, 3, 7, 9, 5, 0, 6, 10, 8, 1, 2)
C11 = B11[np.ix_(r_, c_)].copy()
okC = ((C11 != 0).astype(int) == pat).all()
C11[10, :] *= -1; C11[:, 5] *= -1; C11[:, 6] *= -1
check("9. order 11: 3E = (T - 2I) with rows r, columns c permuted, row 10 and columns 5, 6 negated",
      okC and (C11 == E3).all())


def pattern_automorphisms(Pm):
    """Pairs (row permutation, column permutation) preserving the 0/1 pattern Pm.
    Returns (total, simultaneous, sorted list of numbers of diagonal positions kept on the diagonal)."""
    n = Pm.shape[0]
    rows = [frozenset(np.nonzero(Pm[i])[0].tolist()) for i in range(n)]
    rowset = {s_: i for i, s_ in enumerate(rows)}
    res = [0, 0, set()]

    def bt(pi, used):
        k = len(pi)
        if k == n:
            img = [rowset[frozenset(pi[j] for j in s_)] for s_ in rows]
            res[0] += 1
            res[1] += tuple(img) == tuple(pi)
            res[2].add(sum(1 for i in range(n) if img[i] == pi[i]))
            return
        for x in range(n):
            if x in used:
                continue
            pi.append(x); used.add(x)
            if all(frozenset(pi[j] for j in s_) in rowset for s_ in rows if max(s_) == k):
                bt(pi, used)
            pi.pop(); used.discard(x)
    bt([], set())
    return res[0], res[1], sorted(res[2])


a11 = pattern_automorphisms(pat)
a7 = pattern_automorphisms((T7 + np.eye(7)).astype(int))
check("9. pattern T + I: 660 pairs (55 simultaneous; the other 605 keep exactly one diagonal position) for p = 11; 168 (21) for p = 7",
      a11 == (660, 55, [1, 11]) and a7[:2] == (168, 21), f"{a11} {a7}")

# ---------------------------------------------------------------- 10
# Merging theorem and the decay corollary.
import math


def first_primes(n):
    out, x = [], 2
    while len(out) < n:
        if all(x % q for q in out if q * q <= x):
            out.append(x)
        x += 1
    return out


mu3abc = lambda a, b, c: abs(a * b - c) + abs(a * c - b) + abs(b * c - a)


def best_partition(n):
    ps = first_primes(n); best = None
    for lab in itertools.product(range(3), repeat=n):
        if lab[0] != 0 or len(set(lab)) < 3:
            continue
        pr = [1, 1, 1]
        for p_, l_ in zip(ps, lab):
            pr[l_] *= p_
        A = [x ** -0.5 for x in pr]
        if feasible3(*A) and (best is None or mu3abc(*A) < best[0]):
            best = (mu3abc(*A), sorted(pr), lab)
    return best


quoted = {4: (0.7276, [5, 6, 7]), 5: (0.5996, [11, 14, 15]), 6: (0.4426, [26, 33, 35]),
          9: (0.1169, [595, 598, 627]), 12: (0.0214, [19285, 19499, 19734])}
okm = True
for n_, (val, prods) in quoted.items():
    v, pr, lab = best_partition(n_)
    okm &= pr == prods and v <= val < v + 1e-4
    print(f"     n={n_}: best partition products {pr}, bound {v:.6f}")
check("10. merging bound, best partitions for n = 4, 5, 6, 9, 12 as quoted (rounded up)", okm)

# explicit matrix for n = 6: rotation with diagonal moduli 26^-1/2, 33^-1/2, 35^-1/2, then three splittings
groups = [(2, 13), (3, 11), (5, 7)]
A = [1 / math.sqrt(p_ * q_) for p_, q_ in groups]
U6 = rotation(*A).astype(complex)
for v_, (p_, q_) in enumerate(groups):
    U6, t_ = split(U6, v_, 1 / p_)
d6 = np.abs(np.diag(U6)) ** 2
target = np.array([1 / 2, 1 / 3, 1 / 5, 1 / 13, 1 / 11, 1 / 7])
reach = (np.abs(U6) > 1e-9).astype(int)
conn = np.linalg.matrix_power(reach + np.eye(6, dtype=int), 6)
check("10. n = 6: explicit unitary with diagonal {1/2, 1/3, 1/5, 1/7, 1/11, 1/13}, M <= 0.4426 < 1.403, irreducible support",
      np.abs(U6 @ U6.conj().T - np.eye(6)).max() < 1e-12 and np.abs(d6 - target).max() < 1e-12
      and M(U6) <= 0.4426 and abs(M(U6) - 0.2531) < 5e-5 and int((reach == 0).sum()) == 3 and (conn > 0).all(),
      f"M = {M(U6):.6f}, zero entries: {int((reach == 0).sum())}")
ps6 = first_primes(6); ds = []
for tri in itertools.combinations(range(6), 3):
    if 0 not in tri:
        continue
    oth = [i for i in range(6) if i not in tri]
    A1 = [ps6[i] ** -0.5 for i in tri]; A2 = [ps6[i] ** -0.5 for i in oth]
    if feasible3(*A1) and feasible3(*A2):
        ds.append(mu3abc(*A1) + mu3abc(*A2))
check("10. n = 6: best direct sum of two three-dimensional minimizers gives 1.4034", abs(min(ds) - 1.4034) < 5e-5, f"{min(ds):.6f} ({len(ds)} feasible splittings)")

okd = True
for n_ in (6, 7, 8, 9, 12, 30, 300):
    m_ = n_ // 3; ps = first_primes(n_)
    grp = [ps[k::3] for k in range(3)]
    A = [math.prod(g) ** -0.5 for g in grp]
    bound = 6 / math.sqrt(math.factorial(m_ + 1))
    okd &= all(len(g) >= m_ for g in grp) and max(A) <= 1 / math.sqrt(math.factorial(m_ + 1)) <= 0.5 \
        and feasible3(*A) and mu3abc(*A) <= 2 * sum(A) <= bound
check("10. decay corollary: mu(1/p_1..1/p_n) <= 6/sqrt((floor(n/3)+1)!) for n = 6..300 (residue-class partition)", okd)
okq = True
for _ in range(2000):
    n_ = int(rng.integers(3, 40)); q = rng.uniform(0.01, 0.99); m_ = n_ // 3
    if q ** m_ > 0.25:
        continue
    d = rng.uniform(0, q, n_); grp = [d[k::3] for k in range(3)]
    A = [math.sqrt(np.prod(g)) for g in grp]
    okq &= feasible3(*A) and mu3abc(*A) <= 6 * q ** (m_ / 2) + 1e-15
check("10. decay corollary: mu(d) <= 6 q^(m/2) when all d_i <= q and q^m <= 1/4 (2000 random d)", okq)

print("ALL PASS" if OK else "SOME CHECK FAILED")
sys.exit(0 if OK else 1)
