"""Paper 1: numerical confirmation of the n=3 unistochastic feasible-set geometry.

Verification script for the revision answering the external review (item P1-A).
For the one-parameter family of Section 3.2 it computes

    I = {x : all off-diagonal entries of B(x) >= 0}     (doubly stochastic segment)
    K = {x in I : B(x) is unistochastic}                (chain-link / triangle)

using the exact chain-link (triangle) condition on the three links
L_k = sqrt(B_ik B_jk), k = 1,2,3, of a pair of rows (i,j), written in Heron form
so that it is a polynomial in x:

    G_ij(x) = 2(a1 a2 + a1 a3 + a2 a3) - (a1^2 + a2^2 + a3^2) >= 0,
    a_k = B_ik B_jk = L_k^2.

Since the links are non-negative and the sum of any two of (-L1+L2+L3),
(L1-L2+L3), (L1+L2-L3) equals 2 L_k >= 0, at most one of the three can be
negative; hence the Heron product is >= 0 *exactly* when all three triangle
inequalities hold.  G_ij is a quartic in x, so K is obtained from its real roots
(no grid search needed for the boundaries).

Reported: the connected components of K, min/max of W over K, whether the
symmetry centre S/2 lies in K, agreement of all three row pairs and all three
column pairs (so that one pair suffices at n=3), the constant-family merge at
c = 1/9, cross-checks against the reviewer's c = 0.01 numbers and the Table 1
minima, and the P1-D identity.

Run:  python scripts/paper1_feasible_set_check.py
"""

from __future__ import annotations

import itertools
import numpy as np

np.set_printoptions(precision=12, suppress=False)


# ----------------------------------------------------------------------
# the one-parameter family of Section 3.2, as polynomials in x
# ----------------------------------------------------------------------
def entries(d, x):
    """Return the 3x3 array of B entries evaluated on the array x."""
    d1, d2, d3 = d
    x = np.asarray(x, dtype=float)
    one = np.ones_like(x)
    return {
        (0, 0): d1 * one, (0, 1): x, (0, 2): (1.0 - d1) - x,
        (1, 0): (1.0 - d1 - d2 + d3) - x, (1, 1): d2 * one, (1, 2): (d1 - d3) + x,
        (2, 0): (d2 - d3) + x, (2, 1): (1.0 - d2) - x, (2, 2): d3 * one,
    }


def W_of_x(d, x):
    E = entries(d, x)
    tot = np.zeros_like(np.asarray(x, dtype=float))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        tot = tot + np.sqrt(np.maximum(E[(i, j)] * E[(j, i)], 0.0))
    return tot


def ds_interval(d):
    """Interval I on which all six off-diagonal entries are >= 0."""
    d1, d2, d3 = d
    lo = max(0.0, d3 - d1, d3 - d2)
    hi = min(1.0 - d1, 1.0 - d1 - d2 + d3, 1.0 - d2)
    return lo, hi


def G(d, x, i, j, axis="row"):
    """Heron/chain-link polynomial for the row (or column) pair (i,j)."""
    E = entries(d, x)
    if axis == "row":
        a = [E[(i, k)] * E[(j, k)] for k in range(3)]
    else:
        a = [E[(k, i)] * E[(k, j)] for k in range(3)]
    a1, a2, a3 = a
    return 2.0 * (a1 * a2 + a1 * a3 + a2 * a3) - (a1 * a1 + a2 * a2 + a3 * a3)


def G_roots_in(d, lo, hi, i=0, j=1, axis="row"):
    """Real roots of the quartic G in (lo,hi), found from its coefficients."""
    # sample-and-fit is unnecessary: build coefficients by exact polynomial algebra
    d1, d2, d3 = d
    P = np.polynomial.polynomial  # coefficient arrays in increasing powers
    ent = {
        (0, 0): np.array([d1]), (0, 1): np.array([0.0, 1.0]),
        (0, 2): np.array([1.0 - d1, -1.0]),
        (1, 0): np.array([1.0 - d1 - d2 + d3, -1.0]), (1, 1): np.array([d2]),
        (1, 2): np.array([d1 - d3, 1.0]),
        (2, 0): np.array([d2 - d3, 1.0]), (2, 1): np.array([1.0 - d2, -1.0]),
        (2, 2): np.array([d3]),
    }
    if axis == "row":
        a = [P.polymul(ent[(i, k)], ent[(j, k)]) for k in range(3)]
    else:
        a = [P.polymul(ent[(k, i)], ent[(k, j)]) for k in range(3)]
    poly = np.zeros(1)
    for p, q in itertools.combinations(a, 2):
        poly = P.polyadd(poly, 2.0 * P.polymul(p, q))
    for p in a:
        poly = P.polysub(poly, P.polymul(p, p))
    # trim and take roots
    poly = P.polytrim(poly, tol=1e-14)
    if len(poly) < 2:
        return []
    r = P.polyroots(poly)
    out = sorted(float(z.real) for z in np.atleast_1d(r)
                 if abs(np.imag(z)) < 1e-9 and lo - 1e-12 < z.real < hi + 1e-12)
    return out


def components_of_K(d, i=0, j=1, axis="row"):
    """Connected components of K, from the roots of the quartic G."""
    lo, hi = ds_interval(d)
    if hi <= lo + 1e-15:
        return (lo, hi), []
    brk = [lo] + [r for r in G_roots_in(d, lo, hi, i, j, axis) if lo < r < hi] + [hi]
    brk = sorted(set(round(b, 14) for b in brk))
    comps = []
    for a, b in zip(brk[:-1], brk[1:]):
        if b - a < 1e-14:
            continue
        mid = 0.5 * (a + b)
        if float(G(d, np.array([mid]), i, j, axis)[0]) >= 0.0:
            if comps and abs(comps[-1][1] - a) < 1e-13:
                comps[-1] = (comps[-1][0], b)
            else:
                comps.append((a, b))
    return (lo, hi), comps


def maxW_on(d, a, b, n=200001):
    xs = np.linspace(a, b, n)
    w = W_of_x(d, xs)
    k = int(w.argmax())
    return float(w[k]), float(xs[k])


def report(d, label=""):
    d = tuple(float(t) for t in d)
    S = 1.0 - d[0] - d[1] + d[2]
    (lo, hi), comps = components_of_K(d)
    print("=" * 78)
    print(f"d = ({d[0]:.6g}, {d[1]:.6g}, {d[2]:.6g})   {label}")
    print(f"  S = B12+B21 = {S:.12f}    centre S/2 = {S/2:.12f}")
    if hi <= lo + 1e-15:
        print("  doubly stochastic segment I is EMPTY  ->  d is infeasible")
        return
    print(f"  I = [{lo:.12f}, {hi:.12f}]   (length {hi-lo:.12f})")
    if not comps:
        print("  K is EMPTY")
        return
    print(f"  K has {len(comps)} component(s):")
    for (a, b) in comps:
        print(f"      [{a:.12f}, {b:.12f}]"
              f"   W(a)={float(W_of_x(d,np.array([a]))[0]):.12f}"
              f"   W(b)={float(W_of_x(d,np.array([b]))[0]):.12f}")

    endpoint_vals = [float(W_of_x(d, np.array([t]))[0]) for c in comps for t in c]
    wmin = min(endpoint_vals)
    wmax, xmax = -np.inf, None
    for (a, b) in comps:
        m, xm = maxW_on(d, a, b)
        if m > wmax:
            wmax, xmax = m, xm
    # confirm the minimum really is at an endpoint (concavity on each component)
    wmin_grid, xmin_grid = np.inf, None
    for (a, b) in comps:
        xs = np.linspace(a, b, 200001)
        w = W_of_x(d, xs)
        k = int(w.argmin())
        if w[k] < wmin_grid:
            wmin_grid, xmin_grid = float(w[k]), float(xs[k])

    budget = (3.0 - sum(d)) / 2.0
    incentre = any(a - 1e-12 <= S / 2 <= b + 1e-12 for (a, b) in comps)
    xout = (comps[0][0], comps[-1][1])
    wout = [float(W_of_x(d, np.array([t]))[0]) for t in xout]
    print(f"  min W over K = {wmin:.12f}   (dense-grid min {wmin_grid:.12f}"
          f" at x={xmin_grid:.9f})")
    print(f"  max W over K = {wmax:.12f}   at x = {xmax:.12f}")
    print(f"  budget (3-sum d)/2 = {budget:.12f}    max/budget = {wmax/budget:.9f}")
    print(f"  S/2 in K ?  {incentre}    G(S/2) = {float(G(d,np.array([S/2]),0,1)[0]):+.6e}")
    print(f"  outermost points of K: {xout[0]:.12f}, {xout[1]:.12f}"
          f"   W = {wout[0]:.12f}, {wout[1]:.12f}   (diff {wout[0]-wout[1]:+.2e})")
    print(f"  symmetry of K: x_left + x_right - S = {xout[0]+xout[1]-S:+.3e}")

    # one pair suffices?  compare all 3 row pairs and all 3 column pairs
    xs = np.linspace(lo, hi, 200001)
    base = G(d, xs, 0, 1, "row") >= 0.0
    bad = []
    for axis in ("row", "col"):
        for (i, j) in itertools.combinations(range(3), 2):
            other = G(d, xs, i, j, axis) >= 0.0
            nd = int(np.sum(other != base))
            if nd > 2:
                bad.append((axis, i + 1, j + 1, nd))
    if not bad:
        print("  cross-check: the 3 row pairs and the 3 column pairs cut out"
              " the SAME set (one pair suffices at n=3)")
    else:
        for axis, i, j, nd in bad:
            print(f"  !! {axis} pair ({i},{j}) disagrees at {nd} of {len(xs)} points")


# ----------------------------------------------------------------------
# constant family d=(c,c,c): closed-form geometry and the merge at c = 1/9
# ----------------------------------------------------------------------
def constant_family_check():
    print()
    print("#" * 78)
    print("# constant family d=(c,c,c):  t = sqrt(x(1-c-x)),  W = 3t")
    print("#   triangle <=> w|u-v| <= t <= w(u+v) <=> sqrt(c)-c <= t <= sqrt(c)+c")
    print("#   plus t <= (1-c)/2 (AM-GM);  midpoint feasible iff c >= 1/9")
    print("#" * 78)
    hdr = (f"{'c':>13} {'#comp':>6} {'min W':>16} {'3(sqrt c-c)':>16} "
           f"{'max W':>16} {'3 min(sqc+c,(1-c)/2)':>22} {'S/2 in K':>9}")
    print(hdr)
    for c in [0.005, 0.01, 0.05, 1.0 / 9 - 1e-5, 1.0 / 9, 1.0 / 9 + 1e-5,
              0.15, 0.2, 0.4, 0.6]:
        d = (c, c, c)
        (lo, hi), comps = components_of_K(d)
        S = 1.0 - c
        wmin = min(float(W_of_x(d, np.array([t]))[0]) for cc in comps for t in cc)
        wmax = max(maxW_on(d, a, b, 60001)[0] for (a, b) in comps)
        pred_min = 3.0 * (np.sqrt(c) - c)
        pred_max = 3.0 * min(np.sqrt(c) + c, (1.0 - c) / 2.0)
        inc = any(a - 1e-12 <= S / 2 <= b + 1e-12 for (a, b) in comps)
        print(f"{c:13.9f} {len(comps):6d} {wmin:16.10f} {pred_min:16.10f} "
              f"{wmax:16.10f} {pred_max:22.10f} {str(inc):>9}")

    a, b = 0.01, 0.3
    for _ in range(80):
        m = 0.5 * (a + b)
        _, comps = components_of_K((m, m, m))
        if len(comps) >= 2:
            a = m
        else:
            b = m
    print()
    print(f"  measured merge threshold c* in [{a:.14f}, {b:.14f}]")
    print(f"  1/9 = {1.0/9:.14f}    |c* - 1/9| <= {max(abs(a-1/9),abs(b-1/9)):.3e}")
    print("  algebra: (1-c)/2 <= sqrt(c)+c  <=>  3c + 2 sqrt(c) - 1 >= 0"
          "  <=>  (3 sqrt(c) - 1)(sqrt(c) + 1) >= 0  <=>  c >= 1/9   [exact]")


# ----------------------------------------------------------------------
# reviewer / Table 1 cross-checks
# ----------------------------------------------------------------------
def reviewer_crosscheck():
    print()
    print("#" * 78)
    print("# cross-check: reviewer's c=0.01 numbers, and the Table 1 minima")
    print("#" * 78)
    d = (0.01, 0.01, 0.01)
    _, comps = components_of_K(d)
    print("  c=0.01 components of K:")
    for (a, b) in comps:
        print(f"    [{a:.14f}, {b:.14f}]")
    print("  reviewer:  [0.00825057781236, 0.01237695869343]"
          " u [0.97762304130657, 0.98174942218763]")
    wmin = min(float(W_of_x(d, np.array([t]))[0]) for c in comps for t in c)
    wmax = max(maxW_on(d, a, b, 60001)[0] for (a, b) in comps)
    print(f"  min W = {wmin:.12f}  (reviewer 0.27)"
          f"     max W = {wmax:.12f}  (reviewer 0.33)")
    print(f"  W at the (INFEASIBLE) midpoint x=0.495 = "
          f"{float(W_of_x(d,np.array([0.495]))[0]):.6f}"
          f"   = budget {(3-0.03)/2:.6f}")
    print()
    printed = [
        ((0.5, 1 / 3, 0.2), 0.7490, "von Mangoldt (1/2,1/3,1/5)"),
        ((0.4, 0.4, 0.4), 0.6974, "constant 0.4"),
        ((0.5, 0.45, 0.4), 0.6646, "slightly distinct"),
        ((0.6, 0.4, 0.4), 0.6597, "nearly degenerate"),
    ]
    print(f"  {'case':>28} {'computed min W':>18} {'printed':>9}")
    for dd, val, name in printed:
        _, comps = components_of_K(dd)
        wmin = min(float(W_of_x(dd, np.array([t]))[0]) for c in comps for t in c)
        print(f"  {name:>28} {wmin:18.14f} {val:9.4f}")


# ----------------------------------------------------------------------
# the P1-D identity
# ----------------------------------------------------------------------
def identity_check():
    print()
    print("#" * 78)
    print("# identity  W(B) = (n - sum d)/2 - (1/2) sum_{i<j}(sqrt B_ij - sqrt B_ji)^2")
    print("#" * 78)
    rng = np.random.default_rng(20260910)
    worst = 0.0
    for n in (3, 4, 5, 6, 8):
        for _ in range(400):
            A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            Q, R = np.linalg.qr(A)
            Q = Q * (np.diag(R) / np.abs(np.diag(R)))
            B = np.abs(Q) ** 2
            dg = np.diag(B)
            lhs = sum(np.sqrt(B[i, j] * B[j, i])
                      for i in range(n) for j in range(i + 1, n))
            rhs = (n - dg.sum()) / 2.0 - 0.5 * sum(
                (np.sqrt(B[i, j]) - np.sqrt(B[j, i])) ** 2
                for i in range(n) for j in range(i + 1, n))
            worst = max(worst, abs(lhs - rhs))
    print(f"  max |LHS - RHS| over 2000 random unistochastic B (n=3,4,5,6,8): "
          f"{worst:.3e}")
    # P1-D: the same deficit as a quarter of the FULL squared Frobenius norm
    worst2 = 0.0
    for n in (3, 4, 5, 6, 8):
        for _ in range(200):
            A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            Q, R = np.linalg.qr(A)
            Q = Q * (np.diag(R) / np.abs(np.diag(R)))
            B = np.abs(Q) ** 2
            dg = np.diag(B)
            W = sum(np.sqrt(B[i, j] * B[j, i])
                    for i in range(n) for j in range(i + 1, n))
            deficit = (n - dg.sum()) / 2.0 - W
            pairs = 0.5 * sum((np.sqrt(B[i, j]) - np.sqrt(B[j, i])) ** 2
                              for i in range(n) for j in range(i + 1, n))
            frob = 0.25 * np.linalg.norm(np.sqrt(B) - np.sqrt(B).T, "fro") ** 2
            worst2 = max(worst2, abs(deficit - pairs), abs(deficit - frob))
    print("  deficit = (1/2) sum_{i<j} (sqrt B_ij - sqrt B_ji)^2 "
          "= (1/4) ||sqrt B - sqrt B^T||_F^2   [FULL Frobenius: a QUARTER]")
    print(f"  max discrepancy over 1000 random unistochastic B: {worst2:.3e}")


# ----------------------------------------------------------------------
# P1-B: doubly stochastic feasibility is NOT unitary feasibility
# ----------------------------------------------------------------------
def ds_vs_unitary_check():
    """d = (1/2,1/2,0) passes d_i+d_j <= 1+d_k, yet K is empty.

    Second external review, item P1-B.  The linear inequalities cut out the
    diagonals of doubly stochastic matrices (i.e. I non-empty); unitary
    feasibility is strictly stronger and is decided by K non-empty.  The
    correct n=3 criterion is Thompson's, on the moduli:
    sqrt(d_i) + sqrt(d_j) <= 1 + sqrt(d_k).
    """
    print()
    print("#" * 78)
    print("# P1-B: d=(1/2,1/2,0) is doubly stochastic feasible but NOT unitary")
    print("#" * 78)
    d = (0.5, 0.5, 0.0)
    triples = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
    lin = [d[i] + d[j] <= 1.0 + d[k] for (i, j, k) in triples]
    tho = [np.sqrt(d[i]) + np.sqrt(d[j]) <= 1.0 + np.sqrt(d[k])
           for (i, j, k) in triples]
    lo, hi = ds_interval(d)
    print(f"  linear inequalities d_i+d_j <= 1+d_k : {lin}  (all True)")
    print(f"  I = [{lo:.12f}, {hi:.12f}]   -> the single point x=0,"
          f" so d IS a doubly stochastic diagonal")
    E = entries(d, np.array([0.0]))
    B = np.array([[float(E[(i, j)][0]) for j in range(3)] for i in range(3)])
    print("  unique doubly stochastic B with that diagonal:")
    for row in B:
        print("     ", "  ".join(f"{v: .3f}" for v in row))
    L = [float(np.sqrt(max(B[0, k] * B[1, k], 0.0))) for k in range(3)]
    print(f"  links of rows 1,2:  L = {[round(v, 6) for v in L]}"
          f"   ->  L3 <= L1+L2 is {L[2] <= L[0] + L[1]}")
    print(f"  Heron G_12(B) = {float(G(d, np.array([0.0]), 0, 1)[0]): .6f}"
          f"   (< 0  =>  NOT unistochastic, K is empty)")
    print(f"  supports of rows 1,2: {{1,3}} and {{2,3}} share only column 3,"
          f" so the overlap has one non-zero term and cannot vanish")
    print(f"  Thompson condition sqrt(d_i)+sqrt(d_j) <= 1+sqrt(d_k) : {tho}"
          f"   (fails: {np.sqrt(0.5)+np.sqrt(0.5):.6f} > 1)")

    # Thompson's criterion vs K non-empty, on a grid of diagonals
    rng = np.random.default_rng(20260910)
    tested = mismatches = 0
    for _ in range(3000):
        dd = tuple(float(v) for v in np.round(rng.random(3), 3))
        if not all(dd[i] + dd[j] <= 1.0 + dd[k] for (i, j, k) in triples):
            continue
        a, b = ds_interval(dd)
        if b <= a + 1e-15:
            continue
        th = all(np.sqrt(dd[i]) + np.sqrt(dd[j]) <= 1.0 + np.sqrt(dd[k])
                 for (i, j, k) in triples)
        _, comps = components_of_K(dd)
        tested += 1
        if th != bool(comps):
            mismatches += 1
            if mismatches <= 5:
                print(f"    mismatch at d={dd}: Thompson={th}, K={bool(comps)}")
    print(f"  Thompson criterion vs (K non-empty): {tested} diagonals tested,"
          f" {mismatches} mismatches")

    # P1-E erratum
    print()
    print("  P1-E erratum:  rho(c,c,c) = 2 sqrt(c)/(1+sqrt(c))")
    for c in (0.4, 2.0 / 3.0):
        print(f"    c = {c:.10f}  ->  rho = {2*np.sqrt(c)/(1+np.sqrt(c)):.16f}")
    print("    (0.7748517734455862, not 0.7746...; 0.8989794855663562 -> 0.899)")


# ----------------------------------------------------------------------
# footnote case: d_3 = 0 makes K a finite set of isolated points
# ----------------------------------------------------------------------
def isolated_points_check():
    """d = (0.04699, 0.59058, 0): K = {x*, S - x*} (two isolated points).

    Pre-submission audit 2026-09-24.  With d_3 = 0 the links of rows (1,3)
    are L = (sqrt(d1 B31), sqrt(x B32), 0), so their triangle inequality
    forces L1 = L2, i.e. d1 (d2 + x) = x (1 - d2 - x), i.e.
        x^2 - (1 - d1 - d2) x + d1 d2 = 0,
    whose two roots sum to S = 1 - d1 - d2 (+ d3 = 0): K is symmetric about
    S/2 (Lemma 3), and the footnote must list BOTH points, not one.
    A uniform grid misses them entirely (G < 0 at every grid point).
    """
    print()
    print("#" * 78)
    print("# footnote: d=(0.04699,0.59058,0)  ->  K = {x*, S-x*}, isolated points")
    print("#" * 78)
    d = (0.04699, 0.59058, 0.0)
    d1, d2, _ = d
    S = 1.0 - d1 - d2 + d[2]
    lo, hi = ds_interval(d)
    disc = (1.0 - d1 - d2) ** 2 - 4.0 * d1 * d2
    assert disc > 0.0
    xs_pts = sorted([0.5 * ((1.0 - d1 - d2) - np.sqrt(disc)),
                     0.5 * ((1.0 - d1 - d2) + np.sqrt(disc))])
    th = np.sqrt(d)
    assert th.sum() - 2.0 * th.min() <= 1.0, "Thompson condition must hold"
    print(f"  S = {S:.12f}   I = [{lo:.6f}, {hi:.6f}]")
    Wv = []
    for x in xs_pts:
        assert lo < x < hi
        xa = np.array([x])
        g = [float(G(d, xa, i, j, ax)[0]) for ax in ("row", "col")
             for (i, j) in ((0, 1), (0, 2), (1, 2))]
        w = float(W_of_x(d, xa)[0])
        Wv.append(w)
        # every chain-link condition is tight (degenerate triangle) at x
        assert max(abs(v) for v in g) < 1e-15, g
        # ... and fails strictly on both sides: x is an isolated point of K
        for h in (1e-6, -1e-6):
            assert float(G(d, np.array([x + h]), 0, 2)[0]) < 0.0
        print(f"  x = {x:.12f}   G12 = {g[0]:+.1e}   max|G| (6 pairs) = "
              f"{max(abs(v) for v in g):.1e}   W = {w:.12f}")
    assert abs(xs_pts[0] + xs_pts[1] - S) < 1e-14, "K symmetric about S/2"
    assert abs(xs_pts[0] - 0.109888) < 1e-6 and abs(xs_pts[1] - 0.252542) < 1e-6
    assert abs(Wv[0] - Wv[1]) < 1e-12 and abs(Wv[0] - 1.151851) < 1e-6
    grid = np.linspace(lo, hi, 200001)
    assert not (G(d, grid, 0, 1) >= 0.0).any(), "grid sampling misses K"
    print("  K = {0.109888..., 0.252542...} = {x*, S - x*},  W = 1.151851 at both;"
          "  a 200001-point grid finds K empty   [OK]")


if __name__ == "__main__":
    cases = [
        ((0.01, 0.01, 0.01), "constant c=0.01  (reviewer counterexample)"),
        ((0.05, 0.05, 0.05), "constant c=0.05   (c < 1/9)"),
        ((0.2, 0.2, 0.2), "constant c=0.2    (c > 1/9)"),
        ((0.5, 1.0 / 3.0, 0.2), "von Mangoldt (1/2,1/3,1/5)"),
        ((0.4, 0.4, 0.4), "constant c=0.4   (Table 1)"),
        ((0.5, 0.45, 0.4), "slightly distinct  (Table 1)"),
        ((0.6, 0.4, 0.4), "nearly degenerate  (Table 1)"),
        ((0.8, 0.5, 0.2), "strongly distinct  (Table 1: infeasible)"),
    ]
    for d, lab in cases:
        report(d, lab)
    constant_family_check()
    reviewer_crosscheck()
    identity_check()
    ds_vs_unitary_check()
    isolated_points_check()
