"""Known kernels reproduced exactly with htk (signs included).

Run:  python3 code/tests/test_known.py          (HTK_FAST=1 skips the |I| = 7 oracle runs)
Each test prints the numbers it compares and raises AssertionError on a mismatch.

Sources: S1 = sources/LLC_main.tex, S2 = 6d-holomorphic-twist-note, S3 = Twisted_11d_sugra_loops,
S4 = check_cd.  Relations used (htk/CONVENTIONS.md):
    I^{S1} = (-1)^{n(n+1)/2} K,   I^{S2} = (-1)^{n(n+1)/2} K,   K^{S3} = (-1)^{n(n+1)/2} K,
    I^{S4}_T = (-1)^{nd + C(n,2)C(d,2)} sigma_T^d K.
"""
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import htk
from htk import brute, catalog, oracle, sp_reduce, PiNumber
from htk import conventions as cv
from htk.graph import polygon, theta, banana, fan, box, graph_A, graph_B, ladder, glued_polygons
from htk.linalg import det

FAST = bool(os.environ.get("HTK_FAST"))
f = math.factorial


def cdet(*cols):
    """Determinant of the matrix with the given columns."""
    n = len(cols)
    return det([[cols[j][i] for j in range(n)] for i in range(n)])


def sdet(*cols):
    return sp.Matrix([[cols[j][i] for j in range(len(cols))] for i in range(len(cols))]).det()


def rvec(rng, H, lo=-5, hi=5):
    return [Fr(rng.randint(lo, hi)) for _ in range(H)]


def rsources(rng, G, H, zero=()):
    lam = [rvec(rng, H) for _ in range(G.nv - 1)]
    for v in zero:
        lam[v] = [Fr(0)] * H
    return lam


def rtimes(rng, G):
    return [Fr(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(G.ne)]


def proj_jacobian(tt, xi, sub):
    """det[t ; dt/dxi_1 ; ...] of a symbolic map tt(xi) at the rational point sub."""
    rows = [list(tt)] + [[sp.diff(x, s) for x in tt] for s in xi]
    M = [[Fr(int(sp.Rational(x.subs(sub)).p), int(sp.Rational(x.subs(sub)).q)) for x in row] for row in rows]
    return det(M)


def toF(x):
    x = sp.Rational(x)
    return Fr(int(x.p), int(x.q))


# ====================================================================== S1: triangle and fans on C^2

def s1_coordinates(red, n):
    """Coordinates r_1..r_2n of S1 (eq-fan-local-ratios, eq-1-25) in terms of the ratios of the
    default plan of fan(n): r_{2i-1} = x of the spoke in the series branch, r_{2i} = w of the
    original spoke (i < n); (r_{2n-1}, r_{2n}) = ((1-w)(1-x), (1-w)x) for the last move.
    Returns (r, |d(r_{2n-1}, r_{2n})/d(w, x)|, substitution of the eliminated ratios)."""
    r, jac = {}, None
    for i, mv in enumerate(red.moves, start=1):
        assert sorted(len(c) for c in mv.chains) == [1, 2]
        ci = [len(c) for c in mv.chains].index(2)
        cs = 1 - ci
        if i < n:
            pos = [j for j, e in enumerate(mv.chains[ci]) if e != 2 * (i - 1)][0]
            r[2 * i - 1] = sp.Symbol("x%d_%d_%d" % (i, ci + 1, pos + 1))
            r[2 * i] = sp.Symbol("w%d_%d" % (i, cs + 1))
        else:
            # last move: chains [(n, n+1)] and [(n, n+2)^red, (n+1, n+2)]
            assert mv.chains[cs] == [2 * (n - 1)] and mv.chains[ci][1] == 2 * n
            w, x = sp.Symbol("w%d_%d" % (i, ci + 1)), sp.Symbol("x%d_%d_2" % (i, ci + 1))
            r[2 * n - 1], r[2 * n] = (1 - w) * (1 - x), (1 - w) * x
            jac = 1 - w
    elim = {sp.Symbol(blk.bary[0]): 1 - sum(sp.Symbol(b) for b in blk.indep) for blk in red.blocks}
    return r, jac, elim


def test_S1_fans():
    """S1 Theorem thm-fan-kernel, eq. (eq-1-42):
    I(lam; 0)|_{lam_2..lam_n = 0} = (lam_1 . lam_{n+1})^n / ((n!)^3 (n+1)),  lam_1.lam_2 = det(lam_1, lam_2)."""
    a1, a2, b1, b2 = sp.symbols("a1 a2 b1 b2")
    rng = random.Random(101)
    for n in (1, 2, 3):
        G = fan(n)
        red = sp_reduce(G)
        lam = [[0, 0] for _ in range(G.nv - 1)]
        lam[0] = [a1, a2]
        lam[n] = [b1, b2]
        K = red.kernel(lam)
        I = cv.S1_kernel(K, n)
        expected = (a1 * b2 - a2 * b1) ** n / sp.Integer(f(n) ** 3 * (n + 1))
        assert sp.expand(I - expected) == 0, (n, I)
        # general sources: local determinants, generic Jacobian and the brute-force layer agree
        lamg = rsources(rng, G, 2)
        Kl = red.kernel(lamg, method="local", check=3)
        Kg = red.kernel(lamg, method="generic", check=3)
        assert Kl == Kg
        # density of S1 eq. (eq-1-41) in the coordinates r_1..r_2n of S1 (eq-fan-local-ratios)
        g = red.density(lam, check=2).as_expr(factored=False)
        r, jac, elim = s1_coordinates(red, n)
        s1 = (a1 * b2 - a2 * b1) ** n * r[2 * n] ** (n - 1)
        for j in range(1, n):
            s1 *= r[2 * j - 1] ** (n - j) * r[2 * j] * (1 - r[2 * j]) ** (j - 1)
        assert sp.expand((cv.sign_n(n) * g - s1 * jac).subs(elim)) == 0, n
        print("  fan n=%d: I^{S1} = (-1)^{n(n+1)/2} K = %s ; general sources K = %s" % (n, sp.factor(I), Kl))


def test_S1_fan_general_sources():
    """S1 Theorem thm-fan-kernel with sources on all path vertices, n = 1..4:
    density (eq-1-40), path currents (eq-1-38), spoke currents (eq-1-39), as polynomial identities
    in the coordinates r of S1."""
    rng = random.Random(114)
    d2 = lambda u, v: u[0] * v[1] - u[1] * v[0]
    for n in (1, 2, 3, 4):
        G = fan(n)
        red = sp_reduce(G)
        lam = rsources(rng, G, 2)
        g = red.density(lam, check=2).as_expr(factored=False)
        r, jac, elim = s1_coordinates(red, n)
        L = [sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in v]) for v in lam]
        lamred = {}
        for i in range(1, n + 1):
            acc = sp.zeros(2, 1)
            for j in range(1, i + 1):
                prod = sp.Integer(1)
                for k in range(j, i):
                    prod *= r[2 * k - 1]
                acc += prod * L[j - 1]
            lamred[i] = acc
        y = {n: r[2 * n - 1] * lamred[n] - r[2 * n] * L[n]}
        for i in range(n - 1, 0, -1):
            y[i] = r[2 * i - 1] * r[2 * i] * lamred[i] - (1 - r[2 * i]) * L[i] + (1 - r[2 * i]) * y[i + 1]
        s1 = d2(lamred[n], L[n])
        for i in range(1, n):
            s1 *= r[2 * i] * d2(lamred[i], L[i] - y[i + 1])
        assert sp.expand((cv.sign_n(n) * g - s1 * jac).subs(elim)) == 0, n
        cur = red.currents(lam)
        for i in range(1, n + 1):                       # y_{i,i+1} = edge 2(i-1)
            for a in range(2):
                assert sp.expand((cur[2 * (i - 1)][a].as_expr() - y[i][a]).subs(elim)) == 0
        for i in range(1, n + 2):                       # spokes: edge 2(i-1)+1 (i <= n), edge 2n (i = n+1)
            e = 2 * (i - 1) + 1 if i <= n else 2 * n
            want = L[i - 1] + (y[i - 1] if i > 1 else sp.zeros(2, 1)) - (y[i] if i <= n else sp.zeros(2, 1))
            for a in range(2):
                assert sp.expand((cur[e][a].as_expr() - want[a]).subs(elim)) == 0
    print("  n = 1..4: density (eq-1-40) and currents (eq-1-38), (eq-1-39) reproduced")


def test_S1_triangle_displacements():
    """S1 eq. (eq-1-46): I_triangle(lam; z) = det(lam_1, lam_2) exp(-lam_1.z_13 - lam_2.z_23)
    sum_{i,j} [-lam_1.Z]^i [lam_2.Z]^j / (i+j+2)!,  Z = z_12 + z_23 - z_13; orders 0..3 in z."""
    rng = random.Random(102)
    G = fan(1)                                  # edges (12), (13), (23) of S1 = (0,1), (0,2), (1,2)
    red = sp_reduce(G)
    lam = rsources(rng, G, 2)
    z12, z13, z23 = (rvec(rng, 2), rvec(rng, 2), rvec(rng, 2))
    order = 3
    Ks = red.taylor(lam, [z12, z13, z23], order)
    dot = lambda u, v: sum(x * y for x, y in zip(u, v))
    Z = [z12[i] + z23[i] - z13[i] for i in range(2)]
    eps = sp.Symbol("eps")                       # z -> eps z
    A = -dot(lam[0], z13) - dot(lam[1], z23)
    series = sp.Integer(0)
    for i in range(order + 1):
        for j in range(order + 1 - i):
            series += (-dot(lam[0], Z) * eps) ** i * (dot(lam[1], Z) * eps) ** j / sp.Integer(f(i + j + 2))
    expr = sp.series(cdet(lam[0], lam[1]) * sp.exp(A * eps) * series, eps, 0, order + 1).removeO()
    for m in range(order + 1):
        want = toF(expr.coeff(eps, m))
        got = cv.S1_kernel(Ks[m], 1)
        assert got == want, (m, got, want)
    print("  triangle: Taylor orders 0..%d in z: %s" % (order, [str(cv.S1_kernel(k, 1)) for k in Ks]))


# ====================================================================== S4: polygons, chain move, ladder

def test_S4_polygons():
    """S4 check_single_polygon: (d+1)-gon, cotree e_0: density of d^d y_{e_0} = (-1)^d det(lam_1..lam_d)
    on the simplex, |I(lam; 0)| = |det| / d!.   htk: K = (-1)^d det(lam_0..lam_{d-1}) / d!."""
    rng = random.Random(103)
    for d in (2, 3, 4, 5):
        G = polygon(d + 1)
        while True:
            lam = rsources(rng, G, d)
            if cdet(*lam) != 0:
                break
        lam_root = [-sum(lam[v][i] for v in range(d)) for i in range(d)]
        full = lam + [lam_root]
        red = sp_reduce(G)
        K = red.kernel(lam, check=3)
        assert K == Fr((-1) ** d) * cdet(*lam) / f(d)
        I4 = cv.S4_integral(K, G, [0], d)
        want = Fr((-1) ** d) * cdet(*[full[v] for v in range(1, d + 1)]) / f(d)
        assert I4 == want
        # pointwise statement of S4 on the simplex sum t = 1
        for _ in range(3):
            t = rtimes(rng, G)
            Ddy = brute.jacobian_density(G, t, lam, [0])
            assert Ddy * sum(t) ** (d + 1) == Fr((-1) ** d) * cdet(*[full[v] for v in range(1, d + 1)])
        # oracle
        val, err = oracle.kernel(G, lam, method="gauss")
        assert abs(val - float(K)) <= 1e-10 * max(1.0, abs(float(K)))
        print("  d=%d: K = %s, int d^d y_{e_0} = %s = (-1)^d det(lam_1..lam_d)/d!; oracle(gauss) %.12g" % (d, K, I4, val))


def test_S4_chain_move():
    """S4 check_glued: two (d+1)-gons glued on an edge, cotree (e_1, eps_2):
    density = rho^{d-1} det(y_red, lam_v1..lam_v(d-1)) * (-1)^d det(lam_w^red, lam_p2..lam_pd)
    on (0,1) x {1 > s_1 > .. > s_{d-1} > 0} x {1 > q_1 > .. > q_d > 0}; its integral equals
    the htk value of int d^d y_{c_1} ^ d^d y_{c_2}."""
    rng = random.Random(104)
    for d in ((2, 3, 4) if FAST else (2, 3, 4, 5)):
        G, cot, chain, base = glued_polygons(d)
        lam = [[Fr(rng.randint(-4, 4)) for _ in range(d)] for _ in range(G.nv)]
        tot = [sum(lam[v][c] for v in range(G.nv)) for c in range(d)]
        lam[base[-1]] = [lam[base[-1]][c] - tot[c] for c in range(d)]
        assert base[-1] == G.root
        red = sp_reduce(G)
        K = red.kernel(lam[:-1], check=2)
        I4 = cv.S4_integral(K, G, cot, d)
        # S4 formula
        rho = sp.Symbol("rho")
        sig = sp.symbols("s1:%d" % d)
        q = sp.symbols("q1:%d" % (d + 1))
        vec = lambda v: sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in lam[v]])
        lch = [vec(chain[j]) for j in range(1, d)]
        lw = vec(base[1]) + sum(((1 - sig[j]) * lch[j] for j in range(d - 1)), sp.zeros(d, 1))
        lp = [lw] + [vec(base[j]) for j in range(2, d + 1)]
        yred = -sum((q[j] * lp[j] for j in range(d)), sp.zeros(d, 1))
        pred = sp.expand(rho ** (d - 1) * sp.Matrix.hstack(yred, *lch).det() * (-1) ** d * sp.Matrix.hstack(*lp).det())
        P = sp.Poly(pred, rho, *sig, *q)
        total = sp.Integer(0)
        for mon, c in P.terms():
            a, ms, mq = mon[0], mon[1:d], mon[d:]
            val = sp.Rational(1, a + 1)
            # ordered simplex 1 > s_1 > ... > s_k > 0: integrate from the innermost variable
            for block in (ms, mq):
                acc = 0
                for kk in range(len(block) - 1, -1, -1):
                    acc += block[kk] + 1
                    val *= sp.Rational(1, acc)
            total += c * val
        assert toF(total) == I4, (d, total, I4)
        print("  d=%d: int d^d y ^ d^d y = %s (S4 chain-move density integrated: %s)" % (d, I4, total))


def ladder_closed_form(n):
    return Fr((-1) ** (n - 1) * 2 ** (n + 1), f(n + 2) * f(n) ** 2 * f(n + 1) ** 2)


def test_S4_ladder():
    """S4 check_ladder: n-box ladder on C^3, sources A, B, C at a_0, b_0, b_n, Delta = det(A, B, C):
    I / Delta^n = (-1)^{n-1} 2^{n+1} / ((n+2)! n!^2 (n+1)!^2) = 1/6, -1/432, 1/155520, -1/186624000."""
    rng = random.Random(105)
    values = []
    for n in (1, 2, 3, 4):
        G, cot, a, b = ladder(n)
        while True:
            A, B, C = (rvec(rng, 3, -4, 4) for _ in range(3))
            Delta = cdet(A, B, C)
            if Delta != 0:
                break
        lam = {a(0): A, b(0): B, b(n): C}
        if a(n) != G.root:
            lam[a(n)] = [-(A[i] + B[i] + C[i]) for i in range(3)]
        if b(n) == G.root:
            del lam[b(n)]
            lam[a(n)] = [-(A[i] + B[i] + C[i]) for i in range(3)]
        red = sp_reduce(G)
        K = red.kernel(lam, check=2)
        I4 = cv.S4_integral(K, G, cot, 3)
        assert I4 / Delta ** n == ladder_closed_form(n), (n, I4 / Delta ** n)
        values.append(I4 / Delta ** n)
    assert values == [Fr(1, 6), Fr(-1, 432), Fr(1, 155520), Fr(-1, 186624000)]
    print("  I/Delta^n, n = 1..4:", ", ".join(str(v) for v in values))


# ====================================================================== S2: box, (A), (B) on C^3

def s2_box_times(r):
    full = [1] + list(r) + [0]
    return [full[i] - full[i + 1] for i in range(4)]


def test_S2_box():
    """S2 eq. (eq-box-kernel), (eq-box-zero): I_box(lam; 0) = det(lam_1, lam_2, lam_3) / 6,
    -sum_e y_e.z_e = sum_j lam_j.(r_j z_1 + (1 - r_j) z_2), r_j = sum_{k>j} t_k / T;
    eps_T det B_T = -1 for T = (e_2, e_3, e_4)."""
    rng = random.Random(106)
    G = box()
    assert G.sigma([1, 2, 3]) == -1
    L = [list(sp.symbols("l%d_1:4" % j)) for j in (1, 2, 3)]
    red = sp_reduce(G)
    K = red.kernel(L)
    I = cv.S2_kernel(K, 1)
    assert sp.expand(I - sdet(*L) / 6) == 0
    # displacement map as a function of the Schwinger times
    ts = sp.symbols("t1:5", positive=True)
    Tt = sum(ts)
    r = [sum(ts[k] for k in range(j + 1, 4)) / Tt for j in range(3)]
    ratios = red.ratios_of_t(list(ts))
    sub = {sp.Symbol(c): ratios[c] for c in red.coords}
    X = red.displacement_expr()
    for j in range(3):
        # z_{e_1} = z_1, z_{e_4} = -z_2  =>  X_j = X[j][e_1] z_1 - X[j][e_4] z_2
        assert sp.cancel(X[j][0].subs(sub) - r[j]) == 0
        assert sp.cancel(-X[j][3].subs(sub) - (1 - r[j])) == 0
    # density in the coordinates r of S2:  I-density = det(lam_1, lam_2, lam_3) d^3 r
    rs = sp.symbols("r1:4")
    tt = s2_box_times(rs)
    for _ in range(3):
        lam = rsources(rng, G, 3)
        vals = sorted(rng.sample(range(1, 1000), 3), reverse=True)
        sub_r = {s: sp.Rational(v, 1000) for s, v in zip(rs, vals)}
        t = [toF(x.subs(sub_r)) for x in tt]
        J = proj_jacobian(tt, rs, sub_r)
        dens = cv.S2_kernel(brute.density(G, t, lam).rational() * abs(J), 1)
        assert dens == cdet(*lam)
    # first-order term, S2 eq. (eq-Q1-deriv): int d^3r (r_1, r_2, r_3) = (1/8, 1/12, 1/24)
    lam = rsources(rng, G, 3)
    z1, z2 = rvec(rng, 3), rvec(rng, 3)
    K0, K1 = red.taylor(lam, {0: z1, 3: [-x for x in z2]}, 1)
    c = [Fr(1, 8), Fr(1, 12), Fr(1, 24)]
    dot = lambda u, v: sum(x * y for x, y in zip(u, v))
    want1 = cdet(*lam) * sum(c[j] * dot(lam[j], z1) + (Fr(1, 6) - c[j]) * dot(lam[j], z2) for j in range(3))
    assert cv.S2_kernel(K0, 1) == cdet(*lam) / 6 and cv.S2_kernel(K1, 1) == want1
    print("  box: I^{S2} = -K = %s ; displacement map and first-order term agree" % sp.factor(I))


def s2_times_A(r):
    r1, r2, r3, r4, r5, r6 = r
    tau = s2_box_times((r4, r5, r6))
    T = tau[0] / (1 - r1)
    return [T * (1 - r2), T * (r2 - r3), T * r3, tau[0] / r1, tau[1], tau[2], tau[3]]


def s2_times_B(r):
    r1, r2, r3, r4, r5, r6 = r
    tau = s2_box_times((r4, r5, r6))
    T1, T2 = tau[0] / (1 - r1), tau[0] / r1
    return [T1 * (1 - r2), T1 * r2, T2 * (1 - r3), T2 * r3, tau[1], tau[2], tau[3]]


def s2_ratios_A(t):
    T = t[0] + t[1] + t[2]
    tred = 1 / (1 / T + 1 / t[3])
    Tb = tred + t[4] + t[5] + t[6]
    return [T / (T + t[3]), (t[1] + t[2]) / T, t[2] / T, (t[4] + t[5] + t[6]) / Tb, (t[5] + t[6]) / Tb, t[6] / Tb]


def s2_ratios_B(t):
    T1, T2 = t[0] + t[1], t[2] + t[3]
    tred = 1 / (1 / T1 + 1 / T2)
    Tb = tred + t[4] + t[5] + t[6]
    return [T1 / (T1 + T2), t[1] / T1, t[3] / T2, (t[4] + t[5] + t[6]) / Tb, (t[5] + t[6]) / Tb, t[6] / Tb]


def _s2_two_loop(name, G, s2_tree, s2_sigma, times_fn, ratios_fn, JS2, closed, Xmap, zedges, ordered_first):
    rng = random.Random(107)
    assert G.sigma(s2_tree) == s2_sigma
    red = sp_reduce(G)
    # (i) zero-displacement value with general sources at the four output vertices (lam_3 = 0)
    names = ["l%d" % v for v in (0, 1, 3, 4)]
    L = {v: list(sp.symbols("%s_1:4" % nm)) for v, nm in zip((0, 1, 3, 4), names)}
    K = red.kernel(L)
    I = cv.S2_kernel(K, 2)
    assert sp.expand(I - closed(L)) == 0
    # (ii) density in the coordinates r of S2 with sources at all five vertices
    rs = sp.symbols("r1:7")
    tt = times_fn(rs)
    for _ in range(3):
        lam = rsources(rng, G, 3)
        v = {}
        v[rs[0]] = sp.Rational(rng.randint(1, 99), 100)
        if ordered_first:
            a, b = sorted(rng.sample(range(1, 1000), 2), reverse=True)
            v[rs[1]], v[rs[2]] = sp.Rational(a, 1000), sp.Rational(b, 1000)
        else:
            v[rs[1]], v[rs[2]] = sp.Rational(rng.randint(1, 99), 100), sp.Rational(rng.randint(1, 99), 100)
        c = sorted(rng.sample(range(1, 1000), 3), reverse=True)
        for s, x in zip(rs[3:], c):
            v[s] = sp.Rational(x, 1000)
        t = [toF(x.subs(v)) for x in tt]
        J = proj_jacobian(tt, rs, v)
        dens = cv.S2_kernel(brute.density(G, t, lam).rational() * abs(J), 2)
        want = s2_sigma * toF(JS2(lam, [v[s] for s in rs]))
        assert dens == want, (name, dens, want)
    # (iii) displacement map X_v of S2 eq. (eq-res-X), as functions of the Schwinger times
    ts = sp.symbols("t1:8", positive=True)
    r = ratios_fn(list(ts))
    ratios = red.ratios_of_t(list(ts))
    sub = {sp.Symbol(c): ratios[c] for c in red.coords}
    X = red.displacement_expr()
    z = sp.symbols("z1:4")
    W = (1 - r[0]) * z[0] + r[0] * z[1] - z[2]
    want = Xmap(r, z, W)
    for vtx, expr in want.items():
        mine = sum(X[vtx][e] * ze for e, ze in zedges(z).items()).subs(sub)
        assert sp.cancel(sp.together(mine - expr)) == 0, (name, vtx)
    print("  %s: I^{S2} = -K = %s" % (name, sp.factor(I)))
    return I


def test_S2_graph_A():
    """S2 eq. (eq-IA), (eq-res-J0), (eq-res-X): eps_T det B_T = +1,
    I_A(lam; 0) = (2<124> + <125>)(<145> + 2<245>) / 432,
    density J_A = r_1^2 det(y_red, lam_1, lam_2) (-det(lam_3^red, lam_4, lam_5))."""
    def closed(L):
        l1, l2, l4, l5 = L[0], L[1], L[3], L[4]
        return (2 * sdet(l1, l2, l4) + sdet(l1, l2, l5)) * (sdet(l1, l4, l5) + 2 * sdet(l2, l4, l5)) / 432

    def JA(lam, r):
        r1, r2, r3, r4, r5, r6 = r
        l1, l2, l3, l4, l5 = (sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in v]) for v in lam)
        l3r = l3 + (1 - r2) * l1 + (1 - r3) * l2
        yred = -(r4 * l3r + r5 * l4 + r6 * l5)
        return r1 ** 2 * sp.Matrix.hstack(yred, l1, l2).det() * (-sp.Matrix.hstack(l3r, l4, l5).det())

    def Xmap(r, z, W):
        return {0: r[1] * z[0] + (1 - r[1]) * (z[2] + r[3] * W), 1: r[2] * z[0] + (1 - r[2]) * (z[2] + r[3] * W),
                3: z[2] + r[4] * W, 4: z[2] + r[5] * W}

    # displacements: z_1 on (O->1) = edge 0, z_2 on (O->3) = edge 3, -z_3 on (5->O) = edge 6
    _s2_two_loop("(A)", graph_A(), [1, 2, 3, 5, 6], +1, s2_times_A, s2_ratios_A, JA, closed, Xmap,
                 lambda z: {0: z[0], 3: z[1], 6: -z[2]}, True)


def test_S2_graph_B():
    """S2 eq. (eq-IB), (eq-res-J0), (eq-res-X): eps_T det B_T = -1,
    I_B(lam; 0) = -(2<uwa> + <uwb>)(<uab> + <wab>) / 288,
    density J_B = r_1 (1 - r_1) det(y_red, lam_u, lam_w) (-det(lam_q^red, lam_a, lam_b))."""
    def closed(L):
        lu, lw, la, lb = L[0], L[1], L[3], L[4]
        return -(2 * sdet(lu, lw, la) + sdet(lu, lw, lb)) * (sdet(lu, la, lb) + sdet(lw, la, lb)) / 288

    def JB(lam, r):
        r1, r2, r3, r4, r5, r6 = r
        lu, lw, lq, la, lb = (sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in v]) for v in lam)
        lqr = lq + (1 - r2) * lu + (1 - r3) * lw
        yred = -(r4 * lqr + r5 * la + r6 * lb)
        return r1 * (1 - r1) * sp.Matrix.hstack(yred, lu, lw).det() * (-sp.Matrix.hstack(lqr, la, lb).det())

    def Xmap(r, z, W):
        return {0: r[1] * z[0] + (1 - r[1]) * (z[2] + r[3] * W), 1: r[2] * z[1] + (1 - r[2]) * (z[2] + r[3] * W),
                3: z[2] + r[4] * W, 4: z[2] + r[5] * W}

    # displacements: z_1 on (O->u) = edge 0, z_2 on (O->w) = edge 2, -z_3 on (b->O) = edge 6
    _s2_two_loop("(B)", graph_B(), [1, 2, 3, 5, 6], -1, s2_times_B, s2_ratios_B, JB, closed, Xmap,
                 lambda z: {0: z[0], 2: z[1], 6: -z[2]}, False)


# ====================================================================== S3: R x C^5

C_PQR = {(1, 6, 6): Fr(1, 158400), (2, 5, 6): Fr(1, 285120), (3, 4, 6): Fr(1, 332640),
         (3, 5, 5): Fr(1, 798336), (4, 4, 5): Fr(1, 997920)}
S3_C = [1, 4, 4, 4, 4, 4]
S3_K = [2, 4, 4, 4, 4, 4]


def test_S3_extraction_signs():
    """S3 Appendix A item 3: varsigma^ext = -1 for C_7, +1 for the five theta graphs;
    varsigma = (-1)^n (-1)^{(5n+1)5n/2} varsigma^ext = -1 for all six graphs (Theorems thm-hept, thm-theta)."""
    G = polygon(7)
    ext = cv.extraction_sign(G, 6, c=S3_C, k=S3_K)
    assert ext == -1 == cv.extraction_sign_formula(G, 6)
    assert cv.varsigma_S3(1, ext) == -1 == cv.sign_n(1)
    out = ["C7: %+d" % ext]
    for pqr in C_PQR:
        G, _ = theta(*pqr)
        ext = cv.extraction_sign(G, 6, c=S3_C, k=S3_K)
        assert ext == +1 == cv.extraction_sign_formula(G, 6)
        assert cv.varsigma_S3(2, ext) == -1 == cv.sign_n(2)
        # D = 6 is even: the same sign for S3's cycle basis (A - B, B - C)
        p, q, r = pqr
        C3 = [[1] * p + [-1] * q + [0] * r, [0] * p + [1] * q + [-1] * r]
        assert cv.extraction_sign(G, 6, C=C3, c=S3_C, k=S3_K) == +1
        out.append("theta%s: %+d" % (pqr, ext))
    print("  varsigma^ext:", ", ".join(out))


def test_S3_heptagon():
    """S3 Theorem thm-hept: with q_v = (kappa_v; lam_v),
    U^{-11/2} <psi^tau> ^ prod_i psi^i = det[q_1..q_6] U^{-13/2} iota_E prod dt_e (times exp(sum t w^2)),
    K_7(lam, kappa; 0) = -det[q_1..q_6]/6! (1 + O(kappa^2)); currents eq. (eq-hept-currents)."""
    rng = random.Random(108)
    G = polygon(7)
    lam = rsources(rng, G, 5)
    kappa = rvec(rng, 6)
    q6 = [[kappa[v]] + lam[v] for v in range(6)]
    detq = cdet(*q6)
    for _ in range(3):
        t = rtimes(rng, G)
        U = sum(t)
        D = brute.density(G, t, lam, T=1, kappa=kappa, check_basic=True)
        assert D.coeff * U ** 6 == detq and D.radicand == U
        wbar = [row[0] for row in brute.currents(G, t, [[k] for k in kappa])]
        assert D.expo == sum(te * w * w for te, w in zip(t, wbar))
        # the same density from the T = 0 kernel on C^6 with sources q
        assert brute.density(G, t, q6, T=0).rational() == D.coeff / U
        # currents: J_{e_j} = Lam_j - Lambar, J_{e_7} = -Lambar
        J = brute.currents(G, t, q6)
        Lam = [[sum(q6[v][a] for v in range(j + 1)) for a in range(6)] for j in range(6)]
        Lbar = [sum(t[k] * Lam[k][a] for k in range(6)) / U for a in range(6)]
        assert all(J[j] == [Lam[j][a] - Lbar[a] for a in range(6)] for j in range(6))
        assert J[6] == [-x for x in Lbar]
        # no topological sources: rho_1 = 0
        assert brute.density(G, t, lam, T=1).coeff == 0
    K6 = sp_reduce(G).kernel(q6, check=2)             # int over the simplex of det[q] = det[q]/6!
    assert K6 == detq / 720
    K7 = cv.S3_kernel(K6, 1)                           # varsigma_{C_7} = -1
    assert K7 == -detq / 720
    print("  heptagon: det[q_1..q_6] = %s, K_7 = varsigma * det[q]/6! = %s" % (detq, K7))


def test_S3_thetas():
    """S3 Theorem thm-theta, eq. (R-c): K_theta(lam; 0) = c_pqr pi D_theta(lam),
    c_166 = 1/158400, c_256 = 1/285120, c_346 = 1/332640, c_355 = 1/798336, c_445 = 1/997920;
    integrand = -(1/2) D_theta prod_P w_P^{11/2 - n_P} |dx d^2w|  (Appendix A item 4)."""
    rng = random.Random(109)
    for pqr, c in C_PQR.items():
        G, interior = theta(*pqr)
        lam = rsources(rng, G, 5)
        Dth = brute.theta_D(lam, interior)
        red = sp_reduce(G, prefer="large")
        assert len(red.moves) == 1 and red.moves[0].chains == [list(range(0, pqr[0])),
                                                                list(range(pqr[0], pqr[0] + pqr[1])),
                                                                list(range(pqr[0] + pqr[1], sum(pqr)))]
        dens = red.density(lam, T=1, check=2)
        # monomial density of S3
        assert dens.num.is_ground and Fr(int(dens.num.LC.numerator), int(dens.num.LC.denominator)) == -Dth / 2
        expo = dens.exponents()
        assert all(expo.get("w1_%d" % (i + 1), Fr(0)) == Fr(11, 2) - pqr[i] for i in range(3))
        assert set(expo) <= {"w1_1", "w1_2", "w1_3"}
        K = dens.integrate()
        KS3 = cv.S3_kernel(K, 2)
        assert KS3 == PiNumber({1: c * Dth}), (pqr, KS3)
        # closed form of S3 eq. (eq-Ktheta0)
        cc = sp.gamma(sp.Rational(13, 2) - pqr[0]) * sp.gamma(sp.Rational(13, 2) - pqr[1]) * sp.gamma(sp.Rational(13, 2) - pqr[2]) \
            / (2 * sp.gamma(sp.Rational(13, 2)) * f(pqr[0] - 1) * f(pqr[1] - 1) * f(pqr[2] - 1))
        assert sp.simplify(cc - sp.Rational(c.numerator, c.denominator) * sp.pi) == 0
        print("  theta%s: D_theta = %s, K^{S3} = -K = %s  (c = 1/%d)" % (pqr, Dth, KS3, c.denominator))


def test_S3_theta_currents_and_insertions():
    """S3 eq. (eq-series), (eq-parallel): y^P_j = l_P + Q^P_j, l_P = w_P (lam_N + sum_P' mhat_P') - mhat_P;
    eq. (eq-insertion): insertions are Dirichlet integrals (checked against the direct formula)."""
    rng = random.Random(110)
    pqr = (2, 5, 6)
    G, interior = theta(*pqr)
    lam = rsources(rng, G, 5)
    red = sp_reduce(G, prefer="large")
    pt = red.random_point(rng)
    y = red.currents_at(lam, pt)
    mh, Qs = [], []
    off = 0
    for i, n_ in enumerate(pqr):
        x = [pt["x1_%d_%d" % (i + 1, j + 1)] for j in range(n_)] if n_ > 1 else [Fr(1)]
        Q = [[Fr(0)] * 5]
        for k in range(1, n_):
            Q.append([Q[-1][a] + lam[interior[i][k - 1]][a] for a in range(5)])
        mh.append([sum(x[j] * Q[j][a] for j in range(n_)) for a in range(5)])
        Qs.append(Q)
    for i, n_ in enumerate(pqr):
        w = pt["w1_%d" % (i + 1)]
        ell = [w * (lam[0][a] + sum(m[a] for m in mh)) - mh[i][a] for a in range(5)]
        for j in range(n_):
            assert y[off + j] == [ell[a] + Qs[i][j][a] for a in range(5)]
        off += n_
    # one insertion: int g y^1_{e_0} = -(D_theta/2) int w-monomial * l_A^1 ... compare with brute structure
    dens = red.density(lam, T=1, check=1)
    val = red.integrate(lam, T=1, powers={(0, 0): 1}, dens=dens)
    Dth = brute.theta_D(lam, interior)
    # direct: y^1_{A_1} = w_A (lam_N + sum mhat) - mhat_A, mhat_P = sum_j x_j Q_j; means of x over the fibre: 1/n_P
    from htk.simplex import dirichlet
    base = [Fr(11, 2) - n_ for n_ in pqr]
    I0, k0 = dirichlet(tuple(base))
    I1, k1 = dirichlet((base[0] + 1, base[1], base[2]))
    vol = Fr(1)
    for n_ in pqr:
        vol /= f(n_ - 1)
    mbar = [[sum(Qs[i][j][a] for j in range(pqr[i])) / pqr[i] for a in range(5)] for i in range(3)]
    direct = -Dth / 2 * vol * (I1 * (lam[0][0] + sum(m[0] for m in mbar)) - I0 * mbar[0][0])
    assert k0 == k1 == 1 and val == PiNumber({1: direct})
    print("  theta%s: polynomial currents agree; <y^1_{A_1}> insertion = %s" % (pqr, val))


def test_rho_odd():
    """S3 Proposition prop-rho / Corollary cor-odd: rho_n = 0 for odd n."""
    rng = random.Random(111)
    cases = [(polygon(7), 5), (polygon(3), 1), (fan(3), 1)]
    G33 = catalog.tight_graphs(3, 3)[0]
    cases.append((G33, 2))
    for G, H in cases:
        assert G.loops % 2 == 1
        t = rtimes(rng, G)
        form, _ = brute.topological_form(G, t)
        assert form == {}
        assert brute.density(G, t, rsources(rng, G, H), T=1).coeff == 0
    print("  rho_n = 0 on %d graphs with n = 1, 1, 3, 3" % len(cases))


def test_S3_H0_remark():
    """S3 remark after Corollary cor-odd: H = 0, n = 2, three parallel edges:
    <psi_1 psi_2> = Omega_R / (2U), raw integral (1/4pi) int Omega_R / (2 U^{3/2})
    = (1/4pi) Gamma(1/2)^3 / (2 Gamma(3/2)) = 1/4 = 3 int P^2 dP, P = sgn(tau)/2."""
    rng = random.Random(112)
    G = banana(3)
    C3 = [[1, -1, 0], [0, 1, -1]]                      # cycles A - B, B - C of S3
    for _ in range(3):
        t = rtimes(rng, G)
        U = t[0] * t[1] + t[1] * t[2] + t[2] * t[0]
        D = brute.density(G, t, None, T=1, C=C3, check_basic=True)
        assert D.radicand == U and D.coeff == 1 / (2 * U)          # Omega_R / (2 U^{3/2}) = iota_E vol / (2 U^{3/2})
        assert brute.density(G, t, None, T=1).same_as(D)           # the positive basis gives the same sign here
    red = sp_reduce(G, D=1)
    K = red.kernel(None, T=1, check=3)
    assert K == PiNumber({1: 1})
    raw = K.coeff(1) / 4                                # (1/(4 pi)) * pi
    moyal = 3 * 2 * Fr(1, 3) * Fr(1, 2) ** 3            # 3 int_{-1/2}^{1/2} P^2 dP
    assert raw == Fr(1, 4) == moyal
    g = sp.gamma(sp.Rational(1, 2)) ** 3 / (2 * sp.gamma(sp.Rational(3, 2)))
    assert sp.simplify(g - sp.pi) == 0
    val, err = oracle.kernel(G, None, T=1)
    assert abs(val - math.pi) < 1e-6
    print("  banana, H = 0: K = %s, (1/4pi) K = %s; oracle %.10f +- %.1e" % (K, raw, val, err))


# ====================================================================== numerical oracle

def test_oracle():
    """Independent numerical integration (lattice rule in log coordinates) against the exact values."""
    rng = random.Random(113)
    cases = [(fan(1), 2, 0, 1e-8), (fan(2), 2, 0, 1e-5), (box(), 3, 0, 1e-6),
             (theta(1, 2, 2)[0], 1, 1, 1e-5)]
    if not FAST:
        cases += [(fan(3), 2, 0, 1e-3), (graph_A(), 3, 0, 1e-3), (graph_B(), 3, 0, 1e-3)]
    for G, H, T, rtol in cases:
        while True:
            lam = rsources(rng, G, H)
            exact = sp_reduce(G, prefer="large" if T else "small").kernel(lam, T, check=1)
            ex = float(exact)
            if ex != 0:
                break
        t0 = time.time()
        val, err = oracle.kernel(G, lam, T)
        rel = abs(val - ex) / abs(ex)
        assert rel < rtol, (G.name, val, ex, rel)
        print("  %-12s T=%d |I|=%d: exact %-22s oracle %.10g  rel.err %.1e  (est. %.1e, %.0fs)" % (
            G.name, T, G.ne, exact, val, rel, err / abs(ex), time.time() - t0))
    # one Taylor coefficient: box, coefficient of z_{e_1, 1}
    G = box()
    lam = rsources(rng, G, 3)
    red = sp_reduce(G)
    exact = -red.integrate(lam, powers={(0, 0): 1})
    val, err = oracle.taylor_coefficient(G, lam, {(0, 0): 1})
    assert abs(val - float(exact)) < 1e-6 * max(1.0, abs(float(exact)))
    print("  box: coefficient of z_{e_1}^1: exact %s, oracle %.10g" % (exact, val))


TESTS = [test_S1_fans, test_S1_fan_general_sources, test_S1_triangle_displacements, test_S4_polygons, test_S4_chain_move, test_S4_ladder,
         test_S2_box, test_S2_graph_A, test_S2_graph_B, test_S3_extraction_signs, test_S3_heptagon,
         test_S3_thetas, test_S3_theta_currents_and_insertions, test_rho_odd, test_S3_H0_remark, test_oracle]


def main(argv):
    names = set(argv[1:])
    failed, passed = [], []
    for t in TESTS:
        if names and t.__name__ not in names:
            continue
        t0 = time.time()
        print("%s" % t.__name__, flush=True)
        try:
            t()
            passed.append(t.__name__)
            print("  PASS (%.1fs)" % (time.time() - t0), flush=True)
        except Exception as exc:                                  # noqa: BLE001
            import traceback
            traceback.print_exc()
            failed.append(t.__name__)
            print("  FAIL: %r" % (exc,), flush=True)
    print("passed %d, failed %d" % (len(passed), len(failed)))
    if failed:
        print("failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
