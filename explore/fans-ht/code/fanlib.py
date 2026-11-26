"""Fan kernels at (H,T) = (1,1) and (2,0) in the coordinates r_1..r_2n of S1 (Section "Fan kernels").

S1 labels: path vertices 1..n+1, insertion O = n+2 (root).  Resistances: a_i = t_{i,i+1} (i=1..n),
b_i = t_{i,O} (i=1..n+1).  Engine edge order (htk.fan(n)): a_1, b_1, a_2, b_2, ..., a_n, b_n, b_{n+1}.

u_k = r_{2k-1}, v_k = r_{2k}  (k = 1..n).   Domain: (0,1)^{2n-2} x {u_n, v_n > 0, u_n + v_n < 1}.
Inverse map (gauge a_n + B_n + b_{n+1} = 1, S1 eq-fan-local-ratios, eq-1-25):
    B_n = u_n, b_{n+1} = v_n, a_n = 1 - u_n - v_n;
    i = n-1..1:  R_i = B_{i+1}/(1 - v_i),  b_{i+1} = B_{i+1}/v_i,  B_i = u_i R_i,  a_i = (1 - u_i) R_i;
    b_1 = B_1.
This file is independent of htk except where stated (htk is used only as a second brute-force layer).
"""
from fractions import Fraction as Fr
import itertools
import math
import random

import sympy as sp


# ------------------------------------------------------------------ coordinates

def t_of_r(n, u, v):
    """Schwinger times [a_1, b_1, ..., a_n, b_n, b_{n+1}] from u[1..n], v[1..n] (dicts or lists
    indexed from 1); works for Fractions and sympy expressions."""
    B = {n: u[n]}
    b = {n + 1: v[n]}
    a = {n: 1 - u[n] - v[n]}
    R = {}
    for i in range(n - 1, 0, -1):
        R[i] = B[i + 1] / (1 - v[i])
        b[i + 1] = B[i + 1] / v[i]
        B[i] = u[i] * R[i]
        a[i] = (1 - u[i]) * R[i]
    b[1] = B[1]
    t = []
    for i in range(1, n + 1):
        t += [a[i], b[i]]
    t.append(b[n + 1])
    return t


def r_of_t(n, t):
    """Forward map: the ratios u_k, v_k of S1 (eq-fan-local-ratios) from Schwinger times t (any
    scale).  Returns (u, v) dicts."""
    a = {i: t[2 * (i - 1)] for i in range(1, n + 1)}
    b = {i: t[2 * (i - 1) + 1] for i in range(1, n + 1)}
    b[n + 1] = t[2 * n]
    u, v = {}, {}
    Bred = b[1]
    for i in range(1, n):
        R = a[i] + Bred
        u[i] = Bred / R
        v[i] = R / (R + b[i + 1])
        Bred = 1 / (1 / R + 1 / b[i + 1])
    S = a[n] + Bred + b[n + 1]
    u[n] = Bred / S
    v[n] = b[n + 1] / S
    return u, v


def random_r(n, rng):
    u, v = {}, {}
    for k in range(1, n):
        u[k] = Fr(rng.randint(1, 40), 41)
        v[k] = Fr(rng.randint(1, 40), 41)
    while True:
        x, y = Fr(rng.randint(1, 40), 41), Fr(rng.randint(1, 40), 41)
        if x + y < 1:
            break
    u[n], v[n] = x, y
    return u, v


# ------------------------------------------------------------------ exact linear algebra

def fdet(M):
    M = [[Fr(x) for x in r] for r in M]
    n = len(M)
    d = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            if f:
                for c in range(i, n):
                    M[r][c] -= f * M[i][c]
    return d


def fsolve(A, b):
    n = len(A)
    M = [[Fr(x) for x in A[i]] + [Fr(b[i])] for i in range(n)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [x / piv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]


def finv(A):
    n = len(A)
    cols = [fsolve(A, [Fr(int(i == j)) for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


# ------------------------------------------------------------------ exterior algebra (own code)
# forms: dict {tuple of increasing indices: coefficient}

def wedge(a, b):
    out = {}
    for ka, va in a.items():
        sa = set(ka)
        for kb, vb in b.items():
            if sa.intersection(kb):
                continue
            inv = sum(1 for x in ka for y in kb if x > y)
            key = tuple(sorted(ka + kb))
            val = va * vb * (-1 if inv % 2 else 1)
            out[key] = out.get(key, 0) + val
    return {k: v for k, v in out.items() if v != 0}


def fadd(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def pfaff(A, idx):
    if not idx:
        return {(): Fr(1)}
    out = {}
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        term = wedge(A[idx[0]][idx[j]], pfaff(A, rest))
        out = fadd(out, term, -1 if (j - 1) % 2 else 1)
    return out


# ------------------------------------------------------------------ the fan graph (own encoding)

def fan_edges(n):
    """Vertices 0..n (path, S1 vertex i = i-1), root n+1.  Edge order of htk.fan(n)."""
    E = []
    for i in range(n):
        E += [(i, i + 1), (i, n + 1)]
    E.append((n, n + 1))
    return E


def cycles_triangles(n):
    """Cycle matrix with rows = triangles c_i = a_i + b_{i+1} - b_i (i = 1..n), first row sign
    chosen such that det[B^T; C] > 0 (cf. htk CONVENTIONS Section 1; D = 2 makes K independent of
    the sign anyway)."""
    N = 2 * n + 1
    C = []
    for i in range(1, n + 1):
        row = [0] * N
        row[2 * (i - 1)] = 1                       # a_i : i -> i+1
        row[2 * i + 1 if i < n else 2 * n] = 1     # b_{i+1}: i+1 -> O
        row[2 * (i - 1) + 1] = -1                  # b_i  : back O -> i
        C.append(row)
    Bm = incidence(n)
    M = [[Bm[e][v] for e in range(N)] for v in range(n + 1)] + C
    if fdet(M) < 0:
        C[0] = [-x for x in C[0]]
    return C


def incidence(n):
    E = fan_edges(n)
    g = n + 1
    B = [[0] * g for _ in E]
    for e, (x, y) in enumerate(E):
        if x < g:
            B[e][x] += 1
        if y < g:
            B[e][y] -= 1
    return B


def kirchhoff(n, t, lam):
    """Currents y_e (scalar sources lam[0..n] at the path vertices, H = 1): Ohm t_e y_e =
    phi_{e(0)} - phi_{e(1)}, current law with sources lam.  Also returns the potentials."""
    E = fan_edges(n)
    g = n + 1
    L = [[Fr(0)] * g for _ in range(g)]
    for (x, y), te in zip(E, t):
        c = 1 / Fr(te)
        for p in (x, y):
            if p < g:
                L[p][p] += c
        if x < g and y < g:
            L[x][y] -= c
            L[y][x] -= c
    phi = fsolve(L, lam) + [Fr(0)]
    y = [(phi[x] - phi[yy]) / Fr(te) for (x, yy), te in zip(E, t)]
    return y, phi, L


def density_brute(n, t, lams, T):
    """D(t) with  U^{-H} rho_n(M)^T ^ prod_{i=1..H} prod_r psi^i_r = D(t) iota_E prod_e dt_e,
    lams = list of H source vectors (each over the n+1 path vertices).  Returns (D_rational_part,
    radicand) with D = D_rational_part * radicand^{-1/2} (T = 1: radicand = det M), or
    (D, 1) for T = 0.  Own implementation (no htk)."""
    N = 2 * n + 1
    t = [Fr(x) for x in t]
    C = cycles_triangles(n)
    M = [[sum(Fr(C[r][e] * C[s][e]) * t[e] for e in range(N)) for s in range(n)] for r in range(n)]
    U = fdet(M)
    forms = []
    if T == 1:
        Minv = finv(M)
        dM = [[{(e,): Fr(C[r][e] * C[s][e])} for e in range(N) if C[r][e] * C[s][e] != 0] for r in range(n) for s in range(n)]
        dM = [[{(e,): Fr(C[r][e] * C[s][e]) for e in range(N) if C[r][e] * C[s][e] != 0} for s in range(n)] for r in range(n)]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                acc = {}
                for a_ in range(n):
                    for b_ in range(n):
                        if Minv[a_][b_] != 0:
                            acc = fadd(acc, wedge(dM[r][a_], dM[s][b_]), -Minv[a_][b_] / 2)
                A[r][s] = acc
                A[s][r] = {k: -x for k, x in acc.items()}
        forms.append(pfaff(A, list(range(n))))
    for lam in lams:
        y, _, _ = kirchhoff(n, t, lam)
        for r in range(n):
            forms.append({(e,): C[r][e] * y[e] for e in range(N) if C[r][e] != 0})
    w = {(): Fr(1)}
    for f in forms:
        w = wedge(w, f)
    H = len(lams)
    # basic coefficient: D = (-1)^{e0} w[all but e0] / t_{e0}, e0 = N-1; check one more e0
    vals = []
    for e0 in (N - 1, 0):
        key = tuple(e for e in range(N) if e != e0)
        c = w.get(key, Fr(0))
        vals.append((-c if e0 % 2 else c) / t[e0])
    assert vals[0] == vals[1], vals
    D = vals[0] / U ** H
    return (D, U) if T == 1 else (D, Fr(1))


def jacobian_r(n, u, v):
    """det[t; dt/du_1; dt/dv_1; ...; dt/du_n; dt/dv_n] (rows), exact, in the coordinate order
    r_1, r_2, ..., r_{2n} = u_1, v_1, ..., u_n, v_n."""
    syms_u = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    syms_v = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    texpr = t_of_r(n, syms_u, syms_v)
    order = []
    for k in range(1, n + 1):
        order += [syms_u[k], syms_v[k]]
    sub = {syms_u[k]: sp.Rational(u[k].numerator, u[k].denominator) for k in syms_u}
    sub.update({syms_v[k]: sp.Rational(v[k].numerator, v[k].denominator) for k in syms_v})
    rows = [[te.subs(sub) for te in texpr]]
    for x in order:
        rows.append([sp.diff(te, x).subs(sub) for te in texpr])
    d = sp.Matrix(rows).det()
    return Fr(int(sp.Rational(d).p), int(sp.Rational(d).q))


# ------------------------------------------------------------------ the derived (1,1) density

def Pi_prefactor(n, u, v):
    """Pi(r) = v_n^{-1} prod_{k<n} v_k^{1/2} (1-v_k)^{k/2-1} prod_{k=1}^n u_k^{-(k+1)/2}
    returned as (rational, radicand) with Pi = rational * sqrt(radicand)."""
    rat, rad = Fr(1), Fr(1)

    def mul(x, e2):           # x^(e2/2)
        nonlocal rat, rad
        q, rr = divmod(e2, 2)
        rat *= x ** q
        if rr:
            rad *= x
    mul(v[n], -2)
    for k in range(1, n):
        mul(v[k], 1)
        mul(1 - v[k], k - 2)
    for k in range(1, n + 1):
        mul(u[k], -(k + 1))
    return rat, rad


def LM(n, u, v, lam):
    """L_i (i=1..n) and M_j (j=1..n) of the derivation (lam[0..n] = sources at S1 vertices 1..n+1)."""
    P = {1: Fr(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * u[k - 1]
    R = {n: Fr(1)}
    for k in range(n - 1, 0, -1):
        R[k] = R[k + 1] * u[k + 1] / (1 - v[k])
    c = {k: u[k] * (1 - u[k]) * R[k] for k in range(1, n + 1)}
    L = {}
    acc = Fr(0)
    for i in range(1, n + 1):
        acc += Fr(lam[i - 1]) / P[i]
        L[i] = acc
    M = {}
    for j in range(1, n + 1):
        s = sum((P[k] ** 2 * c[k] * (L[k] - L[j]) for k in range(j + 1, n + 1)), Fr(0))
        M[j] = s + P[n] * u[n] * v[n] * Fr(lam[n])
    return L, M


def matchings(seq):
    if not seq:
        yield []
        return
    a = seq[0]
    for k in range(1, len(seq)):
        rest = seq[1:k] + seq[k + 1:]
        for m in matchings(rest):
            yield [(a, seq[k])] + m


def hafnian_LM(n, L, M):
    tot = 0
    for m in matchings(list(range(1, n + 1))):
        p = 1
        for (i, j) in m:
            p *= L[i] * M[j]
        tot += p
    return tot


def density_11(n, u, v, lam):
    """g(r) = 2^{-n/2} Pi(r) sum_matchings prod L_i M_j  as (rational, radicand)."""
    rat, rad = Pi_prefactor(n, u, v)
    L, M = LM(n, u, v, lam)
    return rat * hafnian_LM(n, L, M) / 2 ** (n // 2), rad


def density_11_endpoint(n, u, v, l1, ln1):
    """Endpoint sources: g = (l1 ln1)^{n/2} (n-1)!!/2^{n/2} v_n^{n/2-1}
    prod_{k<n} v_k^{1/2}(1-v_k)^{k/2-1} prod_{m=1}^n u_m^{(n-1-m)/2}."""
    rat, rad = Fr(1), Fr(1)

    def mul(x, e2):
        nonlocal rat, rad
        q, rr = divmod(e2, 2)
        rat *= x ** q
        if rr:
            rad *= x
    mul(v[n], n - 2)
    for k in range(1, n):
        mul(v[k], 1)
        mul(1 - v[k], k - 2)
    for m in range(1, n + 1):
        mul(u[m], n - 1 - m)
    dfact = 1
    for k in range(n - 1, 0, -2):
        dfact *= k
    return rat * Fr(l1 * ln1) ** (n // 2) * dfact / 2 ** (n // 2), rad


def density_20_endpoint(n, u, v, d):
    """S1 eq-1-41 times (-1)^{n(n+1)/2}: K^{(2,0)} density, d = det(lam_1, lam_{n+1})."""
    g = Fr(d) ** n * v[n] ** (n - 1)
    for j in range(1, n):
        g *= u[j] ** (n - j) * v[j] * (1 - v[j]) ** (j - 1)
    return g * (-1) ** (n * (n + 1) // 2)


def density_20_potential(n, u, v, lam1, lam2):
    """(2,0) density from the potential form: (-1)^{n(n+1)/2} J(r) prod_i det(Phi_i, Phi_{i+1}),
    J = prod_{i<n} v_i R_i/(a_i B_i B_{i+1}) / (a_n B_n b_{n+1}),  Phi = Kirchhoff potentials."""
    t = t_of_r(n, u, v)
    _, phi1, _ = kirchhoff(n, t, lam1)
    _, phi2, _ = kirchhoff(n, t, lam2)
    a = {i: t[2 * (i - 1)] for i in range(1, n + 1)}
    B = {n: u[n]}
    R = {}
    for i in range(n - 1, 0, -1):
        R[i] = B[i + 1] / (1 - v[i])
        B[i] = u[i] * R[i]
    J = Fr(1)
    for i in range(1, n):
        J *= v[i] * R[i] / (a[i] * B[i] * B[i + 1])
    J /= a[n] * B[n] * v[n]
    prod = Fr(1)
    for i in range(1, n + 1):      # S1 vertex i -> index i-1
        prod *= phi1[i - 1] * phi2[i] - phi2[i - 1] * phi1[i]
    return (-1) ** (n * (n + 1) // 2) * J * prod
