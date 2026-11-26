"""Independent verification library for explore/fans-ht (does not import fanlib or htk).

Fan (S1 Theorem thm-fan-kernel): path vertices 1..n+1, insertion O = n+2 (root), edges
a_i = (i, i+1), b_i = (i, O); edge order a_1, b_1, ..., a_n, b_n, b_{n+1} (index 2i-2, 2i-1, 2n).

Brute force: currents by mesh analysis (tree flow + loop currents), cycle basis = fundamental
cycles of a spanning tree chosen here (for the fan: path tree a_1..a_n, b_{n+1}), normalized by
det[B^T; C] > 0; forms = dicts {bitmask: coefficient}; Pfaffian = signed sum over perfect matchings.
"""
from fractions import Fraction as Fr
import itertools
import random

import sympy as sp


# ---------------------------------------------------------------- graphs

class Graph:
    def __init__(self, verts, root, edges):
        self.verts = list(verts)          # non-root vertices, in order
        self.root = root
        self.edges = list(edges)
        self.N = len(edges)
        self.n = self.N - len(self.verts)  # loops = N - (|V|-1)
        self.vpos = {v: k for k, v in enumerate(self.verts)}

    def incidence(self):
        B = [[0] * len(self.verts) for _ in self.edges]
        for e, (x, y) in enumerate(self.edges):
            if x in self.vpos:
                B[e][self.vpos[x]] += 1
            if y in self.vpos:
                B[e][self.vpos[y]] -= 1
        return B


def fan_graph(n):
    O = n + 2
    E = []
    for i in range(1, n + 1):
        E += [(i, i + 1), (i, O)]
    E.append((n + 1, O))
    return Graph(range(1, n + 2), O, E)


def ia(i):
    return 2 * i - 2


def ib(i, n):
    return 2 * i - 1 if i <= n else 2 * n


# ---------------------------------------------------------------- exact linear algebra

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
        inv = 1 / M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] * inv
            if f:
                for c in range(i, n):
                    M[r][c] -= f * M[i][c]
    return d


def fsolve_multi(A, Bcols):
    """Solve A X = B for X (B given as list of columns)."""
    n = len(A)
    m = len(Bcols)
    M = [[Fr(A[i][j]) for j in range(n)] + [Fr(Bcols[k][i]) for k in range(m)] for i in range(n)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [x / piv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [[M[i][n + k] for i in range(n)] for k in range(m)]


def finv(A):
    n = len(A)
    cols = fsolve_multi(A, [[int(i == j) for i in range(n)] for j in range(n)])
    return [[cols[j][i] for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------- trees, cycles, currents

def spanning_tree_edges(G, prefer):
    """Spanning tree containing as many edges of `prefer` (list of edge indices, in order) as
    possible (greedy union-find over prefer, then the remaining edges)."""
    allv = G.verts + [G.root]
    par = {v: v for v in allv}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    T = []
    order = list(prefer) + [e for e in range(G.N) if e not in prefer]
    for e in order:
        x, y = G.edges[e]
        rx, ry = find(x), find(y)
        if rx != ry:
            par[rx] = ry
            T.append(e)
    assert len(T) == len(G.verts)
    return sorted(T)


def tree_path(G, T, p, q):
    """Edges (with signs +1 along, -1 against) of the tree path from p to q."""
    adj = {}
    for e in T:
        x, y = G.edges[e]
        adj.setdefault(x, []).append((y, e, +1))
        adj.setdefault(y, []).append((x, e, -1))
    prev = {p: None}
    stack = [p]
    while stack:
        x = stack.pop()
        for (y, e, s) in adj.get(x, []):
            if y not in prev:
                prev[y] = (x, e, s)
                stack.append(y)
    path = []
    x = q
    while prev[x] is not None:
        px, e, s = prev[x]
        path.append((e, s))
        x = px
    return list(reversed(path))


def fundamental_cycles(G, T):
    cot = [e for e in range(G.N) if e not in T]
    C = []
    for c in cot:
        row = [0] * G.N
        row[c] = 1
        p, q = G.edges[c]
        for (e, s) in tree_path(G, T, q, p):     # back from q to p through the tree
            row[e] += s
        C.append(row)
    return C


def normalize_cycles(G, C):
    B = G.incidence()
    Mt = [[B[e][v] for e in range(G.N)] for v in range(len(G.verts))] + [list(r) for r in C]
    d = fdet(Mt)
    assert d != 0
    if d < 0:
        C = [[-x for x in C[0]]] + [list(r) for r in C[1:]]
    return C


def check_cycles(G, C):
    B = G.incidence()
    for r in C:
        for v in range(len(G.verts)):
            assert sum(r[e] * B[e][v] for e in range(G.N)) == 0


def currents_mesh(G, t, lam, T=None):
    """Kirchhoff currents by mesh analysis: y = y0 + C_T^T l, y0 a tree flow with sources lam
    (lam[k] at vertex G.verts[k]), (C_T t C_T^T) l = -C_T t y0."""
    if T is None:
        T = spanning_tree_edges(G, [])
    B = G.incidence()
    # tree flow: solve B_T^T y_T = lam
    nv = len(G.verts)
    BT = [[B[e][v] for e in T] for v in range(nv)]
    yT = fsolve_multi(BT, [lam])[0]
    y0 = [Fr(0)] * G.N
    for k, e in enumerate(T):
        y0[e] = yT[k]
    CT = fundamental_cycles(G, T)
    n = len(CT)
    t = [Fr(x) for x in t]
    MT = [[sum(CT[r][e] * CT[s][e] * t[e] for e in range(G.N)) for s in range(n)] for r in range(n)]
    rhs = [-sum(CT[r][e] * t[e] * y0[e] for e in range(G.N)) for r in range(n)]
    l = fsolve_multi(MT, [rhs])[0]
    y = [y0[e] + sum(CT[r][e] * l[r] for r in range(n)) for e in range(G.N)]
    # checks: current law and voltage law
    for v in range(nv):
        assert sum(B[e][v] * y[e] for e in range(G.N)) == lam[v]
    for r in range(n):
        assert sum(CT[r][e] * t[e] * y[e] for e in range(G.N)) == 0
    return y


def laplacian(G, t):
    B = G.incidence()
    nv = len(G.verts)
    return [[sum(Fr(B[e][u] * B[e][v]) / Fr(t[e]) for e in range(G.N)) for v in range(nv)] for u in range(nv)]


# ---------------------------------------------------------------- exterior algebra

def _wsign(ma, mb):
    s = 0
    m = mb
    while m:
        low = m & -m
        y = low.bit_length() - 1
        s += bin(ma >> (y + 1)).count("1")
        m ^= low
    return -1 if s & 1 else 1


def wedge(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            if ma & mb:
                continue
            k = ma | mb
            v = ca * cb if _wsign(ma, mb) > 0 else -(ca * cb)
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def fadd(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def one_form(coeffs):
    return {1 << e: c for e, c in enumerate(coeffs) if c != 0}


def perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def matchings(seq):
    if not seq:
        yield []
        return
    a = seq[0]
    for k in range(1, len(seq)):
        rest = seq[1:k] + seq[k + 1:]
        for m in matchings(rest):
            yield [(a, seq[k])] + m


def signed_matchings(n):
    out = []
    for m in matchings(list(range(n))):
        flat = [x for p in m for x in p]
        out.append((m, perm_sign(flat)))
    return out


def pfaffian_forms(A, n):
    """Pf of an antisymmetric n x n matrix of 2-forms: sum over matchings sgn * prod A[i][j]."""
    tot = {}
    for m, s in signed_matchings(n):
        w = {0: Fr(1)}
        for (i, j) in m:
            w = wedge(w, A[i][j])
        tot = fadd(tot, w, s)
    return tot


def top_coefficient(form, N, e0):
    """(-1)^{e0} [form]_{all dt except dt_{e0}}  (to be divided by t_{e0})."""
    key = ((1 << N) - 1) ^ (1 << e0)
    c = form.get(key, Fr(0))
    return -c if e0 % 2 else c


# ---------------------------------------------------------------- brute-force density

def period(G, C, t):
    n = len(C)
    return [[sum(C[r][e] * C[s][e] * Fr(t[e]) for e in range(G.N)) for s in range(n)] for r in range(n)]


def A_forms(G, C, t):
    """A_{rs} = -1/2 sum_ab dM_ra (M^-1)_ab ^ dM_bs as 2-forms, computed from the definition
    dM_ra = sum_e C_re C_ae dt_e."""
    n = len(C)
    M = period(G, C, t)
    Mi = finv(M)
    dM = [[one_form([C[r][e] * C[a][e] for e in range(G.N)]) for a in range(n)] for r in range(n)]
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Mi[a][b] != 0:
                        acc = fadd(acc, wedge(dM[r][a], dM[b][s]), -Mi[a][b] / 2)
            A[r][s] = acc
    for r in range(n):
        for s in range(n):
            chk = fadd(A[r][s], A[s][r])
            assert chk == {}, "A not antisymmetric"
    return A, M


def brute_density(G, C, t, lams, T, e0s=None):
    """D(t): U^{-H} rho_n(M)^T ^ prod_{i} prod_r psi^i_r = D iota_E prod dt.
    lams: list of H source lists (over G.verts).  Returns (Drat, radicand):
    D = Drat (T=0) or D = Drat * radicand^{-1/2} (T=1, radicand = U)."""
    N = G.N
    t = [Fr(x) for x in t]
    n = len(C)
    w = {0: Fr(1)}
    if T == 1:
        A, M = A_forms(G, C, t)
        w = pfaffian_forms(A, n)
    else:
        M = period(G, C, t)
    U = fdet(M)
    for lam in lams:
        y = currents_mesh(G, t, lam)
        for r in range(n):
            w = wedge(w, one_form([C[r][e] * y[e] for e in range(N)]))
    if e0s is None:
        e0s = (N - 1, 0)
    vals = [top_coefficient(w, N, e0) / t[e0] for e0 in e0s]
    assert all(v == vals[0] for v in vals), vals
    H = len(lams)
    D = vals[0] / U ** H
    return D, (U if T == 1 else Fr(1))


def fan_cycles(n, kind="pathtree"):
    G = fan_graph(n)
    if kind == "pathtree":
        T = [ia(i) for i in range(1, n + 1)] + [ib(n + 1, n)]
        C = fundamental_cycles(G, sorted(T))
    elif kind == "triangles":
        C = []
        for i in range(1, n + 1):
            row = [0] * G.N
            row[ia(i)] = 1
            row[ib(i + 1, n)] = 1
            row[ib(i, n)] = -1
            C.append(row)
    check_cycles(G, C)
    return G, normalize_cycles(G, C)


# ---------------------------------------------------------------- the S1 chart (own derivation)

def t_of_r(n, u, v):
    """Inverse of S1 eq-fan-local-ratios with the final triangle (eq-1-25):
    (a_n, B_n, b_{n+1}) = (1-u_n-v_n, u_n, v_n); for i = n-1..1, from
    v_i = R_i/(R_i + b_{i+1}), B_{i+1} = R_i b_{i+1}/(R_i + b_{i+1}):  R_i = B_{i+1}/(1-v_i),
    b_{i+1} = B_{i+1}/v_i; from u_i = B_i/(a_i + B_i), R_i = a_i + B_i: B_i = u_i R_i,
    a_i = (1-u_i) R_i; B_1 = b_1.  u, v indexed 1..n (dicts)."""
    a, b, B = {}, {}, {}
    B[n] = u[n]
    b[n + 1] = v[n]
    a[n] = 1 - u[n] - v[n]
    for i in range(n - 1, 0, -1):
        R = B[i + 1] / (1 - v[i])
        b[i + 1] = B[i + 1] / v[i]
        B[i] = u[i] * R
        a[i] = (1 - u[i]) * R
    b[1] = B[1]
    t = []
    for i in range(1, n + 1):
        t += [a[i], b[i]]
    t.append(b[n + 1])
    return t


def r_of_t(n, t):
    """S1 eq-fan-local-ratios: r_{2i-1} = B_i/(a_i + B_i), r_{2i} = (a_i + B_i)/(a_i + B_i + b_{i+1}),
    B_{i+1} = [1/(a_i + B_i) + 1/b_{i+1}]^{-1}, B_1 = b_1; final triangle (a_n, B_n, b_{n+1})
    normalized to sum 1: (r_{2n-1}, r_{2n}) = (B_n, b_{n+1})/S."""
    a = {i: t[ia(i)] for i in range(1, n + 1)}
    b = {i: t[ib(i, n)] for i in range(1, n + 2)}
    u, v = {}, {}
    Bc = b[1]
    for i in range(1, n):
        R = a[i] + Bc
        u[i] = Bc / R
        v[i] = R / (R + b[i + 1])
        Bc = R * b[i + 1] / (R + b[i + 1])
    S = a[n] + Bc + b[n + 1]
    u[n] = Bc / S
    v[n] = b[n + 1] / S
    return u, v


def rand_r(n, rng, den=97):
    u, v = {}, {}
    for k in range(1, n):
        u[k] = Fr(rng.randint(1, den - 1), den)
        v[k] = Fr(rng.randint(1, den - 1), den)
    while True:
        x, y = Fr(rng.randint(1, den - 1), den), Fr(rng.randint(1, den - 1), den)
        if x + y < 1:
            break
    u[n], v[n] = x, y
    return u, v


_JAC_CACHE = {}


def jac_symbolic(n):
    if n in _JAC_CACHE:
        return _JAC_CACHE[n]
    us = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    vs = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    te = t_of_r(n, us, vs)
    order = []
    for k in range(1, n + 1):
        order += [us[k], vs[k]]
    rows = [[sp.diff(x, s) for x in te] for s in order]
    _JAC_CACHE[n] = (us, vs, te, rows)
    return _JAC_CACHE[n]


def to_fr(x):
    x = sp.Rational(x)
    return Fr(int(x.p), int(x.q))


def jac_det(n, u, v):
    """det[t; dt/du_1; dt/dv_1; ...; dt/du_n; dt/dv_n] at a rational point (exact)."""
    us, vs, te, rows = jac_symbolic(n)
    sub = {us[k]: sp.Rational(u[k].numerator, u[k].denominator) for k in us}
    sub.update({vs[k]: sp.Rational(v[k].numerator, v[k].denominator) for k in vs})
    t = t_of_r(n, u, v)
    mat = [list(t)] + [[to_fr(x.subs(sub)) for x in row] for row in rows]
    return fdet(mat)


# ---------------------------------------------------------------- claimed formulas (re-implemented from the statements)

def half_power(pairs):
    """prod x^(e2/2) for pairs (x, e2) -> (rational, radicand)."""
    rat, rad = Fr(1), Fr(1)
    for x, e2 in pairs:
        q, rr = divmod(e2, 2)
        rat *= Fr(x) ** q
        if rr:
            rad *= x
    return rat, rad


def claim_P(n, u):
    P = {1: Fr(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * u[k - 1]
    return P


def claim_R(n, u, v):
    R = {n: Fr(1)}
    for i in range(n - 1, 0, -1):
        R[i] = R[i + 1] * u[i + 1] / (1 - v[i])
    return R


def claim_gh(n, u, v):
    P = claim_P(n, u)
    R = claim_R(n, u, v)
    c = {k: u[k] * (1 - u[k]) * R[k] for k in range(1, n + 1)}
    g = {w: 1 / P[w] for w in range(1, n + 1)}
    g[n + 1] = (1 - v[n]) / (P[n] * u[n])
    h = {w: sum((P[k] ** 2 * c[k] for k in range(w, n + 1)), Fr(0)) / P[w] for w in range(1, n + 1)}
    h[n + 1] = P[n] * u[n] * v[n]
    return g, h, P, R, c


def claim_LM(n, u, v, lam):
    """lam: dict 1..n+1."""
    g, h, P, R, c = claim_gh(n, u, v)
    L = {}
    acc = Fr(0)
    for i in range(1, n + 1):
        acc += Fr(lam[i]) / P[i]
        L[i] = acc
    M = {}
    for j in range(1, n + 1):
        M[j] = sum((P[k] ** 2 * c[k] * (L[k] - L[j]) for k in range(j + 1, n + 1)), Fr(0)) + P[n] * u[n] * v[n] * Fr(lam[n + 1])
    return L, M


def claim_Pi(n, u, v):
    pairs = [(v[n], -2)]
    for k in range(1, n):
        pairs += [(v[k], 1), (1 - v[k], k - 2)]
    for k in range(1, n + 1):
        pairs.append((u[k], -(k + 1)))
    return half_power(pairs)


def haf_LM(n, L, M):
    tot = Fr(0)
    for m in matchings(list(range(1, n + 1))):
        p = Fr(1)
        for (i, j) in m:
            p *= L[i] * M[j]
        tot += p
    return tot


def claim_density(n, u, v, lam):
    rat, rad = claim_Pi(n, u, v)
    L, M = claim_LM(n, u, v, lam)
    return rat * haf_LM(n, L, M) / 2 ** (n // 2), rad


def same_halfint(x, rx, y, ry):
    """x sqrt(rx) == y sqrt(ry) with rx, ry > 0."""
    if (x > 0) != (y > 0) or (x == 0) != (y == 0):
        return False
    return x * x * rx == y * y * ry
