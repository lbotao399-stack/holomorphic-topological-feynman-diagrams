"""Independent verification core (written from the definitions; does not import the explorer's code).

Graph: nv vertices, edges[e] = (e0, e1) oriented; KCL  sum_{e0=v} y_e - sum_{e1=v} y_e = lam_v (all v, sum lam = 0);
KVL  sum_e Cb_re t_e y_e = 0 for a BFS fundamental cycle basis Cb (not C3).
Forms are represented by their values on vectors: a 1-form is an N-vector a (a(v) = a.v); a 2-form w(u, v).
Omega = U^{-H} rho_2^T ^ prod_{i=1..H} prod_{r=1,2} psi^i_r,  psi^i_r = sum_e C_re y^i_e dt_e,
rho_2 = U^{-1/2} A12, A12 = -1/2 sum_{ab} Minv_ab dM_{1a} ^ dM_{b2},  M = C t C^T.
Density: Omega(v_1..v_{N-1}) = D(t) det[t; v_1; ..; v_{N-1}]  (Omega basic, iota_E vol(v) = det[t; v]).
"""
from fractions import Fraction as Fr
import itertools
import math
import random


# ----------------------------------------------------------------------------- exact linear algebra
def det(M):
    n = len(M)
    if n == 0:
        return Fr(1)
    A = [[Fr(x) for x in row] for row in M]
    s = Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            s = -s
        pv = A[c][c]
        s *= pv
        for r in range(c + 1, n):
            f = A[r][c] / pv
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return s


def solve(A, cols):
    """A X = B with B given as a list of columns; returns list of solution columns (exact)."""
    n = len(A)
    m = len(cols)
    W = [[Fr(x) for x in A[i]] + [Fr(b[i]) for b in cols] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if W[r][c] != 0)
        W[c], W[piv] = W[piv], W[c]
        pv = W[c][c]
        W[c] = [x / pv for x in W[c]]
        for r in range(n):
            if r != c and W[r][c] != 0:
                f = W[r][c]
                W[r] = [x - f * y for x, y in zip(W[r], W[c])]
    return [[W[i][n + j] for i in range(n)] for j in range(m)]


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def sg(k):
    return -1 if k % 2 else 1


# ----------------------------------------------------------------------------- theta graphs
class Theta:
    """Theta graph with explicit path data.
    paths[P] = [(edge, s)] from N to N' (s = +1 if the edge points along N -> N'),
    interior[P] = [vertices] from N."""

    def __init__(self, nv, edges, N, Np, paths, interior):
        self.nv, self.edges, self.N, self.Np = nv, list(edges), N, Np
        self.paths, self.interior = paths, interior
        self.n = {P: len(paths[P]) for P in "ABC"}
        self.pqr = tuple(self.n[P] for P in "ABC")
        self.std = [e for P in "ABC" for e, _ in paths[P]]
        self.order_sign = perm_sign(self.std)

    @property
    def ne(self):
        return len(self.edges)

    def C3(self):
        C = [[0] * self.ne, [0] * self.ne]
        for e, s in self.paths["A"]:
            C[0][e] += s
        for e, s in self.paths["B"]:
            C[0][e] -= s
            C[1][e] += s
        for e, s in self.paths["C"]:
            C[1][e] -= s
        return C


def theta_std(p, q, r):
    """N = 0, N' = nv-1, interior numbered along A, B, C, edge order A, B, C, all edges N -> N'
    (same labelling as htk.graph.theta and S3 forms_check.theta)."""
    nv = p + q + r - 1
    Np = nv - 1
    edges, paths, interior, nxt = [], {}, {}, 1
    for P, n_ in zip("ABC", (p, q, r)):
        prev, ints, pe = 0, [], []
        for j in range(n_):
            if j == n_ - 1:
                cur = Np
            else:
                cur = nxt
                nxt += 1
                ints.append(cur)
            pe.append((len(edges), 1))
            edges.append((prev, cur))
            prev = cur
        paths[P], interior[P] = pe, ints
    return Theta(nv, edges, 0, Np, paths, interior)


def theta_from_graph(nv, edges, N, Np, order=None):
    """Trace the three paths N -> N' in an arbitrary labelled graph; order = permutation assigning
    the traced paths (sorted by (length, smallest edge)) to A, B, C."""
    inc = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        inc[a].append(e)
        inc[b].append(e)
    found = []
    for e0 in inc[N]:
        pe, ints, cur, e = [], [], N, e0
        while True:
            a, b = edges[e]
            s, nxt = (1, b) if a == cur else (-1, a)
            pe.append((e, s))
            if nxt == Np:
                break
            ints.append(nxt)
            cur = nxt
            e = next(f for f in inc[cur] if f != e)
        found.append((pe, ints))
    found.sort(key=lambda x: (len(x[0]), min(e for e, _ in x[0])))
    if order is not None:
        found = [found[i] for i in order]
    return Theta(nv, edges, N, Np, {P: found[i][0] for i, P in enumerate("ABC")},
                 {P: found[i][1] for i, P in enumerate("ABC")})


def transform(th, vperm=None, eperm=None, reverse=()):
    """New Theta: vertex v -> vperm[v]; new edge j = old edge eperm[j]; old edges in `reverse` flipped."""
    nv, ne = th.nv, th.ne
    vperm = list(range(nv)) if vperm is None else vperm
    eperm = list(range(ne)) if eperm is None else eperm
    inv = [0] * ne
    for j, e in enumerate(eperm):
        inv[e] = j
    rev = set(reverse)
    edges = []
    for j in range(ne):
        a, b = th.edges[eperm[j]]
        if eperm[j] in rev:
            a, b = b, a
        edges.append((vperm[a], vperm[b]))
    paths = {P: [(inv[e], -s if e in rev else s) for e, s in th.paths[P]] for P in "ABC"}
    interior = {P: [vperm[u] for u in th.interior[P]] for P in "ABC"}
    return Theta(nv, edges, vperm[th.N], vperm[th.Np], paths, interior)


def theta_list(D):
    return [(p, q, r) for p in range(1, D + 1) for q in range(p, D + 1) for r in range(q, D + 1)
            if p + q + r == 2 * D + 1]


# ----------------------------------------------------------------------------- Kirchhoff (exact)
def incidence_full(nv, edges):
    """Bf[e][v] = +1 (v = e0), -1 (v = e1)."""
    Bf = [[0] * nv for _ in edges]
    for e, (a, b) in enumerate(edges):
        Bf[e][a] += 1
        Bf[e][b] -= 1
    return Bf


def bfs_cycles(nv, edges, root):
    """Fundamental cycles of a BFS tree from `root` (rows over edges)."""
    adj = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        adj[a].append((b, e, 1))
        adj[b].append((a, e, -1))
    par = {root: None}
    order = [root]
    tree = set()
    for v in order:
        for u, e, s in adj[v]:
            if u not in par:
                par[u] = (v, e, s)
                tree.add(e)
                order.append(u)

    def path(v):  # signed edge vector of the tree path root -> v
        vec = [0] * len(edges)
        while par[v] is not None:
            u, e, s = par[v]
            vec[e] += s
            v = u
        return vec
    C = []
    for e, (a, b) in enumerate(edges):
        if e in tree:
            continue
        pa, pb = path(a), path(b)
        row = [x - y for x, y in zip(pa, pb)]
        row[e] += 1
        C.append(row)
    return C


def is_cycle_basis_row(nv, edges, row):
    Bf = incidence_full(nv, edges)
    return all(sum(row[e] * Bf[e][v] for e in range(len(edges))) == 0 for v in range(nv))


def kirchhoff(nv, edges, t, lam, root=None):
    """Exact currents y[e][i] from KCL (all vertices but `root`) + KVL (BFS cycles).
    lam: full list over vertices of H-vectors with sum zero."""
    root = nv - 1 if root is None else root
    H = len(lam[0])
    ne = len(edges)
    Bf = incidence_full(nv, edges)
    Cb = bfs_cycles(nv, edges, root)
    A, rhs = [], [[] for _ in range(H)]
    for v in range(nv):
        if v == root:
            continue
        A.append([Bf[e][v] for e in range(ne)])
        for i in range(H):
            rhs[i].append(Fr(lam[v][i]))
    for row in Cb:
        A.append([row[e] * Fr(t[e]) for e in range(ne)])
        for i in range(H):
            rhs[i].append(Fr(0))
    if H == 0:
        return [[] for _ in range(ne)]
    sol = solve(A, rhs)
    y = [[sol[i][e] for i in range(H)] for e in range(ne)]
    # checks: KCL at every vertex, KVL
    for v in range(nv):
        for i in range(H):
            assert sum(Bf[e][v] * y[e][i] for e in range(ne)) == lam[v][i]
    for row in Cb:
        for i in range(H):
            assert sum(row[e] * Fr(t[e]) * y[e][i] for e in range(ne)) == 0
    return y


def kirchhoff_dt(nv, edges, t, y, root=None):
    """dy[e][i][j] = d y^i_e / d t_j by implicit differentiation of the KCL+KVL system."""
    root = nv - 1 if root is None else root
    ne = len(edges)
    H = len(y[0]) if y else 0
    Bf = incidence_full(nv, edges)
    Cb = bfs_cycles(nv, edges, root)
    A = [[Bf[e][v] for e in range(ne)] for v in range(nv) if v != root]
    nk = len(A)
    for row in Cb:
        A.append([row[e] * Fr(t[e]) for e in range(ne)])
    cols = []
    for j in range(ne):
        for i in range(H):
            b = [Fr(0)] * nk + [-Fr(row[j]) * y[j][i] for row in Cb]
            cols.append(b)
    sol = solve(A, cols)
    dy = [[[None] * ne for _ in range(H)] for _ in range(ne)]
    k = 0
    for j in range(ne):
        for i in range(H):
            for e in range(ne):
                dy[e][i][j] = sol[k][e]
            k += 1
    return dy


# ----------------------------------------------------------------------------- forms (exact)
def matvec(F, v):
    return sum(Fr(a) * Fr(b) for a, b in zip(F, v))


def omega_on(C, t, y, H, T, vecs):
    """Omega(vecs) for T = 0, and sqrt(U) * Omega(vecs) for T = 1 (rational); also returns U."""
    ne = len(t)
    t = [Fr(x) for x in t]
    n = len(C)
    M = [[sum(Fr(C[a][e] * C[b][e]) * t[e] for e in range(ne)) for b in range(n)] for a in range(n)]
    U = det(M)
    F = [[Fr(C[r][e]) * y[e][i] for e in range(ne)] for i in range(H) for r in range(n)]
    FV = [[matvec(Fk, v) for v in vecs] for Fk in F]
    if T == 0:
        return det(FV) / U ** H, U
    assert n == 2
    Minv = [[M[1][1] / U, -M[0][1] / U], [-M[1][0] / U, M[0][0] / U]]
    dM = [[[Fr(C[a][e] * C[b][e]) for e in range(ne)] for b in range(n)] for a in range(n)]
    dMv = [[[matvec(dM[a][b], v) for v in vecs] for b in range(n)] for a in range(n)]

    def A12(k, l):
        s = Fr(0)
        for a in range(n):
            for b in range(n):
                s += Minv[a][b] * (dMv[0][a][k] * dMv[b][1][l] - dMv[0][a][l] * dMv[b][1][k])
        return -s / 2
    m = len(vecs)
    tot = Fr(0)
    for k in range(m):
        for l in range(k + 1, m):
            w = A12(k, l)
            if w == 0:
                continue
            rest = [c for c in range(m) if c not in (k, l)]
            d = det([[row[c] for c in rest] for row in FV]) if F else Fr(1)
            tot += sg(k + l + 1 + 2) * w * d  # 1-based (k+1)+(l+1)+1
    return tot / U ** H, U


def density_bf(th, t, lam, H, T, C=None, nvec=2, rng=None, root=None):
    """Brute-force D(t) (T = 0) or sqrt(U) D(t) (T = 1) relative to iota_E(dt_0..dt_{N-1}) in the graph's
    edge order, from Omega evaluated on random rational vectors (two independent vector sets must agree)."""
    rng = rng or random.Random(1)
    C = th.C3() if C is None else C
    ne = th.ne
    y = kirchhoff(th.nv, th.edges, t, lam, root=root) if H else [[] for _ in range(ne)]
    vals = []
    for _ in range(nvec):
        while True:
            vecs = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(ne)] for _ in range(ne - 1)]
            dd = det([list(map(Fr, t))] + vecs)
            if dd != 0:
                break
        om, U = omega_on(C, t, y, H, T, vecs)
        vals.append(om / dd)
    assert all(v == vals[0] for v in vals), "Omega not proportional to iota_E vol"
    return vals[0], U, y


# ----------------------------------------------------------------------------- closed forms (claims)
def col(top, bot, H):
    return [Fr(x) for x in top] + [Fr(x) for x in bot]


def D_theta(th, lam, H):
    if H == 0:
        return Fr(1)
    z = [0] * H
    cols = [col(lam[u], z, H) for u in th.interior["A"]] + [col(lam[u], lam[u], H) for u in th.interior["B"]] + \
           [col(z, lam[u], H) for u in th.interior["C"]]
    assert len(cols) == 2 * H
    return det([[c[i] for c in cols] for i in range(2 * H)])


def Qform(th, lam, H, v):
    z = [0] * H
    cols = [col(v, z, H), col(z, v, H)] + [col(lam[u], z, H) for u in th.interior["A"]] + \
           [col(lam[u], lam[u], H) for u in th.interior["B"]] + [col(z, lam[u], H) for u in th.interior["C"]]
    assert len(cols) == 2 * H
    return det([[c[i] for c in cols] for i in range(2 * H)])


def Rsum(th, t):
    return {P: sum(Fr(t[e]) for e, _ in th.paths[P]) for P in "ABC"}


def sigma_t(th, t):
    R = Rsum(th, t)
    out = {}
    for P in "ABC":
        es = [e for e, _ in th.paths[P]]
        out[P] = [sum(Fr(t[e]) for e in es[k:]) / R[P] for k in range(1, len(es))]
    return out


def L_t(th, t, lam, H):
    sig = sigma_t(th, t)
    L = [Fr(x) for x in lam[th.N]]
    for P in "ABC":
        for s, u in zip(sig[P], th.interior[P]):
            L = [a + s * Fr(b) for a, b in zip(L, lam[u])]
    return L


def density_closed(th, t, lam, H, T):
    """claimed density for C3 in the graph's edge order: D (T = 0), sqrt(U) D (T = 1)."""
    R = Rsum(th, t)
    U = R["A"] * R["B"] + R["B"] * R["C"] + R["C"] * R["A"]
    if T == 1:
        return th.order_sign * sg(H * (H + 1) // 2) * D_theta(th, lam, H) / (2 * U ** (H + 1)), U
    L = L_t(th, t, lam, H)
    return th.order_sign * sg(H * (H + 1) // 2 + 1) * R["A"] * R["B"] * R["C"] * Qform(th, lam, H, L) / U ** (H + 2), U


def full_sources(nv, root, vecs):
    """vecs: dict v -> H-vector for v != root; root gets minus the sum."""
    H = len(next(iter(vecs.values())))
    lam = [[Fr(0)] * H for _ in range(nv)]
    for v, x in vecs.items():
        if v != root:
            lam[v] = [Fr(a) for a in x]
    lam[root] = [-sum(lam[v][i] for v in range(nv) if v != root) for i in range(H)]
    return lam


def rand_sources(nv, root, H, rng, lo=-9, hi=9, den=4):
    return full_sources(nv, root, {v: [Fr(rng.randint(lo, hi), rng.randint(1, den)) for _ in range(H)]
                                   for v in range(nv) if v != root}) if H else [[] for _ in range(nv)]


def rand_t(ne, rng):
    return [Fr(rng.randint(1, 60), rng.randint(1, 13)) for _ in range(ne)]


# ----------------------------------------------------------------------------- exact z = 0 values (own integration)
def ghalf(m):
    """Gamma(m + 1/2)/sqrt(pi) for integer m (any sign), from Gamma(x + 1) = x Gamma(x)."""
    g = Fr(1)
    if m >= 0:
        for j in range(m):
            g *= Fr(2 * j + 1, 2)
    else:
        for j in range(m, 0):
            g /= Fr(2 * j + 1, 2)
    return g


def sigma_moments(n):
    """Exact first and second moments of sigma_k = sum_{j>k} x_j, x uniform on Delta^{n-1}
    (normalized measure), from E[x_j] = 1/n, E[x_j x_l] = (1 + delta_jl)/(n(n+1))."""
    Ex = Fr(1, n)
    Exx = lambda j, l: Fr(1 + (j == l), n * (n + 1))
    m1 = [sum(Ex for j in range(k + 1, n + 1)) for k in range(1, n)]
    m2 = [[sum(Exx(j, l) for j in range(k + 1, n + 1) for l in range(kk + 1, n + 1)) for kk in range(1, n)]
          for k in range(1, n)]
    return m1, m2
