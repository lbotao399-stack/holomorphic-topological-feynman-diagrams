"""Independent brute-force layer for the verification of ht-families (written from CONVENTIONS.md
sections 1-2 only; no code of the explorer or of htk is imported here).

Exact layer (Fractions):
  Net(nv, edges): root = nv-1, B_ev = delta_{v,e0} - delta_{v,e1} (v < root).
  DFS spanning tree from the root (adjacency scanned in reverse edge order), fundamental cycles,
  first row flipped so that det[B^T; C] > 0  (positive basis of CONVENTIONS 1).
  Kirchhoff: Lambda = B^T t^{-1} B, phi = Lambda^{-1} lambda, y = t^{-1} B phi  (B^T y = lambda).
  Forms: dict {bitmask: Fraction}; 1-forms dt_e given as vectors over an arbitrary set of generators.
  Omega = U^{-H} det(M)^{-1/2} Pf(-1/2 dM M^{-1} dM) ^ prod_{i} prod_{r} psi^i_r.
Numeric layer (numpy, batched): same objects, Pf(A) ^ psi via A_rs = -1/2 sum_k beta_rk ^ beta_sk,
  beta_rk = sum_q dM_rq L_qk, M^{-1} = L L^T, and determinants.
"""
from fractions import Fraction as Fr
import itertools
import math
import numpy as np


# ---------------------------------------------------------------- exact linear algebra
def det_fr(A):
    A = [[Fr(x) for x in row] for row in A]
    n = len(A)
    d = Fr(1)
    for c in range(n):
        p = None
        for r in range(c, n):
            if A[r][c] != 0:
                p = r
                break
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        piv = A[c][c]
        d *= piv
        for r in range(c + 1, n):
            f = A[r][c] / piv
            if f != 0:
                Ar, Ac = A[r], A[c]
                for k in range(c, n):
                    Ar[k] -= f * Ac[k]
    return d


def solve_fr(A, Bm):
    """A X = Bm, A n x n, Bm n x m (lists)."""
    n = len(A)
    m = len(Bm[0])
    Mx = [[Fr(x) for x in A[i]] + [Fr(x) for x in Bm[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Mx[r][c] != 0)
        Mx[c], Mx[p] = Mx[p], Mx[c]
        piv = Mx[c][c]
        Mx[c] = [x / piv for x in Mx[c]]
        for r in range(n):
            if r != c and Mx[r][c] != 0:
                f = Mx[r][c]
                Mx[r] = [a - f * b for a, b in zip(Mx[r], Mx[c])]
    return [row[n:] for row in Mx]


def inv_fr(A):
    n = len(A)
    return solve_fr(A, [[Fr(int(i == j)) for j in range(n)] for i in range(n)])


# ---------------------------------------------------------------- graph
class Net:
    def __init__(self, nv, edges):
        self.nv = nv
        self.edges = [tuple(e) for e in edges]
        self.N = len(self.edges)
        self.root = nv - 1
        self.n = self.N - nv + 1

    def Bmat(self):
        B = [[0] * (self.nv - 1) for _ in range(self.N)]
        for e, (a, b) in enumerate(self.edges):
            if a < self.root:
                B[e][a] += 1
            if b < self.root:
                B[e][b] -= 1
        return B

    def dfs_tree(self):
        adj = {v: [] for v in range(self.nv)}
        for e, (a, b) in enumerate(self.edges):
            adj[a].append((b, e))
            adj[b].append((a, e))
        for v in adj:
            adj[v].sort(key=lambda x: -x[1])
        parent = {self.root: None}

        def dfs(u):
            for (w, e) in adj[u]:
                if w not in parent:
                    parent[w] = (u, e)
                    dfs(w)

        dfs(self.root)
        tree = set(pe[1] for v, pe in parent.items() if pe is not None)
        assert len(tree) == self.nv - 1
        return tree, parent

    def cycle_basis(self):
        tree, parent = self.dfs_tree()

        def up(v):  # signed edges from v to the root along the tree
            out = []
            while parent[v] is not None:
                u, e = parent[v]
                s = 1 if self.edges[e] == (v, u) else -1   # +1 if traversed tail -> head
                out.append((v, e, s))
                v = u
            return out

        rows = []
        for c in range(self.N):
            if c in tree:
                continue
            a, b = self.edges[c]
            row = [0] * self.N
            row[c] = 1
            # from b back to a in the tree
            pb, pa = up(b), up(a)
            sb = set(e for _, e, _ in pb)
            sa = set(e for _, e, _ in pa)
            common = sb & sa
            for v, e, s in pb:
                if e not in common:
                    row[e] += s
            for v, e, s in pa:
                if e not in common:
                    row[e] -= s
            rows.append(row)
        B = self.Bmat()
        BT = [[B[e][v] for e in range(self.N)] for v in range(self.nv - 1)]
        d = det_fr(BT + rows)
        assert d != 0
        if d < 0:
            rows[0] = [-x for x in rows[0]]
        # each row must be a cycle: B^T c = 0
        for row in rows:
            for v in range(self.nv - 1):
                assert sum(B[e][v] * row[e] for e in range(self.N)) == 0
        return rows

    def currents(self, t, lam):
        """exact Kirchhoff currents y[e][i]; lam: list over v < root of H-lists"""
        B = self.Bmat()
        V = self.nv - 1
        Lam = [[sum(Fr(B[e][a] * B[e][b]) / t[e] for e in range(self.N)) for b in range(V)] for a in range(V)]
        phi = solve_fr(Lam, [[Fr(x) for x in l] for l in lam])
        H = len(lam[0])
        y = [[sum(B[e][v] * phi[v][i] for v in range(V)) / t[e] for i in range(H)] for e in range(self.N)]
        # checks: B^T y = lam, cycle law sum_e C_re t_e y_e = 0 is implied by y = t^{-1} B phi
        for v in range(V):
            for i in range(H):
                assert sum(B[e][v] * y[e][i] for e in range(self.N)) == lam[v][i]
        return y


# ---------------------------------------------------------------- exterior algebra
def _sign(a, b):
    """sign of dx_a ^ dx_b -> dx_{a|b} (a, b disjoint bitmasks)"""
    s = 0
    bb = b
    while bb:
        j = (bb & -bb).bit_length() - 1
        s += bin(a >> (j + 1)).count("1")
        bb &= bb - 1
    return -1 if s & 1 else 1


def wedge(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            if a & b:
                continue
            k = a | b
            out[k] = out.get(k, 0) + _sign(a, b) * x * y
    return {k: v for k, v in out.items() if v != 0}


def fadd(f, g, c=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def fscale(f, c):
    return {k: c * v for k, v in f.items() if c * v != 0}


def oneform(vec):
    """vec: list over generators -> 1-form"""
    return {1 << j: Fr(c) for j, c in enumerate(vec) if c != 0}


def pfaffian_forms(A):
    """Pf of an antisymmetric matrix of even forms, expansion along the first row:
    Pf(A) = sum_{j>=2} (-1)^j A_{1j} Pf(A_{^1^j}) (1-based)."""
    m = len(A)
    if m == 0:
        return {0: Fr(1)}
    if m % 2:
        return {}
    out = {}
    for j in range(1, m):
        if not A[0][j]:
            continue
        idx = [k for k in range(1, m) if k != j]
        sub = [[A[a][b] for b in idx] for a in idx]
        term = wedge(A[0][j], pfaffian_forms(sub))
        out = fadd(out, term, (-1) ** (j + 1))   # (-1)^j with 1-based j+1
    return out


def rho_forms(C, t, dts):
    """returns (Pf(-1/2 dM M^{-1} dM) as a form, U); dts[e] = 1-form of dt_e."""
    n = len(C)
    N = len(t)
    M = [[sum(C[r][e] * C[s][e] * Fr(t[e]) for e in range(N)) for s in range(n)] for r in range(n)]
    U = det_fr(M)
    Mi = inv_fr(M)
    dM = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            f = {}
            for e in range(N):
                c = C[r][e] * C[s][e]
                if c:
                    f = fadd(f, dts[e], c)
            dM[r][s] = f
    # X = M^{-1} dM
    X = [[{} for _ in range(n)] for _ in range(n)]
    for p in range(n):
        for s in range(n):
            f = {}
            for q in range(n):
                if Mi[p][q] != 0:
                    f = fadd(f, dM[q][s], Mi[p][q])
            X[p][s] = f
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            f = {}
            for p in range(n):
                f = fadd(f, wedge(dM[r][p], X[p][s]))
            A[r][s] = fscale(f, Fr(-1, 2))
    return pfaffian_forms(A), U, M


def psi_forms(C, y, dts):
    """psi^i_r = sum_e C_re y^i_e dt_e, order i outer, r inner"""
    n = len(C)
    N = len(y)
    H = len(y[0])
    out = []
    for i in range(H):
        for r in range(n):
            f = {}
            for e in range(N):
                if C[r][e] and y[e][i] != 0:
                    f = fadd(f, dts[e], C[r][e] * y[e][i])
            out.append(f)
    return out


def omega(net, C, t, lam, dts, T=1, powers=None):
    """Omega as a form (without the factor det(M)^{-1/2} when T = 1): returns (form, U, y)."""
    y = net.currents(t, lam)
    H = len(lam[0])
    n = len(C)
    if T == 1:
        pf, U, M = rho_forms(C, t, dts)
        f = pf
    else:
        M = [[sum(C[r][e] * C[s][e] * Fr(t[e]) for e in range(len(t))) for s in range(n)] for r in range(n)]
        U = det_fr(M)
        f = {0: Fr(1)}
    for p in psi_forms(C, y, dts):
        f = wedge(f, p)
    f = fscale(f, Fr(1) / U ** H)
    return f, U, y


def density_slice(net, C, t, lam, T=1, e0=None):
    """D(t) = coef * U^{-1/2} (T=1) or coef (T=0), from D = (-1)^{e0}/t_{e0} [Omega]_{omit dt_e0}."""
    N = net.N
    if e0 is None:
        e0 = N - 1
    gens = [e for e in range(N) if e != e0]
    dts = []
    for e in range(N):
        vec = [0] * (N - 1)
        if e != e0:
            vec[gens.index(e)] = 1
        dts.append(oneform(vec))
    f, U, y = omega(net, C, t, lam, dts, T=T)
    top = f.get((1 << (N - 1)) - 1, Fr(0))
    coef = Fr((-1) ** e0) * top / Fr(t[e0])
    return coef, U, y


# ---------------------------------------------------------------- numeric (batched)
def pairings(n):
    """perfect matchings of range(n) with the Pfaffian sign"""
    if n == 0:
        yield 1, []
        return
    for j in range(1, n):
        rest = [k for k in range(1, n) if k != j]
        for s, m in pairings(len(rest)):
            yield (-1) ** (j + 1) * s, [(0, j)] + [(rest[a], rest[b]) for a, b in m]


def density_batch(nv, edges, C, t, lam, T=1, return_y=False):
    """t: (B, N) array; lam: (nv-1, H) array. Returns D(t) (slice formula, e0 = N-1)."""
    N = len(edges)
    V = nv - 1
    Bm = np.zeros((N, V))
    for e, (a, b) in enumerate(edges):
        if a < V:
            Bm[e, a] += 1
        if b < V:
            Bm[e, b] -= 1
    C = np.array(C, dtype=float)
    n = C.shape[0]
    lam = np.asarray(lam, dtype=float)
    H = lam.shape[1]
    it = 1.0 / t
    Lam = np.einsum("ev,be,ew->bvw", Bm, it, Bm)
    phi = np.linalg.solve(Lam, np.broadcast_to(lam, (t.shape[0],) + lam.shape))
    y = np.einsum("ev,bvh->beh", Bm, phi) * it[:, :, None]
    M = np.einsum("re,be,se->brs", C, t, C)
    U = np.linalg.det(M)
    # psi rows (i outer, r inner), columns e = 0..N-2
    P = np.einsum("re,beh->bhre", C, y).reshape(t.shape[0], H * n, N)[:, :, : N - 1]
    if T == 0:
        top = np.linalg.det(P)
        D = (-1) ** (N - 1) * top / t[:, N - 1] / U ** H
    else:
        Minv = np.linalg.inv(M)
        ev, Vv = np.linalg.eigh(Minv)                      # Minv = L L^T, L = V sqrt(ev)
        L = Vv * np.sqrt(np.clip(ev, 0, None))[:, None, :]
        CL = np.einsum("qe,bqk->bek", C, L)                # (B, N, k)
        beta = C[None, :, :, None] * CL[:, None, :, :]     # (B, r, e, k)
        beta = beta[:, :, : N - 1, :]
        top = np.zeros(t.shape[0])
        npair = n // 2
        for sgn, m in pairings(n):
            for ks in itertools.product(range(n), repeat=npair):
                rows = []
                for (r, s), k in zip(m, ks):
                    rows.append(beta[:, r, :, k])
                    rows.append(beta[:, s, :, k])
                Mat = np.concatenate([np.stack(rows, axis=1), P], axis=1)
                top += sgn * np.linalg.det(Mat)
        top *= (-0.5) ** npair
        D = (-1) ** (N - 1) * top / t[:, N - 1] / U ** H / np.sqrt(U)
    if return_y:
        return D, y
    return D
