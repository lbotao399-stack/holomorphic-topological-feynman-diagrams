"""Independent layer for the verification of msym-D8 (does not import htk or the explorer's code).

Graphs (built from the descriptions in the claims, H = 4, T = 0):
  pentagon: O -> 1 -> 2 -> 3 -> 4 -> O; outputs 1..4 = vertices 0..3, O = vertex 4 (root);
            edges e0=(O,0), e1=(0,1), e2=(1,2), e3=(2,3), e4=(3,O).
  theta(p,q,r): O = root (last vertex), N = vertex 0; paths P1,P2,P3 from O to N with p,q,r edges, oriented away
            from O, edges listed P1, P2, P3; interior vertices 1,2,... numbered along P1, then P2, then P3.

Kernel (T = 0, H = 4): K = int_{t>0, sum t = 1} D(t) exp(-sum_e y_e.z_e) dt_0..dt_{N-2},
  wedge_r d^4 l_r = D(t) iota_E (dt_0 ^ .. ^ dt_{N-1}),  l_r loop currents of the fundamental cycles of a BFS tree
  (any tree: l' = A l + const, det A = +-1, (det A)^4 = 1).
Mesh solver: y = y0 + C^T l, (C t C^T) l = -C t y0;  d l / d t_j = -M^{-1} C[:, j] y_j.
D(t) = (-1)^{e0} / t_{e0} det[ d l_{r,a} / d t_j ]_{(r,a) (r outer); j != e0},  e0 = N-1.
"""
from fractions import Fraction as Fr
from collections import deque
import itertools
import numpy as np

H = 4


def pentagon():
    return 5, [(4, 0), (0, 1), (1, 2), (2, 3), (3, 4)]


def theta(p, q, r):
    nv = p + q + r - 1
    O, N = nv - 1, 0
    edges, paths, nxt = [], [], 1
    for L in (p, q, r):
        prev, eids, ints = O, [], []
        for j in range(L):
            cur = N if j == L - 1 else nxt
            if j < L - 1:
                ints.append(cur)
                nxt += 1
            eids.append(len(edges))
            edges.append((prev, cur))
            prev = cur
        paths.append((eids, ints))
    return nv, edges, paths


# ------------------------------------------------------------------ exact linear algebra

def solve_exact(A, B):
    """A X = B, A square (list of lists of Fractions), B list of columns."""
    n = len(A)
    M = [list(A[i]) + [b[i] for b in B] for i in range(n)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[p] = M[p], M[i]
        pv = M[i][i]
        M[i] = [x / pv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [[M[i][n + k] for i in range(n)] for k in range(len(B))]


def det_exact(M):
    M = [list(r) for r in M]
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
            if M[r][i] != 0:
                f = M[r][i] / M[i][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return d


# ------------------------------------------------------------------ mesh-current network

class Network:
    def __init__(self, nv, edges):
        self.nv, self.edges, self.N = nv, edges, len(edges)
        self.root = nv - 1
        adj = {v: [] for v in range(nv)}
        for e, (a, b) in enumerate(edges):
            adj[a].append((e, b, +1))     # traversing e from a to b: along the orientation
            adj[b].append((e, a, -1))
        parent = {self.root: None}
        order = [self.root]
        dq = deque([self.root])
        while dq:
            v = dq.popleft()
            for e, w, s in adj[v]:
                if w not in parent:
                    parent[w] = (e, v, s)
                    order.append(w)
                    dq.append(w)
        self.parent = parent
        self.order = order
        self.tree = sorted(parent[v][0] for v in parent if parent[v] is not None)
        self.cotree = [e for e in range(self.N) if e not in self.tree]
        self.n = len(self.cotree)
        # fundamental cycle of cotree edge c = (a -> b): c, then the tree path b -> a
        C = [[0] * self.N for _ in self.cotree]
        for r, c in enumerate(self.cotree):
            a, b = edges[c]
            C[r][c] = 1
            for e, s in self.tree_path(b, a):
                C[r][e] += s
        self.C = C

    def path_to_root(self, v):
        out = []
        while self.parent[v] is not None:
            e, u, s = self.parent[v]          # u -> v via e with sign s (s = +1 if e is oriented u -> v)
            out.append((e, -s))               # traversing v -> u
            v = u
        return out

    def tree_path(self, a, b):
        """Edges (e, sign) of the tree path a -> b (sign +1 when traversed along e's orientation)."""
        pa = self.path_to_root(a)
        pb = self.path_to_root(b)
        ea = [e for e, _ in pa]
        eb = [e for e, _ in pb]
        while ea and eb and ea[-1] == eb[-1]:
            ea.pop(); eb.pop(); pa.pop(); pb.pop()
        return pa + [(e, -s) for e, s in reversed(pb)]

    def tree_currents(self, lam):
        """y0 supported on the tree with B^T y0 = lam (outflow - inflow = lam_v at v != root)."""
        y0 = [[Fr(0)] * H for _ in range(self.N)]
        sub = {v: [Fr(x) for x in lam.get(v, [0] * H)] for v in range(self.nv)}
        for v in reversed(self.order):
            if self.parent[v] is None:
                continue
            e, u, s = self.parent[v]
            # the subtree of v injects sub[v]; it flows out of v towards u through e
            for a in range(H):
                y0[e][a] = s * (-sub[v][a]) if s == +1 else sub[v][a]
            for a in range(H):
                sub[u][a] += sub[v][a]
        # orientation: e = (u -> v) (s=+1): current along e is from u to v, i.e. into the subtree: -sub[v]
        return y0

    def solve(self, t, lam):
        t = [Fr(x) for x in t]
        y0 = self.tree_currents(lam)
        Mm = [[sum(self.C[r][e] * t[e] * self.C[s][e] for e in range(self.N)) for s in range(self.n)] for r in range(self.n)]
        rhs = [[-sum(self.C[r][e] * t[e] * y0[e][a] for e in range(self.N)) for r in range(self.n)] for a in range(H)]
        lsol = solve_exact(Mm, rhs)
        ell = [[lsol[a][r] for a in range(H)] for r in range(self.n)]
        y = [[y0[e][a] + sum(self.C[r][e] * ell[r][a] for r in range(self.n)) for a in range(H)] for e in range(self.N)]
        return y, ell, Mm

    def check_kcl(self, y, lam):
        for v in range(self.nv - 1):
            for a in range(H):
                s = sum((y[e][a] if p == v else 0) - (y[e][a] if q == v else 0) for e, (p, q) in enumerate(self.edges))
                if s != Fr(lam.get(v, [0] * H)[a]):
                    return False
        for r in range(self.n):
            for a in range(H):
                pass
        return True

    def density(self, t, lam):
        """D(t) exactly, and the currents y."""
        t = [Fr(x) for x in t]
        y, ell, Mm = self.solve(t, lam)
        cols = [[self.C[r][j] for r in range(self.n)] for j in range(self.N)]
        Minvcols = solve_exact(Mm, cols)                     # M^{-1} C[:, j]
        e0 = self.N - 1
        rows = []
        for r in range(self.n):
            for a in range(H):
                rows.append([-Minvcols[j][r] * y[j][a] for j in range(self.N) if j != e0])
        D = Fr((-1) ** e0) / t[e0] * det_exact(rows)
        return D, y

    # -------------------------------------------------------------- float, vectorized
    def density_float(self, T, lam):
        """T: P x N positive times; lam: dict v -> 4-vector. Returns D (P,), y (P, N, H)."""
        P = T.shape[0]
        C = np.array(self.C, dtype=float)                     # n x N
        y0 = np.array([[float(x) for x in row] for row in self.tree_currents(lam)])   # N x H
        Mm = np.einsum('re,pe,se->prs', C, T, C)
        rhs = -np.einsum('re,pe,ea->pra', C, T, y0)
        ell = np.linalg.solve(Mm, rhs)                        # P x n x H
        y = y0[None, :, :] + np.einsum('re,pra->pea', C, ell)
        Minv = np.linalg.inv(Mm)                              # P x n x n
        MC = np.einsum('prs,sj->prj', Minv, C)               # P x n x N
        e0 = self.N - 1
        J = -MC[:, :, None, :] * np.transpose(y, (0, 2, 1))[:, None, :, :]   # P x n x H x N
        rows = J.reshape(P, self.n * H, self.N)[:, :, :e0]
        D = (-1) ** e0 / T[:, e0] * np.linalg.det(rows)
        return D, y


def bracket(lam, S):
    from sympy import Matrix
    return int(Matrix([list(lam[v]) for v in S]).det())


def perm_sign(seq):
    s = list(seq)
    sign = 1
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            if s[i] > s[j]:
                sign = -sign
    return sign


# ------------------------------------------------------------------ quadrature on products of simplices

def gauss_simplex(k, n):
    """Nodes (M x k) and weights for the open simplex {b_1..b_k > 0, sum b < 1} (Lebesgue measure db_1..db_k),
    Duffy (collapsed) map from the cube with Gauss-Legendre nodes, n points per direction."""
    if k == 0:
        return np.zeros((1, 0)), np.ones(1)
    x, w = np.polynomial.legendre.leggauss(n)
    x = (x + 1) / 2
    w = w / 2
    grids = np.array(list(itertools.product(range(n), repeat=k)))
    U = x[grids]
    W = np.prod(w[grids], axis=1)
    # b_1 = u_1, b_2 = (1-u_1) u_2, ..., b_k = prod_{i<k}(1-u_i) u_k ; Jacobian prod_i (1-u_i)^{k-i}
    B = np.zeros_like(U)
    rem = np.ones(U.shape[0])
    jac = np.ones(U.shape[0])
    for i in range(k):
        B[:, i] = rem * U[:, i]
        jac *= rem
        rem = rem * (1 - U[:, i])
    return B, W * jac


def product_rule(dims, n):
    """Tensor product of gauss_simplex rules; returns list of node arrays per block and the weights."""
    rules = [gauss_simplex(k, n) for k in dims]
    idx = list(itertools.product(*[range(len(r[1])) for r in rules]))
    idx = np.array(idx)
    nodes = [rules[b][0][idx[:, b]] for b in range(len(dims))]
    W = np.ones(len(idx))
    for b in range(len(dims)):
        W = W * rules[b][1][idx[:, b]]
    return nodes, W
