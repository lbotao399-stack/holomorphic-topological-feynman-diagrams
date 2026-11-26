"""Independent kernel layer for the verification of explore/nonSP (no import of htk, no explorer code).

Definition used (CONVENTIONS.md Sec. 1-2 restated from the task text):
  B_{ev} = delta_{v,e0} - delta_{v,e1} (v < root = nv-1);  t y = B phi,  B^T y = lam.
  C = fundamental cycle matrix of a BFS spanning tree, first row flipped so that det[B^T; C] > 0.
  M = C t C^T, U = det M, Psi_{(i,r),e} = C_re y^i_e  (i outer, r inner).
  Omega = U^{-H} prod_{i,r} psi^i_r = D iota_E vol,   D = det[1; Psi] / (U^H sum_e t_e).
  K = int_{t>0, sum t = 1} D dt_0..dt_{N-2}  (= int_{t_k = 1} D, any k, D homogeneous of degree -N).
"""
from fractions import Fraction as Fr
import itertools
import numpy as np


class G:
    def __init__(self, nv, edges):
        self.nv, self.edges = nv, [tuple(e) for e in edges]
        self.N = len(edges)
        self.root = nv - 1
        self.n = self.N - nv + 1
        self.B = [[(1 if v == a else 0) - (1 if v == b else 0) for v in range(nv - 1)] for (a, b) in self.edges]
        self.C = self._cycles()

    def _tree(self):
        adj = {v: [] for v in range(self.nv)}
        for e, (a, b) in enumerate(self.edges):
            adj[a].append((b, e)); adj[b].append((a, e))
        seen = {0: None}; order = [0]; tree = []
        while order:
            v = order.pop(0)
            for (w, e) in adj[v]:
                if w not in seen:
                    seen[w] = (v, e); tree.append(e); order.append(w)
        return tree, seen

    def _cycles(self):
        tree, par = self._tree()
        cot = [e for e in range(self.N) if e not in tree]
        def path_to_root(v):
            p = []
            while par[v] is not None:
                u, e = par[v]; p.append((v, u, e)); v = u
            return p
        C = []
        for c in cot:
            a, b = self.edges[c]
            row = [0] * self.N
            row[c] = 1
            # cycle: a -> b along c, then b -> a in the tree: b -> root then root -> a
            for (v, u, e) in path_to_root(b):   # traverse v -> u
                row[e] += 1 if self.edges[e] == (v, u) else -1
            for (v, u, e) in path_to_root(a):   # traverse u -> v (reverse)
                row[e] += 1 if self.edges[e] == (u, v) else -1
            C.append(row)
        Mx = [[self.B[e][v] for e in range(self.N)] for v in range(self.nv - 1)] + C
        if fdet(Mx) < 0:
            C[0] = [-x for x in C[0]]
        Mx = [[self.B[e][v] for e in range(self.N)] for v in range(self.nv - 1)] + C
        assert fdet(Mx) > 0
        # each row is a closed cycle
        for row in C:
            for v in range(self.nv - 1):
                assert sum(row[e] * self.B[e][v] for e in range(self.N)) == 0
        return C


def fdet(A):
    A = [[Fr(x) for x in r] for r in A]
    n = len(A); d = Fr(1)
    for i in range(n):
        p = next((k for k in range(i, n) if A[k][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            A[i], A[p] = A[p], A[i]; d = -d
        d *= A[i][i]
        for k in range(i + 1, n):
            f = A[k][i] / A[i][i]
            if f:
                A[k] = [A[k][j] - f * A[i][j] for j in range(n)]
    return d


def fsolve(A, b):
    n = len(A)
    M = [[Fr(x) for x in A[i]] + [Fr(x) for x in b[i]] for i in range(n)]
    m = len(b[0])
    for i in range(n):
        p = next(k for k in range(i, n) if M[k][i] != 0)
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [x / piv for x in M[i]]
        for k in range(n):
            if k != i and M[k][i] != 0:
                f = M[k][i]
                M[k] = [M[k][j] - f * M[i][j] for j in range(n + m)]
    return [M[i][n:] for i in range(n)]


def currents_exact(g, t, lam):
    """lam: list over non-root vertices of H-vectors. returns y[e][i]."""
    H = len(lam[0])
    nv1 = g.nv - 1
    L = [[sum(Fr(g.B[e][u] * g.B[e][v]) / t[e] for e in range(g.N)) for v in range(nv1)] for u in range(nv1)]
    phi = fsolve(L, [[Fr(x) for x in lam[v]] for v in range(nv1)])
    y = [[sum(g.B[e][v] * phi[v][i] for v in range(nv1)) / t[e] for i in range(H)] for e in range(g.N)]
    # Kirchhoff checks
    for v in range(nv1):
        for i in range(H):
            assert sum(g.B[e][v] * y[e][i] for e in range(g.N)) == lam[v][i]
    for row in g.C:
        for i in range(H):
            assert sum(row[e] * t[e] * y[e][i] for e in range(g.N)) == 0
    return y


def density_exact(g, t, lam):
    H = len(lam[0])
    y = currents_exact(g, t, lam)
    n, N = g.n, g.N
    M = [[sum(g.C[r][e] * g.C[s][e] * t[e] for e in range(N)) for s in range(n)] for r in range(n)]
    U = fdet(M)
    Psi = [[g.C[r][e] * y[e][i] for e in range(N)] for i in range(H) for r in range(n)]
    return fdet([[1] * N] + Psi) / (U ** H * sum(t))


# ---------------- numpy batch ----------------
class NG:
    def __init__(self, g, lams):
        """lams: array (S, nv-1, H) of source configurations evaluated simultaneously."""
        self.g = g
        self.B = np.array(g.B, dtype=float)          # (N, nv-1)
        self.C = np.array(g.C, dtype=float)          # (n, N)
        self.lams = np.asarray(lams, dtype=float)    # (S, nv-1, H)
        self.H = self.lams.shape[2]

    def D(self, t):
        """t: (P, N) -> D: (P, S)"""
        B, C = self.B, self.C
        P, N = t.shape
        n, H = self.g.n, self.H
        it = 1.0 / t
        L = np.einsum('ev,pe,ew->pvw', B, it, B)
        S = self.lams.shape[0]
        rhs = np.transpose(self.lams, (1, 0, 2)).reshape(1, B.shape[1], S * H)
        phi = np.linalg.solve(L, np.broadcast_to(rhs, (P,) + rhs.shape[1:]))   # (P, nv-1, S*H)
        y = np.einsum('ev,pvk->pek', B, phi) * it[:, :, None]                     # (P, N, S*H)
        y = y.reshape(P, N, S, H)
        M = np.einsum('re,pe,se->prs', C, t, C)
        U = np.linalg.det(M)
        # Psi[p, s, i*n + r, e] = C[r,e] y[p,e,s,i]
        Psi = np.einsum('re,pesi->psire', C, y).reshape(P, S, H * n, N)
        A = np.concatenate([np.ones((P, S, 1, N)), Psi], axis=2)
        d = np.linalg.det(A)
        return d / (U[:, None] ** H * t.sum(1)[:, None])
