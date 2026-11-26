"""Independent brute force for two-loop T = 1 kernels (no htk import).

Definition (task setting, htk CONVENTIONS Section 2):
  B_{ev} = delta_{v,e(0)} - delta_{v,e(1)} (v != root = nv-1),  B^T y = lambda,  t_e y_e = phi_{e(0)} - phi_{e(1)},
  C: integer cycle basis with det[B^T ; C] > 0,  M = C t C^T,  U = det M,
  rho_2 = U^{-1/2} * ( -1/2 sum_{ab} dM_{1a} (M^{-1})_{ab} ^ dM_{b2} ),
  psi^i_r = sum_e C_re y^i_e dt_e,
  Omega = U^{-H} rho_2 ^ psi^1_1 ^ psi^1_2 ^ ... ^ psi^H_1 ^ psi^H_2 = D(t) iota_E (dt_0 ^ ... ^ dt_{N-1}),
  K(lam; z) = int_{sum t = 1} D(t) exp(-sum_e y_e.z_e) dt_0..dt_{N-2}.
Here D(t) = qD * U^{-1/2} with qD rational; returned as (qD, U).
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations
import random


# ---------------------------------------------------------------- exact linear algebra
def det(A):
    A = [list(map(Fr, r)) for r in A]
    n = len(A)
    s = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            s = -s
        s *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c] / A[c][c]
                A[r] = [A[r][k] - f * A[c][k] for k in range(n)]
    return s


def solve(A, b):
    """A x = b, b a list of column vectors given as list of rows (n x m)."""
    n = len(A)
    m = len(b[0])
    Aug = [list(map(Fr, A[i])) + list(map(Fr, b[i])) for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Aug[r][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        piv = Aug[c][c]
        Aug[c] = [v / piv for v in Aug[c]]
        for r in range(n):
            if r != c and Aug[r][c]:
                f = Aug[r][c]
                Aug[r] = [Aug[r][k] - f * Aug[c][k] for k in range(n + m)]
    return [row[n:] for row in Aug]


def inv(A):
    n = len(A)
    return solve(A, [[Fr(int(i == j)) for j in range(n)] for i in range(n)])


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


# ---------------------------------------------------------------- exterior algebra on dt_0..dt_{N-1}
def fwedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            if set(ka) & set(kb):
                continue
            seq = ka + kb
            key = tuple(sorted(seq))
            out[key] = out.get(key, 0) + perm_sign(seq) * va * vb
    return {k: v for k, v in out.items() if v != 0}


def fadd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


# ---------------------------------------------------------------- graph
class Net:
    def __init__(self, nv, edges):
        self.nv, self.edges = nv, list(edges)
        self.N = len(edges)
        self.root = nv - 1
        self.n = self.N - nv + 1
        self.B = [[(1 if a == v else 0) - (1 if b == v else 0) for v in range(nv - 1)] for (a, b) in edges]
        self.C = self._cycles()

    def _cycles(self):
        # BFS spanning tree from the root (independent of Kruskal); fundamental cycles of the cotree edges
        adj = {v: [] for v in range(self.nv)}
        for e, (a, b) in enumerate(self.edges):
            adj[a].append((e, b)); adj[b].append((e, a))
        par = {self.root: None}
        order = [self.root]
        tree = set()
        for v in order:
            for e, u in adj[v]:
                if u not in par:
                    par[u] = (e, v); tree.add(e); order.append(u)

        def path_to_root(v):
            out = []   # list of (edge, sign) for traversing from v toward root
            while par[v] is not None:
                e, p = par[v]
                a, b = self.edges[e]
                out.append((e, +1 if a == v else -1))
                v = p
            return out
        C = []
        for e, (a, b) in enumerate(self.edges):
            if e in tree:
                continue
            row = [0] * self.N
            row[e] += 1                      # traverse a -> b along e, then b -> root, root -> a
            for f, s in path_to_root(b):
                row[f] += s
            for f, s in path_to_root(a):
                row[f] -= s
            C.append(row)
        # check circulation: B^T c = 0
        for row in C:
            for v in range(self.nv - 1):
                assert sum(self.B[e][v] * row[e] for e in range(self.N)) == 0
        Mtop = [[self.B[e][v] for e in range(self.N)] for v in range(self.nv - 1)] + C
        d = det(Mtop)
        assert d != 0
        if d < 0:
            C[0] = [-x for x in C[0]]
        return C

    def currents(self, t, lam):
        """lam: list over non-root vertices of H-vectors; returns y[e][i]."""
        g = self.nv - 1
        H = len(lam[0])
        L = [[sum(Fr(self.B[e][u] * self.B[e][v]) / t[e] for e in range(self.N)) for v in range(g)] for u in range(g)]
        phi = solve(L, [[Fr(c) for c in lam[v]] for v in range(g)])
        y = [[sum(self.B[e][v] * phi[v][i] for v in range(g)) / t[e] for i in range(H)] for e in range(self.N)]
        # Kirchhoff check
        for v in range(g):
            for i in range(H):
                assert sum(self.B[e][v] * y[e][i] for e in range(self.N)) == lam[v][i]
        for row in self.C:
            for i in range(H):
                assert sum(row[e] * t[e] * y[e][i] for e in range(self.N)) == 0
        return y

    def density(self, t, lam):
        """(qD, U): D(t) = qD * U^{-1/2}."""
        assert self.n == 2
        H = len(lam[0])
        N, C = self.N, self.C
        y = self.currents(t, lam)
        M = [[sum(C[r][e] * C[s][e] * t[e] for e in range(N)) for s in range(2)] for r in range(2)]
        U = det(M)
        Mi = inv(M)
        dMf = [[{} for _ in range(2)] for _ in range(2)]
        for r in range(2):
            for s in range(2):
                f = {}
                for e in range(N):
                    if C[r][e] * C[s][e]:
                        f[(e,)] = Fr(C[r][e] * C[s][e])
                dMf[r][s] = f
        A12 = {}
        for a in range(2):
            for b in range(2):
                A12 = fadd(A12, {k: v * Mi[a][b] for k, v in fwedge(dMf[0][a], dMf[b][1]).items()})
        R2 = {k: Fr(-1, 2) * v for k, v in A12.items()}          # rho_2 = U^{-1/2} R2
        form = R2
        for i in range(H):
            for r in range(2):
                psi = {(e,): C[r][e] * y[e][i] for e in range(N) if C[r][e] * y[e][i] != 0}
                form = fwedge(form, psi)
        top = form.get(tuple(range(N - 1)), Fr(0))
        qD = top / ((-1) ** (N - 1) * t[N - 1]) / U ** H
        # basicness check: coefficient extraction with another omitted edge gives the same D
        for e0 in range(N - 1):
            key = tuple(x for x in range(N) if x != e0)
            alt = form.get(key, Fr(0)) / ((-1) ** e0 * t[e0]) / U ** H
            assert alt == qD, (alt, qD)
        return qD, U, y


# ---------------------------------------------------------------- theta graphs
def theta_rootN(p, q, r):
    """Explorer's theta_rootN data (kernels_check.py): interior vertices along A, B, C numbered first,
    N' = nv-2, N = root = nv-1; edges oriented along the paths from N; edge order A, B, C."""
    nv = p + q + r - 1
    N, Np = nv - 1, nv - 2
    edges, paths, nxt = [], [], 0
    for n_ in (p, q, r):
        prev, pe = N, []
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                nxt += 1
            pe.append(len(edges))
            edges.append((prev, cur))
            prev = cur
        paths.append(pe)
    return nv, edges, paths


# S2 Section 7 data:  (A) vertices 1..5 -> 0..4, O = 5; (B) u,w,q,a,b -> 0..4, O = 5
S2_A = (6, [(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)], [[3], [0, 1, 2], [6, 5, 4]])
S2_B = (6, [(5, 0), (0, 2), (5, 1), (1, 2), (2, 3), (3, 4), (4, 5)], [[0, 1], [2, 3], [6, 5, 4]])


class Theta:
    """Theta graph with N = root; paths (edge lists) from N to N'."""

    def __init__(self, nv, edges, paths, H):
        self.net = Net(nv, edges)
        self.H = H
        self.paths = paths
        self.nP = [len(p) for p in paths]
        self.sg, self.inter = [], []
        for pe in paths:
            v, sg, ints = nv - 1, [], []
            for pos, e in enumerate(pe):
                a, b = edges[e]
                if a == v:
                    sg.append(1); v = b
                else:
                    assert b == v
                    sg.append(-1); v = a
                if pos < len(pe) - 1:
                    ints.append(v)
            self.sg.append(sg); self.inter.append(ints)
        self.Np = v
        self.alpha = [Fr(2 * H + 1, 2) - n for n in self.nP]

    # coordinates xi = (wB, wC, x^A_2.., x^B_2.., x^C_2..)
    def point(self, rng):
        def simp(k):
            a = [rng.randint(1, 40) for _ in range(k)]
            return [Fr(v, sum(a)) for v in a]
        return simp(3), [simp(n) for n in self.nP]

    def t_of(self, w, x):
        t = [None] * self.net.N
        for P in range(3):
            for j, e in enumerate(self.paths[P]):
                t[e] = x[P][j] / w[P]
        return t

    def jac(self, w, x):
        """|det[t ; dt/dxi]| for t^P_j = x^P_j / w_P, w_A = 1 - wB - wC, x^P_1 = 1 - sum_{j>=2} x^P_j."""
        N = self.net.N
        rows = [self.t_of(w, x)]
        # d/d wB, d/d wC
        for Q in (1, 2):
            row = [Fr(0)] * N
            for P in range(3):
                dwP = (1 if P == Q else 0) - (1 if P == 0 else 0)
                for j, e in enumerate(self.paths[P]):
                    row[e] = -x[P][j] / w[P] ** 2 * dwP
            rows.append(row)
        for P in range(3):
            for l in range(1, self.nP[P]):
                row = [Fr(0)] * N
                for j, e in enumerate(self.paths[P]):
                    row[e] = (Fr(int(j == l)) - Fr(int(j == 0))) / w[P]
                rows.append(row)
        return abs(det(rows))

    def Dtheta(self, lam):
        H = self.H
        A, B, C = self.inter
        rows = []
        for i in range(H):
            rows.append([lam[v][i] for v in A] + [lam[v][i] for v in B] + [0] * len(C))
        for i in range(H):
            rows.append([0] * len(A) + [lam[v][i] for v in B] + [lam[v][i] for v in C])
        return det(rows)

    def claimed(self, w, lam):
        """(coefficient, radicand) of (D/2) prod_P w_P^{H+1/2-n_P} = coeff * radicand^{-1/2}."""
        c = self.Dtheta(lam) / 2
        rad = Fr(1)
        for P in range(3):
            c *= w[P] ** (self.H + 1 - self.nP[P])
            rad *= w[P]
        return c, rad

    def check(self, npts, seed, lam_fixed=None, zero_Np=False):
        """returns the set of signs s with  g = s * claimed  (squares compared exactly)."""
        rng = random.Random(seed)
        g = self.net.nv - 1
        signs = set()
        for _ in range(npts):
            lam = lam_fixed or [[Fr(rng.randint(-7, 7)) for _ in range(self.H)] for _ in range(g)]
            if zero_Np and self.Np != self.net.root:
                lam[self.Np] = [Fr(0)] * self.H
            w, x = self.point(rng)
            t = self.t_of(w, x)
            qD, U, y = self.net.density(t, lam)
            J = self.jac(w, x)
            lhs_c, lhs_r = qD * J, U
            c, rad = self.claimed(w, lam)
            if c == 0:
                assert lhs_c == 0
                continue
            assert lhs_c ** 2 / lhs_r == c ** 2 / rad, "magnitude mismatch"
            signs.add((1 if lhs_c > 0 else -1) * (1 if c > 0 else -1))
        return signs

    def X_rule(self, w, x, zP):
        """X_v for interior vertices and N' (S3 eq-series/eq-parallel based rule):
        X_{P_k} = sigma^P_k z_P + (1 - sigma^P_k) Zbar, sigma^P_k = sum_{j>k} x^P_j, X_{N'} = Zbar,
        Zbar = sum_P w_P z_P, z_P the displacement of the first edge of P in the direction from N."""
        H = self.H
        Zb = [sum(w[P] * zP[P][i] for P in range(3)) for i in range(H)]
        X = {}
        for P in range(3):
            for k, v in enumerate(self.inter[P], start=1):
                sig = sum(x[P][j] for j in range(k, self.nP[P]))
                X[v] = [sig * zP[P][i] + (1 - sig) * Zb[i] for i in range(H)]
        X[self.Np] = Zb
        return X


# ---------------------------------------------------------------- Dirichlet integrals with half-integers
from math import factorial


def gamma_half(a):
    """Gamma(a) for a in (1/2) Z_{>0}: returns (rational q, k) with Gamma(a) = q * sqrt(pi)^k."""
    a = Fr(a)
    if a.denominator == 1:
        return Fr(factorial(int(a) - 1)), 0
    assert a.denominator == 2 and a > 0
    m = int(a - Fr(1, 2))
    return Fr(factorial(2 * m), 4 ** m * factorial(m)), 1


def dirichlet(exps):
    """int_{Delta^m} prod_{j=0}^m b_j^{a_j} db_1..db_m = prod Gamma(a_j+1) / Gamma(m+1+sum a_j);
    returns (q, k): q * sqrt(pi)^k."""
    q, k = Fr(1), 0
    for a in exps:
        g, kk = gamma_half(Fr(a) + 1)
        q *= g; k += kk
    g, kk = gamma_half(len(exps) + sum(Fr(a) for a in exps))
    return q / g, k - kk
