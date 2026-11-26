"""tb.py -- exact brute force for the topological factor, general T.  Own code; htk is not imported.

Conventions (identical to code/htk/CONVENTIONS.md sections 1-2):
  Graph(nv, edges): vertices 0..nv-1, root = nv-1, edges[e] = (e0, e1).
  B[e][v] = [v==e0] - [v==e1]  (v < root).
  C = fundamental cycles of the Kruskal tree (edge order); row r belongs to the r-th cotree edge
      (increasing), C[r][c_r] = +1; first row multiplied by sgn det[B^T; C_T]  ->  det[B^T; C] > 0.
  One topological direction: potentials tau_v (tau_root = 0), w_e = tau_e/(2 t_e),
      eta_e = tau_e/(2 sqrt t_e), Theta_e = pi^{-1/2} exp(-eta_e^2) d eta_e, outputs exp(kappa . tau).
  Gaussian weight exp(-tau^T Lam tau / 4 + kappa . tau), Lam = B^T t^{-1} B:
      mean tau = 2 Lam^{-1} kappa, cov tau = 2 Lam^{-1};
      wbar = Kirchhoff currents with sources kappa;  P = sum_e t_e wbar_e^2 = kappa Lam^{-1} kappa;
      <dw_e dw_f> = (B Lam^{-1} B^T)_{ef} / (2 t_e t_f).
  Topological factor with insertion F (polynomial in the w_e):
      Tf[F; kappa] = U^{-1/2} e^{P} < F prod_r psi_r >,   psi_r = sum_e C_re w_e dt_e,   U = det(C t C^T).
  Returned as (Phi, U, P): Phi = <F prod_r psi_r> a form with rational coefficients.

Forms: dict {bitmask: Fraction}; bit e <-> dt_e; monomial = wedge of its dt_e in increasing e.
"""
from fractions import Fraction as Fr
from itertools import combinations
from functools import lru_cache
import random

# ---------------------------------------------------------------- linear algebra over Q

def det(M):
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
        piv = M[i][i]
        d *= piv
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / piv
                for c in range(i, n):
                    M[r][c] -= f * M[i][c]
    return d


def solve(A, Bm):
    n = len(A)
    k = len(Bm[0]) if n else 0
    aug = [[Fr(x) for x in A[i]] + [Fr(x) for x in Bm[i]] for i in range(n)]
    for i in range(n):
        p = next((r for r in range(i, n) if aug[r][i] != 0), None)
        if p is None:
            raise ZeroDivisionError("singular")
        aug[i], aug[p] = aug[p], aug[i]
        piv = aug[i][i]
        aug[i] = [x / piv for x in aug[i]]
        for r in range(n):
            if r != i and aug[r][i] != 0:
                f = aug[r][i]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[i])]
    return [row[n:] for row in aug]


def inverse(A):
    n = len(A)
    return solve(A, [[Fr(int(i == j)) for j in range(n)] for i in range(n)])


def matmul(A, Bm):
    Bt = list(zip(*Bm))
    return [[sum((a * b for a, b in zip(r, c)), Fr(0)) for c in Bt] for r in A]


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        while seq[i] != i:
            j = seq[i]
            seq[i], seq[j] = seq[j], seq[i]
            s = -s
    return s


# ---------------------------------------------------------------- exterior algebra

def _wsign(a, b):
    """sign of sorting (bits of a)(bits of b), a, b disjoint"""
    s = 0
    bb = b
    while bb:
        low = bb & -bb
        j = low.bit_length() - 1
        s += (a >> (j + 1)).bit_count()
        bb ^= low
    return -1 if s & 1 else 1


def wedge(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            if a & b:
                continue
            m = a | b
            v = out.get(m, 0) + _wsign(a, b) * x * y
            if v == 0:
                out.pop(m, None)
            else:
                out[m] = v
    return out


def wedge_all(fs):
    out = {0: Fr(1)}
    for f in fs:
        out = wedge(out, f)
    return out


def fadd(f, g, c=1):
    out = dict(f)
    for m, v in g.items():
        w = out.get(m, 0) + c * v
        if w == 0:
            out.pop(m, None)
        else:
            out[m] = w
    return out


def fscale(f, c):
    c = Fr(c)
    return {m: c * v for m, v in f.items()} if c != 0 else {}


def one_form(coeffs):
    return {1 << e: Fr(c) for e, c in enumerate(coeffs) if c != 0}


def pfaffian(A):
    """Pfaffian of an antisymmetric n x n matrix of even forms; Pf = sum_j (-1)^(j+1) A[0][j] Pf(minor)."""
    n = len(A)
    if n == 0:
        return {0: Fr(1)}
    if n % 2:
        return {}
    out = {}
    for j in range(1, n):
        if not A[0][j]:
            continue
        idx = [k for k in range(1, n) if k != j]
        sub = [[A[a][b] for b in idx] for a in idx]
        term = wedge(A[0][j], pfaffian(sub))
        out = fadd(out, term, 1 if j % 2 == 1 else -1)
    return out


def mask_of(S):
    m = 0
    for e in S:
        m |= 1 << e
    return m


def bits(m):
    out = []
    while m:
        low = m & -m
        out.append(low.bit_length() - 1)
        m ^= low
    return out


# ---------------------------------------------------------------- graphs

class Graph:
    def __init__(self, nv, edges, name=None):
        self.nv = nv
        self.edges = [tuple(e) for e in edges]
        self.name = name or "G"
        for a, b in self.edges:
            assert a != b and 0 <= a < nv and 0 <= b < nv

    @property
    def ne(self):
        return len(self.edges)

    @property
    def root(self):
        return self.nv - 1

    @property
    def loops(self):
        return self.ne - self.nv + 1

    def incidence(self):
        g = self.nv - 1
        B = [[0] * g for _ in range(self.ne)]
        for e, (a, b) in enumerate(self.edges):
            if a < g:
                B[e][a] += 1
            if b < g:
                B[e][b] -= 1
        return B

    def kruskal(self):
        par = list(range(self.nv))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        tree = []
        for e, (a, b) in enumerate(self.edges):
            ra, rb = find(a), find(b)
            if ra != rb:
                par[ra] = rb
                tree.append(e)
        assert len(tree) == self.nv - 1, "not connected"
        return tree

    def tree_path(self, tree, src, dst):
        """edges of the tree path src -> dst with signs (+1 along edge orientation)"""
        adj = {v: [] for v in range(self.nv)}
        for e in tree:
            a, b = self.edges[e]
            adj[a].append((b, e, +1))
            adj[b].append((a, e, -1))
        prev = {src: None}
        stack = [src]
        while stack:
            v = stack.pop()
            for (u, e, s) in adj[v]:
                if u not in prev:
                    prev[u] = (v, e, s)
                    stack.append(u)
        path = []
        v = dst
        while prev[v] is not None:
            p, e, s = prev[v]
            path.append((e, s))
            v = p
        return path[::-1]

    def fundamental_cycles(self, tree=None):
        tree = self.kruskal() if tree is None else list(tree)
        tset = set(tree)
        cot = [e for e in range(self.ne) if e not in tset]
        C = []
        for c in cot:
            a, b = self.edges[c]
            row = [0] * self.ne
            row[c] = 1
            for e, s in self.tree_path(tree, b, a):   # back from b to a through the tree
                row[e] += s
            C.append(row)
        return C, cot

    def cycle_matrix(self, tree=None):
        C, cot = self.fundamental_cycles(tree)
        if not C:
            return C
        Bt = list(map(list, zip(*self.incidence())))
        s = det(Bt + C)
        assert s != 0
        if s < 0:
            C[0] = [-x for x in C[0]]
        return C


# ---------------------------------------------------------------- circuit quantities

def laplacian(G, t):
    g = G.nv - 1
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(G.edges, t):
        c = 1 / Fr(te)
        if a < g:
            L[a][a] += c
        if b < g:
            L[b][b] += c
        if a < g and b < g:
            L[a][b] -= c
            L[b][a] -= c
    return L


def response(G, t, Linv=None):
    """R[e][v] = d wbar_e / d kappa_v (Kirchhoff currents), v < root"""
    g = G.nv - 1
    Linv = inverse(laplacian(G, t)) if Linv is None else Linv
    out = []
    for (a, b), te in zip(G.edges, t):
        ra = Linv[a] if a < g else [Fr(0)] * g
        rb = Linv[b] if b < g else [Fr(0)] * g
        out.append([(x - y) / Fr(te) for x, y in zip(ra, rb)])
    return out


def currents(G, t, src, Linv=None):
    """Kirchhoff currents for scalar sources src[v], v < root (or None)"""
    if src is None:
        return [Fr(0)] * G.ne
    R = response(G, t, Linv)
    return [sum((R[e][v] * Fr(src[v]) for v in range(G.nv - 1)), Fr(0)) for e in range(G.ne)]


def power(G, t, wbar):
    return sum((Fr(te) * w * w for te, w in zip(t, wbar)), Fr(0))


def period_matrix(G, t, C):
    return [[sum((Fr(ra[e] * rb[e]) * Fr(t[e]) for e in range(G.ne)), Fr(0)) for rb in C] for ra in C]


def w_covariance(G, t, Linv=None):
    """c[e][f] = <dw_e dw_f> = (B Lam^{-1} B^T)_{ef} / (2 t_e t_f)"""
    g = G.nv - 1
    Linv = inverse(laplacian(G, t)) if Linv is None else Linv
    B = G.incidence()
    BL = [[sum((B[e][a] * Linv[a][b] for a in range(g) if B[e][a]), Fr(0)) for b in range(g)] for e in range(G.ne)]
    c = [[sum((BL[e][b] * B[f][b] for b in range(g) if B[f][b]), Fr(0)) / (2 * Fr(t[e]) * Fr(t[f]))
          for f in range(G.ne)] for e in range(G.ne)]
    return c


def gaussian_moment(idx, mean, cov):
    """< prod_{k} X_{idx[k]} > for a Gaussian vector X with given mean and covariance (Isserlis)."""
    idx = tuple(sorted(idx))

    @lru_cache(maxsize=None)
    def mom(tup):
        if not tup:
            return Fr(1)
        first, rest = tup[0], tup[1:]
        tot = mean[first] * mom(rest)
        for k in range(len(rest)):
            c = cov[first][rest[k]]
            if c != 0:
                tot += c * mom(rest[:k] + rest[k + 1:])
        return tot
    return mom(idx)


def cotree_dets(C, ne):
    """{mask(S): det C[:, S]} over n-subsets S with nonzero det"""
    n = len(C)
    out = {}
    for S in combinations(range(ne), n):
        d = det([[C[r][e] for e in S] for r in range(n)])
        if d != 0:
            out[mask_of(S)] = d
    return out


class TopoData:
    """Precomputed data of one graph at one rational point t (one topological direction)."""

    def __init__(self, G, t, C=None):
        self.G = G
        self.t = [Fr(x) for x in t]
        self.C = G.cycle_matrix() if C is None else C
        self.n = len(self.C)
        self.Linv = inverse(laplacian(G, self.t)) if G.nv > 1 else []
        self.cov = w_covariance(G, self.t, self.Linv)
        self.M = period_matrix(G, self.t, self.C)
        self.U = det(self.M) if self.n else Fr(1)
        self._cd = None

    @property
    def cdets(self):
        if self._cd is None:
            self._cd = cotree_dets(self.C, self.G.ne)
        return self._cd

    def wbar(self, kappa):
        return currents(self.G, self.t, kappa, self.Linv)

    def topo_form(self, kappa=None, F=None):
        """(Phi, U, P):  Tf[F;kappa] = U^{-1/2} e^P Phi,  Phi = < F prod_r psi_r >.
        F: dict {tuple of edge indices (with repetition): coefficient} (polynomial in w), default 1."""
        wb = self.wbar(kappa)
        P = power(self.G, self.t, wb)
        F = {(): Fr(1)} if F is None else F
        Phi = {}
        cache = {}
        for m, d in self.cdets.items():
            S = bits(m)
            val = Fr(0)
            for mono, cf in F.items():
                key = tuple(sorted(tuple(mono) + tuple(S)))
                if key not in cache:
                    cache[key] = gaussian_moment(key, wb, self.cov)
                val += Fr(cf) * cache[key]
            if val != 0:
                Phi[m] = d * val
        return Phi, self.U, P

    def topo_form_cut(self, kappa=None, F=None, A=()):
        """(Phi, U, P) with Phi = < F d_A prod_r psi_r >,  d_A = odd derivative with respect to the
        one-forms zeta_e = w_e dt_e, e in A:  prod psi = sum_A zeta_{a1}^...^zeta_{am} ^ d_A(prod psi) + ..."""
        wb = self.wbar(kappa)
        P = power(self.G, self.t, wb)
        F = {(): Fr(1)} if F is None else F
        Am = mask_of(A)
        Phi = {}
        cache = {}
        for m, d in self.cdets.items():
            if m & Am != Am:
                continue
            R = m ^ Am
            S = bits(R)
            val = Fr(0)
            for mono, cf in F.items():
                key = tuple(sorted(tuple(mono) + tuple(S)))
                if key not in cache:
                    cache[key] = gaussian_moment(key, wb, self.cov)
                val += Fr(cf) * cache[key]
            if val != 0:
                Phi[R] = Phi.get(R, 0) + d * _wsign(Am, R) * val
        return {k: v for k, v in Phi.items() if v != 0}, self.U, P

    def rho_pfaffian(self):
        """Pf(-1/2 dM M^{-1} dM) as a form (rho_n = U^{-1/2} * this)."""
        n, C, ne = self.n, self.C, self.G.ne
        if n == 0:
            return {0: Fr(1)}
        Minv = inverse(self.M)
        dM = [[one_form([C[r][e] * C[a][e] for e in range(ne)]) for a in range(n)] for r in range(n)]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                acc = {}
                for a in range(n):
                    for b in range(n):
                        if Minv[a][b] != 0:
                            acc = fadd(acc, wedge(dM[r][a], dM[s][b]), -Minv[a][b] / 2)
                A[r][s] = acc
                A[s][r] = fscale(acc, -1)
        return pfaffian(A)

    def transfer_K(self):
        """K = C^T M^{-1} C  (transfer admittance)"""
        Minv = inverse(self.M)
        C = self.C
        n, ne = self.n, self.G.ne
        CM = [[sum((C[a][e] * Minv[a][b] for a in range(n)), Fr(0)) for b in range(n)] for e in range(ne)]
        return [[sum((CM[e][b] * C[b][f] for b in range(n)), Fr(0)) for f in range(ne)] for e in range(ne)]

    def rho_hafnian(self):
        """(-2)^{-n/2} sum_S det(C_S) Haf(K_SS) dt_S  (rho_n = U^{-1/2} * this)"""
        n = self.n
        if n % 2:
            return {}
        K = self.transfer_K()
        out = {}
        for m, d in self.cdets.items():
            S = bits(m)
            h = hafnian([[K[a][b] for b in S] for a in S])
            if h != 0:
                out[m] = d * h * Fr(-2) ** (-(n // 2))
        return out

    def gkw_fiber(self, kappa=None):
        """Spanning-tree expansion of the fiber integral of exp(kappa.tau) prod_e Theta_e over the vertex
        potentials (prod_v d tau_v in front, vertex order).  Returns Phi_GKW with
        fiber integral = (2 sqrt(pi))^{-n} U^{-1/2} e^{P} Phi_GKW."""
        G = self.G
        g = G.nv - 1
        B = G.incidence()
        wb = self.wbar(kappa)
        out = {}
        for T in combinations(range(G.ne), g):
            dB = det([B[e] for e in T])
            if dB == 0:
                continue
            Tm = mask_of(T)
            cot = [e for e in range(G.ne) if not (Tm >> e) & 1]
            # sign of moving the tree factors (d tau_e, e in T) to the front, keeping relative orders
            order = list(T) + cot
            eps = perm_sign([order.index(e) for e in range(G.ne)])
            # wedge_{e in T} d tau_e = det(B_T) prod_v d tau_v ;  cotree factors (-w_e dt_e)
            val = gaussian_moment(tuple(cot), wb, self.cov) * (-1) ** len(cot)
            if val != 0:
                m = mask_of(cot)
                out[m] = out.get(m, 0) + eps * dB * val
        return {m: v for m, v in out.items() if v != 0}


def hafnian(A):
    n = len(A)
    if n == 0:
        return Fr(1)
    if n % 2:
        return Fr(0)
    tot = Fr(0)
    for j in range(1, n):
        if A[0][j] == 0:
            continue
        idx = [k for k in range(1, n) if k != j]
        tot += A[0][j] * hafnian([[A[a][b] for b in idx] for a in idx])
    return tot


def rand_t(ne, rng, lo=1, hi=40):
    return [Fr(rng.randint(lo, hi), rng.randint(1, 12)) for _ in range(ne)]


def det_G_sign(G, C):
    """det G, G = (P; C) with P B = 1, P = B_T^{-1} on the Kruskal tree T (CONVENTIONS sec 4)."""
    T = G.kruskal()
    B = G.incidence()
    BT = [B[e] for e in T]          # g x g
    BTinv = inverse(BT)             # g x g ; P = BTinv on tree columns, 0 elsewhere
    g = G.nv - 1
    Pm = [[Fr(0)] * G.ne for _ in range(g)]
    for i in range(g):
        for k, e in enumerate(T):
            Pm[i][e] = BTinv[i][k]
    return det(Pm + [[Fr(x) for x in r] for r in C])


def pullback(form, images):
    """images[e] = one-form (dict) replacing dt_e; returns the pulled-back form"""
    out = {}
    for m, c in form.items():
        f = {0: Fr(c)}
        for e in bits(m):
            f = wedge(f, images[e])
            if not f:
                break
        out = fadd(out, f)
    return out
