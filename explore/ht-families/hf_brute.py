"""Independent exact brute-force layer for explore/ht-families (does not import htk).

Graph: nv vertices 0..nv-1, root = nv-1, edges = list of (tail, head).
Conventions of code/htk/CONVENTIONS.md Sections 1-2:
  B_ev = delta_{v,e0} - delta_{v,e1} (v < root);  t y = B phi,  B^T y = lambda;
  C: cycle matrix with det [B^T ; C] > 0;  M = C t C^T, U = det M;
  psi^i_r = sum_e C_re y^i_e dt_e;  rho_n = U^{-1/2} Pf(-1/2 dM M^{-1} dM);
  U^{-H} rho_n^T ^ prod_{i=1..H} prod_{r=1..n} psi^i_r = D(t) iota_E (dt_0 ... dt_{N-1}),
  D = (-1)^{N-1} [Omega]_{dt_0..dt_{N-2}} / t_{N-1}.
Forms: dict {sorted tuple of generator indices: Fraction}.
"""
from fractions import Fraction as Fr
import itertools


# ------------------------------------------------------------------ exact linear algebra
def fdet(A):
    A = [[Fr(x) for x in row] for row in A]
    n = len(A)
    if n == 0:
        return Fr(1)
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] * inv
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return d


def finv(A):
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


# ------------------------------------------------------------------ exterior algebra
def merge_sign(a, b):
    """sign of sorting the concatenation a + b (a, b sorted, disjoint); 0 if they overlap"""
    if set(a) & set(b):
        return 0
    inv = 0
    j = 0
    for x in a:
        inv += sum(1 for y in b if y < x)
    return -1 if inv & 1 else 1


def wedge(f, g):
    out = {}
    for ka, va in f.items():
        for kb, vb in g.items():
            s = merge_sign(ka, kb)
            if s == 0:
                continue
            key = tuple(sorted(ka + kb))
            out[key] = out.get(key, 0) + s * va * vb
    return {k: v for k, v in out.items() if v != 0}


def fadd(f, g, c=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def pfaff(A, idx):
    """Pfaffian of the antisymmetric matrix A of 2-forms (A[r][s] dicts), rows idx."""
    if not idx:
        return {(): Fr(1)}
    if len(idx) % 2:
        return {}
    out = {}
    i0 = idx[0]
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        term = wedge(A[i0][idx[j]], pfaff(A, rest))
        out = fadd(out, term, -1 if (j - 1) % 2 else 1)
    return out


# ------------------------------------------------------------------ graph data
class G:
    def __init__(self, nv, edges):
        self.nv, self.edges = nv, [tuple(e) for e in edges]
        self.ne = len(edges)
        self.root = nv - 1
        self.n = self.ne - nv + 1

    def incidence(self):
        g = self.nv - 1
        B = [[0] * g for _ in self.edges]
        for e, (a, b) in enumerate(self.edges):
            if a < g:
                B[e][a] += 1
            if b < g:
                B[e][b] -= 1
        return B

    def bfs_tree(self):
        """spanning tree by BFS from the root (different from the engine's Kruskal tree)"""
        adj = {v: [] for v in range(self.nv)}
        for e, (a, b) in enumerate(self.edges):
            adj[a].append((e, b))
            adj[b].append((e, a))
        seen, tree, queue = {self.root}, [], [self.root]
        while queue:
            v = queue.pop(0)
            for e, u in sorted(adj[v]):
                if u not in seen:
                    seen.add(u)
                    tree.append(e)
                    queue.append(u)
        return sorted(tree)

    def cycles(self, tree=None):
        """fundamental cycles of the tree, first row flipped so that det [B^T ; C] > 0"""
        tree = self.bfs_tree() if tree is None else tree
        tset = set(tree)
        adj = {v: [] for v in range(self.nv)}
        for e in tree:
            a, b = self.edges[e]
            adj[a].append((e, b, +1))
            adj[b].append((e, a, -1))

        def path(u, v):  # signed edge path in the tree from u to v
            prev = {u: None}
            stack = [u]
            while stack:
                x = stack.pop()
                for e, y, s in adj[x]:
                    if y not in prev:
                        prev[y] = (x, e, s)
                        stack.append(y)
            out = {}
            x = v
            while prev[x] is not None:
                px, e, s = prev[x]
                out[e] = out.get(e, 0) + s
                x = px
            return out

        C = []
        for c in range(self.ne):
            if c in tset:
                continue
            a, b = self.edges[c]
            row = [0] * self.ne
            row[c] = 1
            for e, s in path(b, a).items():   # back from head to tail through the tree
                row[e] += s
            C.append(row)
        B = self.incidence()
        Mt = [[B[e][v] for e in range(self.ne)] for v in range(self.nv - 1)] + C
        d = fdet(Mt)
        assert d != 0
        if d < 0:
            C[0] = [-x for x in C[0]]
        return C


def currents(Gr, t, lam):
    """y[e][i] by direct solution of Kirchhoff's laws (exact)"""
    g = Gr.nv - 1
    H = len(lam[0])
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(Gr.edges, t):
        c = 1 / Fr(te)
        if a < g:
            L[a][a] += c
        if b < g:
            L[b][b] += c
        if a < g and b < g:
            L[a][b] -= c
            L[b][a] -= c
    Li = finv(L)
    phi = [[sum(Li[v][u] * Fr(lam[u][i]) for u in range(g)) for i in range(H)] for v in range(g)] + [[Fr(0)] * H]
    return [[(phi[a][i] - phi[b][i]) / Fr(te) for i in range(H)] for (a, b), te in zip(Gr.edges, t)]


def check_kcl(Gr, t, lam, y, C):
    B = Gr.incidence()
    H = len(lam[0])
    for v in range(Gr.nv - 1):
        for i in range(H):
            assert sum(B[e][v] * y[e][i] for e in range(Gr.ne)) == Fr(lam[v][i])
    for row in C:
        for i in range(H):
            assert sum(row[e] * Fr(t[e]) * y[e][i] for e in range(Gr.ne)) == 0


def period(Gr, t, C):
    return [[sum(Fr(a[e] * b[e]) * Fr(t[e]) for e in range(Gr.ne)) for b in C] for a in C]


def rho_form(Gr, t, C, gens=None):
    """Pf(-1/2 dM M^{-1} dM) as an n-form in the dt_e (e in gens; default all edges).
    rho_n = U^{-1/2} * this."""
    n = len(C)
    gens = list(range(Gr.ne)) if gens is None else gens
    M = period(Gr, t, C)
    Mi = finv(M)
    dM1 = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for a in range(n):
            f = {}
            for e in gens:
                c = C[r][e] * C[a][e]
                if c:
                    f[(e,)] = Fr(c)
            dM1[r][a] = f
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(r + 1, n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Mi[a][b] != 0:
                        acc = fadd(acc, wedge(dM1[r][a], dM1[s][b]), -Mi[a][b] / 2)
            A[r][s] = acc
            A[s][r] = {k: -v for k, v in acc.items()}
    return pfaff(A, list(range(n))), fdet(M)


def density(Gr, t, lam, T=1, C=None):
    """(coeff, U): D(t) = coeff * U^{-T/2}.  T in {0,1}."""
    t = [Fr(x) for x in t]
    C = Gr.cycles() if C is None else C
    n, N = len(C), Gr.ne
    H = len(lam[0]) if lam else 0
    assert n * (H + T) == N - 1
    y = currents(Gr, t, lam) if H else [[] for _ in range(N)]
    if H:
        check_kcl(Gr, t, lam, y, C)
    gens = list(range(N - 1))                      # drop dt_{N-1}
    U = fdet(period(Gr, t, C))
    rows = [[C[r][e] * y[e][i] for e in gens] for i in range(H) for r in range(n)]
    if T == 1:
        rho, _ = rho_form(Gr, t, C, gens)
    else:
        rho = {(): Fr(1)}
    top = Fr(0)
    full = tuple(gens)
    for S, val in rho.items():
        rest = [e for e in gens if e not in S]
        s = merge_sign(S, tuple(rest))
        if s == 0:
            continue
        minor = fdet([[row[gens.index(e)] for e in rest] for row in rows]) if rows else Fr(1)
        top += s * val * minor
    coeff = (-1) ** (N - 1) * top / (t[N - 1] * U ** H)
    return coeff, U
