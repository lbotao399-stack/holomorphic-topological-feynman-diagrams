"""Independent kernel layer for VERIFY.md (imports neither htk nor the explorer's code).

Graph: nv vertices 0..nv-1, root = nv-1, edges[e] = (e0, e1).
B[e][v] = d(v,e0) - d(v,e1) (v < root).  Kirchhoff: t_e y_e = phi_e0 - phi_e1, phi_root = 0,
B^T y = lam.  Cycle matrix: fundamental cycles of a spanning tree chosen by the caller (default:
DFS tree from vertex 0), first row negated if det[B^T; C] < 0.
Omega = U^{-H} rho^T ^ prod_{i=1..H} prod_{r=1..n} psi^i_r,  psi^i_r = sum_e C_re y^i_e dt_e,
rho_2 = U^{-1/2} A_12, A = -1/2 dM M^{-1} dM  (T = 1, n = 2 only).
Omega = D(t) iota_E (dt_0 ... dt_{N-1}),  D = (-1)^{e0} [Omega]_{dt_0..^e0..dt_{N-1}} / t_e0.
Returned as (q, U): D = q (T = 0), D = q U^{-1/2} (T = 1).
Derivatives d/dt_j of currents: exact forward mode, L(t + eps e_j) = L0 + eps L1.
"""
from fractions import Fraction as Fr
import itertools, random


def det(A):
    A = [[Fr(x) for x in r] for r in A]
    n = len(A)
    if n == 0:
        return Fr(1)
    s = Fr(1)
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
            s = -s
        piv = A[c][c]
        s *= piv
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / piv
                Ar, Ac = A[r], A[c]
                for k in range(c, n):
                    Ar[k] -= f * Ac[k]
    return s


def solve(A, b):
    """A x = b, b a list of column vectors (list of lists, each length n)."""
    n = len(A)
    m = len(b)
    M = [[Fr(x) for x in A[i]] + [Fr(b[j][i]) for j in range(m)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                Mr, Mc = M[r], M[c]
                M[r] = [x - f * y for x, y in zip(Mr, Mc)]
    return [[M[i][n + j] for i in range(n)] for j in range(m)]


def inv(A):
    n = len(A)
    cols = solve(A, [[int(i == j) for i in range(n)] for j in range(n)])
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def incidence(nv, edges):
    g = nv - 1
    return [[(1 if v == a else 0) - (1 if v == b else 0) for v in range(g)] for (a, b) in edges]


def dfs_tree(nv, edges, start=0, order=None):
    """spanning tree by DFS from `start`, neighbours visited in `order` of edge ids (default reversed)."""
    adj = {v: [] for v in range(nv)}
    eids = list(range(len(edges))) if order is None else order
    for e in eids:
        a, b = edges[e]
        adj[a].append((b, e))
        adj[b].append((a, e))
    seen = {start}
    tree = []
    stack = [start]
    while stack:
        u = stack.pop()
        for (w, e) in adj[u]:
            if w not in seen:
                seen.add(w)
                tree.append(e)
                stack.append(w)
    assert len(tree) == nv - 1
    return sorted(tree)


def random_tree(nv, edges, rng):
    order = list(range(len(edges)))
    rng.shuffle(order)
    return dfs_tree(nv, edges, start=rng.randrange(nv), order=order)


def tree_path(nv, edges, tree, a, b):
    """signed edge list of the tree path from a to b: list of (e, s), s = +1 if traversed e0->e1."""
    adj = {v: [] for v in range(nv)}
    for e in tree:
        x, y = edges[e]
        adj[x].append((y, e, +1))
        adj[y].append((x, e, -1))
    prev = {a: None}
    q = [a]
    while q:
        u = q.pop()
        for (w, e, s) in adj[u]:
            if w not in prev:
                prev[w] = (u, e, s)
                q.append(w)
    out = []
    v = b
    while prev[v] is not None:
        u, e, s = prev[v]
        out.append((e, s))
        v = u
    return out[::-1]


def fundamental_cycles(nv, edges, tree):
    N = len(edges)
    C = []
    cot = [e for e in range(N) if e not in tree]
    for c in cot:
        a, b = edges[c]
        row = [0] * N
        row[c] = 1
        for (e, s) in tree_path(nv, edges, tree, b, a):  # back from head to tail
            row[e] += s
        C.append(row)
    return C, cot


def sigma_tree(nv, edges, tree):
    """sigma_T = eps_T det B_T."""
    N = len(edges)
    cot = [e for e in range(N) if e not in tree]
    perm = list(tree) + cot
    inv_ = sum(1 for i in range(N) for j in range(i + 1, N) if perm[i] > perm[j])
    B = incidence(nv, edges)
    return (-1) ** inv_ * int(det([B[e] for e in tree]))


def positive_cycles(nv, edges, tree):
    C, cot = fundamental_cycles(nv, edges, tree)
    B = incidence(nv, edges)
    N = len(edges)
    BT = [[B[e][v] for e in range(N)] for v in range(nv - 1)]
    dd = det(BT + C)
    assert dd != 0
    if dd < 0:
        C = [[-x for x in C[0]]] + C[1:]
    return C


def check_cycles(nv, edges, C):
    """C B = 0 (cycles are cycles)."""
    B = incidence(nv, edges)
    return all(sum(C[r][e] * B[e][v] for e in range(len(edges))) == 0 for r in range(len(C)) for v in range(nv - 1))


def laplacian(nv, edges, t):
    g = nv - 1
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(edges, t):
        c = 1 / Fr(te)
        if a < g:
            L[a][a] += c
        if b < g:
            L[b][b] += c
        if a < g and b < g:
            L[a][b] -= c
            L[b][a] -= c
    return L


def currents(nv, edges, t, lam, H):
    """lam: dict v -> H-vector for v < root (missing = 0). Returns y[e][i]."""
    g = nv - 1
    L = laplacian(nv, edges, t)
    rhs = [[Fr(lam.get(v, [0] * H)[i]) for v in range(g)] for i in range(H)]
    phis = solve(L, rhs)  # phis[i][v]
    y = []
    for (a, b), te in zip(edges, t):
        pa = [phis[i][a] if a < g else Fr(0) for i in range(H)]
        pb = [phis[i][b] if b < g else Fr(0) for i in range(H)]
        y.append([(pa[i] - pb[i]) / Fr(te) for i in range(H)])
    return y


def current_derivs(nv, edges, t, lam, H):
    """dy[j][e][i] = d y^i_e / d t_j, exact forward mode."""
    g = nv - 1
    t = [Fr(x) for x in t]
    L0 = laplacian(nv, edges, t)
    rhs = [[Fr(lam.get(v, [0] * H)[i]) for v in range(g)] for i in range(H)]
    phi0 = solve(L0, rhs)
    out = []
    for j, (a, b) in enumerate(edges):
        # d(1/t_j) = -1/t_j^2 ; L1 = -(1/t_j^2) b_j b_j^T
        c = -1 / t[j] ** 2
        bj = [(1 if v == a else 0) - (1 if v == b else 0) for v in range(g)]
        rhs1 = []
        for i in range(H):
            L1phi = [c * bj[v] * sum(bj[u] * phi0[i][u] for u in range(g)) for v in range(g)]
            rhs1.append([-x for x in L1phi])
        phi1 = solve(L0, rhs1)
        dyj = []
        for e, (x, y_) in enumerate(edges):
            row = []
            for i in range(H):
                p0 = (phi0[i][x] if x < g else 0) - (phi0[i][y_] if y_ < g else 0)
                p1 = (phi1[i][x] if x < g else 0) - (phi1[i][y_] if y_ < g else 0)
                val = p1 / t[e]
                if e == j:
                    val -= p0 / t[e] ** 2
                row.append(val)
            dyj.append(row)
        out.append(dyj)
    return out


def period(C, t):
    n = len(C)
    N = len(t)
    return [[sum(Fr(C[r][e] * C[s][e]) * t[e] for e in range(N)) for s in range(n)] for r in range(n)]


def density(nv, edges, t, lam, H, T=0, C=None, e0=0):
    """(q, U) with D = q (T=0) or q U^{-1/2} (T=1)."""
    t = [Fr(x) for x in t]
    N = len(edges)
    if C is None:
        C = positive_cycles(nv, edges, dfs_tree(nv, edges))
    n = len(C)
    y = currents(nv, edges, t, lam, H)
    M = period(C, t)
    U = det(M)
    Psi = [[C[r][e] * y[e][i] for e in range(N)] for i in range(H) for r in range(n)]
    S = [e for e in range(N) if e != e0]
    if T == 0:
        assert len(Psi) == N - 1
        cf = det([[row[e] for e in S] for row in Psi])
        return (-1) ** e0 * cf / t[e0] / U ** H, U
    assert T == 1 and n == 2 and len(Psi) == N - 3
    Mi = inv(M)
    # A_12 = -1/2 sum_ab Mi[a][b] dM_{1a} ^ dM_{b2};  dM_{rs} = sum_e C_re C_se dt_e
    rho = {}
    for a in range(2):
        for b in range(2):
            if Mi[a][b] == 0:
                continue
            u = [C[0][e] * C[a][e] for e in range(N)]
            v = [C[b][e] * C[1][e] for e in range(N)]
            for e in range(N):
                for f in range(e + 1, N):
                    val = u[e] * v[f] - u[f] * v[e]
                    if val:
                        rho[(e, f)] = rho.get((e, f), Fr(0)) - Fr(1, 2) * Mi[a][b] * val
    cf = Fr(0)
    pos = {e: k for k, e in enumerate(S)}
    for (e, f), r in rho.items():
        if e == e0 or f == e0 or r == 0:
            continue
        rest = [x for x in S if x not in (e, f)]
        sgn = (-1) ** (pos[e] + pos[f] - 1)
        cf += sgn * r * det([[row[x] for x in rest] for row in Psi])
    return (-1) ** e0 * cf / t[e0] / U ** H, U


def dy_density(nv, edges, t, lam, H, tree):
    """density of wedge_{r} d^H y_{c_r} (cotree of `tree` in increasing order), w.r.t. iota_E vol,
    using exact forward-mode derivatives:  [form]_{dt_1..dt_{N-1}} / t_0."""
    N = len(edges)
    cot = [e for e in range(N) if e not in tree]
    dy = current_derivs(nv, edges, t, lam, H)
    rows = [[dy[j][c][i] for j in range(1, N)] for c in cot for i in range(H)]
    return det(rows) / Fr(t[0])


def eps_nH(n, H):
    return (-1) ** (n * H + (n * (n - 1) // 2) * (H * (H - 1) // 2))


# ---------------------------------------------------------------- Grundmann-Moeller rules
from math import factorial


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for k in range(total + 1):
        for rest in compositions(total - k, parts - 1):
            yield (k,) + rest


def gm_rule(nd, s):
    """Grundmann-Moeller rule of index s on the nd-simplex {x_i >= 0, sum x_i <= 1}, exact for
    polynomials of degree <= 2s+1.  Returns list of (weight, point) with exact Fractions; weights
    sum to the volume 1/nd!.  Point = barycentric coordinates (nd+1 of them)."""
    d = 2 * s + 1
    out = []
    for i in range(s + 1):
        w = Fr((-1) ** i * (d + nd - 2 * i) ** d, 2 ** (2 * s) * factorial(i) * factorial(d + nd - i))
        for beta in compositions(s - i, nd + 1):
            pt = tuple(Fr(2 * b + 1, d + nd - 2 * i) for b in beta)
            out.append((w, pt))
    # merge duplicates
    acc = {}
    for w, p in out:
        acc[p] = acc.get(p, Fr(0)) + w
    return [(w, p) for p, w in acc.items() if w != 0]


def dirichlet(exps):
    """int over the simplex (barycentric b_0..b_m, Lebesgue measure on m free coords) of prod b_j^a_j
    for nonnegative integers a_j."""
    m = len(exps) - 1
    num = 1
    for a in exps:
        num *= factorial(a)
    return Fr(num, factorial(m + sum(exps)))
