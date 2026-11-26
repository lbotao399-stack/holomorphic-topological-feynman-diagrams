"""Independent brute-force layer (does not import htk).

Graph: nv vertices 0..nv-1, root = nv-1, edges[e] = (e0, e1).
Kirchhoff:  t_e y_e = phi_{e0} - phi_{e1}, phi_root = 0;  sum_{e0=v} y_e - sum_{e1=v} y_e = lam_v.
Cycle matrix: fundamental cycles of a BFS tree, first row multiplied by -1 if det[B^T; C] < 0.
Form:  Omega = U^{-H} rho_n(M)^T ^ prod_{i=1..H} prod_{r=1..n} psi^i_r,
       rho_n(M) = det(M)^{-1/2} Pf(-1/2 dM M^{-1} dM),
       Omega = D(t) iota_E dt_0...dt_{N-1},  D(t) = [Omega]_{dt_1...dt_{N-1}} / t_0.
D(t) is returned as (q, U) with D = q (T = 0) or D = q / sqrt(U) (T = 1).
"""
from fractions import Fraction as Fr
from itertools import combinations


def mdet(A):
    A = [[Fr(x) for x in row] for row in A]
    n = len(A)
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return d


def minv(A):
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def incidence(nv, edges):
    return [[(1 if v == a else 0) - (1 if v == b else 0) for v in range(nv - 1)] for (a, b) in edges]


def cycle_matrix(nv, edges):
    N = len(edges)
    adj = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        adj[a].append((b, e, +1))
        adj[b].append((a, e, -1))
    root = nv - 1
    parent = {root: None}
    order = [root]
    k = 0
    while k < len(order):
        u = order[k]
        k += 1
        for (w, e, s) in sorted(adj[u], key=lambda x: x[1]):
            if w not in parent:
                parent[w] = (u, e, s)
                order.append(w)
    tree = set(p[1] for p in parent.values() if p is not None)

    def path_to_root(v):
        out = []
        while parent[v] is not None:
            u, e, s = parent[v]
            out.append((e, s))  # traversing u->v has sign s; v->u has sign -s
            v = u
        return out

    C = []
    for e, (a, b) in enumerate(edges):
        if e in tree:
            continue
        row = [0] * N
        row[e] += 1
        # go a -> b along e, then back b -> root -> a
        for (f, s) in path_to_root(b):
            row[f] += -s
        for (f, s) in path_to_root(a):
            row[f] += s
        C.append(row)
    B = incidence(nv, edges)
    BT = [[B[e][v] for e in range(N)] for v in range(nv - 1)]
    if C and mdet(BT + C) < 0:
        C[0] = [-x for x in C[0]]
    return C


def currents(nv, edges, t, lam):
    """lam: dict v -> H-vector (non-root vertices), missing = 0."""
    g = nv - 1
    H = len(next(iter(lam.values())))
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
    Li = minv(L)
    lamv = [[Fr(x) for x in lam.get(v, [0] * H)] for v in range(g)]
    phi = [[sum(Li[v][u] * lamv[u][i] for u in range(g)) for i in range(H)] for v in range(g)] + [[Fr(0)] * H]
    y = [[(phi[a][i] - phi[b][i]) / Fr(te) for i in range(H)] for (a, b), te in zip(edges, t)]
    return y


# forms: dict sorted-tuple -> Fraction
def wedge(f, g):
    out = {}
    for I, a in f.items():
        for J, b in g.items():
            if set(I) & set(J):
                continue
            K = I + J
            # sign of sorting permutation
            inv = sum(1 for i in range(len(K)) for j in range(i + 1, len(K)) if K[i] > K[j])
            KS = tuple(sorted(K))
            out[KS] = out.get(KS, Fr(0)) + (-a * b if inv % 2 else a * b)
    return {k: v for k, v in out.items() if v != 0}


def one_form(coeffs):
    return {(e,): Fr(c) for e, c in enumerate(coeffs) if c != 0}


def fadd(f, g, s=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, Fr(0)) + s * v
    return {k: v for k, v in out.items() if v != 0}


def pfaffian(A):
    n = len(A)
    if n == 0:
        return {(): Fr(1)}
    if n % 2:
        return {}
    out = {}
    for j in range(1, n):
        idx = [k for k in range(1, n) if k != j]
        sub = [[A[a][b] for b in idx] for a in idx]
        term = wedge(A[0][j], pfaffian(sub))
        out = fadd(out, term, -1 if (j - 1) % 2 else 1)
    return out


def density(nv, edges, t, lam, T=0, C=None):
    t = [Fr(x) for x in t]
    N = len(edges)
    C = cycle_matrix(nv, edges) if C is None else C
    n = len(C)
    y = currents(nv, edges, t, lam)
    H = len(y[0])
    M = [[sum(Fr(C[r][e] * C[s][e]) * t[e] for e in range(N)) for s in range(n)] for r in range(n)]
    U = mdet(M)
    form = {(): Fr(1)}
    if T == 1:
        Mi = minv(M)
        dM = [[one_form([C[r][e] * C[s][e] for e in range(N)]) for s in range(n)] for r in range(n)]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(n):
                acc = {}
                for a in range(n):
                    for b in range(n):
                        if Mi[a][b] != 0:
                            acc = fadd(acc, {k: v * (-Mi[a][b] / 2) for k, v in wedge(dM[r][a], dM[b][s]).items()})
                A[r][s] = acc
        form = pfaffian(A)
    for i in range(H):
        for r in range(n):
            form = wedge(form, one_form([C[r][e] * y[e][i] for e in range(N)]))
    key = tuple(range(1, N))
    q = form.get(key, Fr(0)) / t[0] / U ** H
    return q, U


def response_and_derivs(nv, edges, t):
    """R[e][v] (y_e = sum_v R[e][v] lam_v) and dR[j][e][v] = d R[e][v] / d t_j, exact."""
    g = nv - 1
    N = len(edges)
    t = [Fr(x) for x in t]
    B = incidence(nv, edges)
    L = [[sum(Fr(B[e][a] * B[e][b]) / t[e] for e in range(N)) for b in range(g)] for a in range(g)]
    Li = minv(L)
    # R = t^{-1} B L^{-1}
    BLi = [[sum(B[e][u] * Li[u][v] for u in range(g)) for v in range(g)] for e in range(N)]
    R = [[BLi[e][v] / t[e] for v in range(g)] for e in range(N)]
    dR = []
    for j in range(N):
        # d(t^{-1})_jj = -1/t_j^2 ; dL = B^T d(t^{-1}) B = -(1/t_j^2) b_j b_j^T
        # d L^{-1} = -L^{-1} dL L^{-1} = (1/t_j^2) (L^{-1} b_j)(b_j^T L^{-1})
        Lb = [sum(Li[v][u] * B[j][u] for u in range(g)) for v in range(g)]
        dLi = [[Lb[a] * Lb[b] / t[j] ** 2 for b in range(g)] for a in range(g)]
        dRj = []
        for e in range(N):
            row = []
            for v in range(g):
                val = sum(B[e][u] * dLi[u][v] for u in range(g)) / t[e]
                if e == j:
                    val -= BLi[e][v] / t[e] ** 2
                row.append(val)
            dRj.append(row)
        dR.append(dRj)
    return R, dR
