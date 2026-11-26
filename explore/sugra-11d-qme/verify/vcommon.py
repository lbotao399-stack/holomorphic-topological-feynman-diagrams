"""Common helpers for the independent verification (no import of htk or of the explorer's code).

Theta graph theta_pqr (S3 Section sec-Q2): vertices N = 0, interior vertices of A, B, C numbered
along the paths, N' = 11 (root, last).  Paths A, B, C from N to N', edges oriented N -> N',
edge order (A_1..A_p, B_1..B_q, C_1..C_r).
"""
from fractions import Fraction as Fr
import random

PQR = [(1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5)]
H = 5


def theta_graph(p, q, r):
    nv = p + q + r - 1
    root = nv - 1
    edges, paths, nxt = [], [], 1
    for n_ in (p, q, r):
        verts = [0]
        for j in range(n_ - 1):
            verts.append(nxt)
            nxt += 1
        verts.append(root)
        pe = []
        for j in range(n_):
            pe.append(len(edges))
            edges.append((verts[j], verts[j + 1]))
        paths.append((verts[1:-1], pe))
    return nv, edges, paths


def gauss_det(M):
    M = [[Fr(x) for x in row] for row in M]
    n = len(M)
    d = Fr(1)
    for i in range(n):
        piv = None
        for r in range(i, n):
            if M[r][i] != 0:
                piv = r
                break
        if piv is None:
            return Fr(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / M[i][i]
                for c in range(i, n):
                    M[r][c] -= f * M[i][c]
    return d


def gauss_solve(A, b):
    """Solve A x = b (b: list of vectors -> returns list of vectors), exact."""
    n = len(A)
    m = len(b[0])
    M = [[Fr(x) for x in A[i]] + [Fr(x) for x in b[i]] for i in range(n)]
    for i in range(n):
        piv = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[piv] = M[piv], M[i]
        pv = M[i][i]
        M[i] = [x / pv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [row[n:] for row in M]


def kirchhoff(nv, edges, t, lam):
    """Own Kirchhoff solver.  lam: dict v -> H-vector for v != root (root = nv-1 absorbs minus the sum).
    Potentials phi (phi_root = 0), t_e y_e = phi_{e(0)} - phi_{e(1)}, sum_{e(0)=v} y_e - sum_{e(1)=v} y_e = lam_v."""
    g = nv - 1
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(edges, t):
        c = 1 / Fr(te)
        for u in (a, b):
            if u < g:
                L[u][u] += c
        if a < g and b < g:
            L[a][b] -= c
            L[b][a] -= c
    rhs = [list(lam.get(v, [0] * H)) for v in range(g)]
    phi = gauss_solve(L, rhs) + [[Fr(0)] * H]
    return [[(phi[a][i] - phi[b][i]) / Fr(te) for i in range(H)] for (a, b), te in zip(edges, t)]


def D_theta(lam, paths):
    """S3 eq-Dtheta: det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]]."""
    A, B, C = [pv for pv, pe in paths]
    rows = []
    for i in range(H):
        rows.append([lam[v][i] for v in A] + [lam[v][i] for v in B] + [Fr(0)] * len(C))
    for i in range(H):
        rows.append([Fr(0)] * len(A) + [lam[v][i] for v in B] + [lam[v][i] for v in C])
    return gauss_det(rows)


def perm_parity_sign(seq):
    """Sign of the permutation that sorts seq (distinct entries), by inversion count."""
    inv = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                inv += 1
    return -1 if inv % 2 else 1


def rq(rng, a=1, b=11):
    return Fr(rng.randint(-b, b), rng.randint(a, b))


def rsimplex(rng, n):
    a = [Fr(rng.randint(1, 13)) for _ in range(n)]
    s = sum(a)
    return [x / s for x in a]
