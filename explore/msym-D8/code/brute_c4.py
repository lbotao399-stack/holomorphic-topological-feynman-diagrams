"""Independent brute-force layer for T = 0 kernels (written for this task; does not import htk).

Network: Ohm t_e y_e = phi_{e(0)} - phi_{e(1)}, phi_root = 0; Kirchhoff sum_{e(0)=v} y_e - sum_{e(1)=v} y_e
= lambda_v (v != root).  Kernel (S2 eq-kernel-def with d^4 y, H = 4; sign factor (eps_T det B_T)^4 = 1):

    I_Gamma = int_{RP^{m-1}_{>0}, omega_m} exp(-sum_e y_e . z_e) wedge_{c in cotree} d^4 y_c,
    omega_m = sum_j (-1)^{j-1} t_j dt_1 .. ^dt_j .. dt_m  (= iota_E dt_1 ... dt_m).

wedge_c d^4 y_c = D(t) omega_m,  D(t) = (-1)^{j0} / t_{j0} * det[ d y_{c,a} / d t_j ]_{(c,a); j != j0}  (j0 0-based).
In ratio coordinates xi with a section t(xi):  I(lam; 0) = int g(xi) d xi,  g = D(t(xi)) |det[t; dt/dxi]|.
"""
from fractions import Fraction as Fr
import numpy as np


def solve(A, cols):
    n = len(A)
    k = len(cols)
    aug = [list(A[i]) + [cols[c][i] for c in range(k)] for i in range(n)]
    for i in range(n):
        p = next(r for r in range(i, n) if aug[r][i] != 0)
        aug[i], aug[p] = aug[p], aug[i]
        piv = aug[i][i]
        aug[i] = [x / piv for x in aug[i]]
        for r in range(n):
            if r != i and aug[r][i] != 0:
                f = aug[r][i]
                aug[r] = [x - f * z for x, z in zip(aug[r], aug[i])]
    return [[aug[i][n + c] for i in range(n)] for c in range(k)]


def fdet(M):
    M = [list(r) for r in M]
    n = len(M)
    det = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            det = -det
        det *= M[i][i]
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / M[i][i]
                M[r] = [x - f * z for x, z in zip(M[r], M[i])]
    return det


def kirchhoff(nv, edges, t, lam, H=4):
    """Currents y[e][a] and derivatives dy[e][a][j] = d y_{e,a} / d t_j (exact).
    lam: dict vertex -> H-vector (root excluded)."""
    g = nv - 1
    m = len(edges)
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(edges, t):
        c = 1 / Fr(te)
        for (i, j, s) in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
            if i < g and j < g:
                L[i][j] += s * c

    def bvec(e):
        a, b = edges[e]
        v = [Fr(0)] * g
        if a < g:
            v[a] += 1
        if b < g:
            v[b] -= 1
        return v

    rhs = [[Fr(lam.get(v, [0] * H)[c]) for v in range(g)] for c in range(H)] + [bvec(e) for e in range(m)]
    sol = solve(L, rhs)
    phi = [s + [Fr(0)] for s in sol[:H]]
    Lib = [s + [Fr(0)] for s in sol[H:]]
    y = [[(phi[c][a] - phi[c][b]) / t[e] for c in range(H)] for e, (a, b) in enumerate(edges)]
    dy = []
    for e, (a, b) in enumerate(edges):
        row = []
        for c in range(H):
            rr = []
            for j in range(m):
                val = (y[j][c] / t[j]) * (Lib[j][a] - Lib[j][b]) / t[e]
                if j == e:
                    val -= y[e][c] / t[e]
                rr.append(val)
            row.append(rr)
        dy.append(row)
    return y, dy


def check_kirchhoff(nv, edges, t, lam, y, H=4):
    g = nv - 1
    for v in range(g):
        for c in range(H):
            s = sum((y[e][c] if a == v else 0) - (y[e][c] if b == v else 0) for e, (a, b) in enumerate(edges))
            assert s == Fr(lam.get(v, [0] * H)[c])
    return True


def density_D(nv, edges, t, lam, cotree, H=4, j0=None):
    m = len(edges)
    j0 = m - 1 if j0 is None else j0
    y, dy = kirchhoff(nv, edges, t, lam, H)
    rows = [[dy[c][a][j] for j in range(m) if j != j0] for c in cotree for a in range(H)]
    return Fr((-1) ** j0) / t[j0] * fdet(rows), y


def cotree_of(nv, edges):
    """Kruskal tree in edge order; returns the cotree (increasing)."""
    parent = list(range(nv))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree = []
    for e, (a, b) in enumerate(edges):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            tree.append(e)
    return [e for e in range(len(edges)) if e not in tree]


# ------------------------------------------------------------------ floating point (vectorized)

def density_float_batch(nv, edges, T, lam, cotree, H=4):
    """D(t) for a batch T (P x m) of positive Schwinger times, float64.  lam: dict vertex -> H-vector."""
    P, m = T.shape
    g = nv - 1
    B = np.zeros((m, g))
    for e, (a, b) in enumerate(edges):
        if a < g:
            B[e, a] += 1
        if b < g:
            B[e, b] -= 1
    Lam = np.array([lam.get(v, [0] * H) for v in range(g)], dtype=float)       # g x H
    inv = 1.0 / T                                                              # P x m
    Lp = np.einsum('ev,pe,ew->pvw', B, inv, B)                                 # P x g x g
    phi = np.linalg.solve(Lp, np.broadcast_to(Lam, (P, g, H)))                 # P x g x H
    Bt = np.broadcast_to(B.T, (P, g, m))
    Lib = np.linalg.solve(Lp, Bt)                                              # P x g x m  (Lambda^{-1} b_j)
    y = np.einsum('ev,pvc->pec', B, phi) * inv[:, :, None]                     # P x m x H
    BL = np.einsum('ev,pvj->pej', B, Lib)                                      # (b_e . Lambda^{-1} b_j)
    # dy[p,e,c,j] = y[p,j,c]/t_j * BL[p,e,j]/t_e - delta_ej y[p,e,c]/t_e
    dy = (y.transpose(0, 2, 1)[:, None, :, :] * inv[:, None, None, :]) * BL[:, :, None, :] * inv[:, :, None, None]
    idx = np.arange(m)
    dy[:, idx, :, idx] -= (y * inv[:, :, None]).transpose(1, 0, 2)
    j0 = m - 1
    rows = dy[:, cotree, :, :j0].reshape(P, len(cotree) * H, j0)
    return (-1) ** j0 / T[:, j0] * np.linalg.det(rows)
