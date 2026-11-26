"""Explicit evaluation of trace-basis cochains on matrices (independent of the formal rules)."""
import random
import itertools
from fractions import Fraction as Fr
import ce


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def ident(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def tau_value(tau, mats):
    v = Fr(1)
    n = len(mats[0])
    for c in tau:
        M = ident(n)
        for i in c:
            M = matmul(M, mats[i])
        v *= trace(M)
    return v


def mono_value(m, moms):
    v = Fr(1)
    for s in range(len(m) // 2):
        for c in range(2):
            if m[2 * s + c]:
                v *= Fr(moms[s][c]) ** m[2 * s + c]
    return v


def eval_cochain(c, mats, moms, N):
    """value of canonical cochain at x_i = mats[i] e^{moms[i] z}; moms must satisfy sum = 0."""
    k = len(mats)
    assert all(sum(Fr(mm[cc]) for mm in moms) == 0 for cc in range(2))
    tot = Fr(0)
    for (tau, m), npl in c.items():
        tot += ce.np_eval(npl, N) * tau_value(tau, mats) * mono_value(m, moms)
    return tot


def rand_mat(N, R=5, traceless=False, rng=random):
    A = [[Fr(rng.randint(-R, R)) for _ in range(N)] for _ in range(N)]
    if traceless:
        t = trace(A)
        A[N - 1][N - 1] -= t
    return A


def rand_moms(k, R=5, rng=random):
    ms = [[Fr(rng.randint(-R, R)), Fr(rng.randint(-R, R))] for _ in range(k - 1)]
    ms.append([-sum(m[0] for m in ms), -sum(m[1] for m in ms)])
    return ms


def basis_and_dual(N, sl=False):
    """basis X_e of gl_N or sl_N and dual basis X^e w.r.t. tr(XY)."""
    def E(i, j):
        M = [[Fr(0)] * N for _ in range(N)]
        M[i][j] = Fr(1)
        return M
    if not sl:
        B = [E(i, j) for i in range(N) for j in range(N)]
        D = [E(j, i) for i in range(N) for j in range(N)]
        return B, D
    B, D = [], []
    for i in range(N):
        for j in range(N):
            if i != j:
                B.append(E(i, j))
                D.append(E(j, i))
    Hs = []
    for k in range(N - 1):
        M = E(k, k)
        M[k + 1][k + 1] = Fr(-1)
        Hs.append(M)
    G = [[trace(matmul(a, b)) for b in Hs] for a in Hs]
    # inverse of G
    n = len(G)
    Aug = [list(G[i]) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Aug[r][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        pv = Aug[c][c]
        Aug[c] = [x / pv for x in Aug[c]]
        for r in range(n):
            if r != c and Aug[r][c] != 0:
                f = Aug[r][c]
                Aug[r] = [x - f * y for x, y in zip(Aug[r], Aug[c])]
    Gi = [row[n:] for row in Aug]
    for a in range(n):
        M = [[sum(Gi[a][b] * Hs[b][i][j] for b in range(n)) for j in range(N)] for i in range(N)]
        B.append(Hs[a])
        D.append(M)
    return B, D


def bracket_explicit(a, m, b, n, mats, moms, N, sl=False):
    """{a,b}(x_1..x_K) by the shuffle formula with explicit dual bases."""
    K = m + n - 2
    Bs, Ds = basis_and_dual(N, sl)
    tot = Fr(0)
    for comb in itertools.combinations(range(K), m - 1):
        restl = [i for i in range(K) if i not in comb]
        sg = ce.perm_sign(list(comb) + restl)
        xs = [mats[i] for i in comb]
        xm = [moms[i] for i in comb]
        ys = [mats[i] for i in restl]
        ym = [moms[i] for i in restl]
        qx = [sum(Fr(p[c]) for p in xm) for c in range(2)]
        for Xe, Xd in zip(Bs, Ds):
            va = eval_cochain(a, xs + [Xe], xm + [[-qx[0], -qx[1]]], N)
            if va == 0:
                continue
            vb = eval_cochain(b, [Xd] + ys, [qx] + ym, N)
            tot += sg * va * vb
    return tot
