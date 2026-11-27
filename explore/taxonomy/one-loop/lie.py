"""Exact matrix Lie algebras and representations (Fractions) for the supertrace checks."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, combinations


def mat(n, m=None):
    m = n if m is None else m
    return [[F(0)] * m for _ in range(n)]


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    out = mat(n, m)
    for i in range(n):
        Ai = A[i]
        for l in range(k):
            a = Ai[l]
            if a:
                Bl = B[l]
                row = out[i]
                for j in range(m):
                    if Bl[j]:
                        row[j] += a * Bl[j]
    return out


def madd(A, B, c=1):
    return [[a + c * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, c):
    return [[c * a for a in r] for r in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def mpow(A, k):
    n = len(A)
    out = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    for _ in range(k):
        out = mmul(out, A)
    return out


def tr_pow(A, k):
    return trace(mpow(A, k))


def block_diag(*Ms):
    n = sum(len(M) for M in Ms)
    out = mat(n)
    o = 0
    for M in Ms:
        for i in range(len(M)):
            for j in range(len(M)):
                out[o + i][o + j] = M[i][j]
        o += len(M)
    return out


def comm(A, B):
    return madd(mmul(A, B), mmul(B, A), -1)


# ---------------------------------------------------------------- algebras as lists of basis matrices (defining rep)

def E(n, i, j):
    M = mat(n)
    M[i][j] = F(1)
    return M


def basis_sl(n):
    B = [E(n, i, j) for i in range(n) for j in range(n) if i != j]
    for i in range(n - 1):
        B.append(madd(E(n, i, i), E(n, i + 1, i + 1), -1))
    return B


def basis_so(n):
    return [madd(E(n, i, j), E(n, j, i), -1) for i in range(n) for j in range(i + 1, n)]


def basis_sp(n2):
    """sp(2m) preserving J = [[0,1],[-1,0]] blocks: X^T J + J X = 0"""
    m = n2 // 2
    J = mat(n2)
    for i in range(m):
        J[i][m + i] = F(1)
        J[m + i][i] = F(-1)
    # X = J S with S symmetric
    B = []
    for i in range(n2):
        for j in range(i, n2):
            S = mat(n2)
            S[i][j] = F(1)
            S[j][i] = F(1)
            B.append(mmul(J, S))
    return B, J


def coords(B, X):
    """solve X = sum c_a B_a (least squares not needed: exact elimination on flattened entries)"""
    n = len(X)
    rows = [[B[a][i][j] for a in range(len(B))] + [X[i][j]] for i in range(n) for j in range(n)]
    # gaussian elimination
    m = len(B)
    R = [r[:] for r in rows]
    piv = []
    r0 = 0
    for c in range(m):
        p = next((r for r in range(r0, len(R)) if R[r][c] != 0), None)
        if p is None:
            continue
        R[r0], R[p] = R[p], R[r0]
        inv = 1 / R[r0][c]
        R[r0] = [x * inv for x in R[r0]]
        for r in range(len(R)):
            if r != r0 and R[r][c] != 0:
                f = R[r][c]
                R[r] = [x - f * y for x, y in zip(R[r], R[r0])]
        piv.append(c)
        r0 += 1
    assert all(all(x == 0 for x in R[r][:m]) and R[r][m] == 0 for r in range(r0, len(R))), "not in span"
    sol = [F(0)] * m
    for k, c in enumerate(piv):
        sol[c] = R[k][m]
    return sol


def adjoint(B, X):
    """matrix of ad_X in basis B"""
    m = len(B)
    out = mat(m)
    for b in range(m):
        c = coords(B, comm(X, B[b]))
        for a in range(m):
            out[a][b] = c[a]
    return out


def sym2(X):
    """X acting on Sym^2 C^n, basis e_i e_j (i <= j)"""
    n = len(X)
    idx = list(combinations_with_replacement(range(n), 2))
    pos = {p: k for k, p in enumerate(idx)}
    out = mat(len(idx))
    for b, (i, j) in enumerate(idx):
        # X(e_i e_j) = (X e_i) e_j + e_i (X e_j)
        for k in range(n):
            for (u, v, c) in ((k, j, X[k][i]), (i, k, X[k][j])):
                if c:
                    out[pos[tuple(sorted((u, v)))]][b] += c
    return out


def alt2(X):
    n = len(X)
    idx = list(combinations(range(n), 2))
    pos = {p: k for k, p in enumerate(idx)}
    out = mat(len(idx))
    for b, (i, j) in enumerate(idx):
        for k in range(n):
            for (u, v, c) in ((k, j, X[k][i]), (i, k, X[k][j])):
                if c and u != v:
                    s = 1 if u < v else -1
                    out[pos[(min(u, v), max(u, v))]][b] += s * c
    return out


def dual(M):
    return mscale(transpose(M), -1)


def killing_form(B):
    m = len(B)
    ads = [adjoint(B, b) for b in B]
    return [[trace(mmul(ads[a], ads[b])) for b in range(m)] for a in range(m)]
