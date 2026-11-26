"""Lie algebras with integer structure constants and the trace form.

gl_N: basis E_ij (a = i N + j).  sl_N: E_ij (i != j), H_k = E_kk - E_{k+1,k+1}.
[X_a, X_b] = sum_c fu[a, b, c] X_c          (fu: structure constants, upper index last)
kappa_ab   = tr(X_a X_b)                     (trace form of the defining representation)
fl[a,b,c]  = kappa([X_a, X_b], X_c)          (totally antisymmetric, integer)
kinv_scaled = s * kappa^{-1}, integer, with the integer s returned.
"""
import numpy as np
from fractions import Fraction as Fr


class LieAlg:
    def __init__(self, name, mats, coords):
        self.name = name
        self.mats = mats            # list of N x N integer numpy arrays
        self.coords = coords        # function: matrix -> integer coordinate vector
        self.dim = len(mats)
        d = self.dim
        fu = np.zeros((d, d, d), dtype=np.int64)
        for a in range(d):
            for b in range(d):
                c = self.mats[a] @ self.mats[b] - self.mats[b] @ self.mats[a]
                fu[a, b] = self.coords(c)
        self.fu = fu
        kap = np.array([[int(np.trace(self.mats[a] @ self.mats[b])) for b in range(d)] for a in range(d)],
                       dtype=np.int64)
        self.kappa = kap
        self.fl = np.einsum('abe,ec->abc', fu, kap)
        # exact inverse of kappa, scaled to integers
        K = [[Fr(int(x)) for x in row] for row in kap]
        inv = _inv(K)
        den = 1
        for row in inv:
            for x in row:
                den = den * x.denominator // np.gcd(den, x.denominator)
        self.kinv_scale = den
        self.kinv = np.array([[int(x * den) for x in row] for row in inv], dtype=np.int64)
        # checks
        assert np.array_equal(self.fl, -np.transpose(self.fl, (1, 0, 2)))
        assert np.array_equal(self.fl, np.transpose(self.fl, (1, 2, 0)))
        assert np.array_equal(self.kappa @ self.kinv, den * np.eye(d, dtype=np.int64))
        # Jacobi
        J = (np.einsum('abe,ecf->abcf', fu, fu) + np.einsum('bce,eaf->abcf', fu, fu)
             + np.einsum('cae,ebf->abcf', fu, fu))
        assert not J.any()

    def center_dim(self):
        # dim of {x : [x, X_b] = 0 for all b}
        M = self.fu.transpose(0, 2, 1).reshape(self.dim, -1).astype(float)
        return self.dim - np.linalg.matrix_rank(M)


def _inv(A):
    n = len(A)
    M = [list(row) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
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


def E(N, i, j):
    m = np.zeros((N, N), dtype=np.int64)
    m[i, j] = 1
    return m


def gl(N):
    mats = [E(N, i, j) for i in range(N) for j in range(N)]

    def coords(m):
        return np.array([m[i, j] for i in range(N) for j in range(N)], dtype=np.int64)
    return LieAlg("gl%d" % N, mats, coords)


def sl(N):
    off = [(i, j) for i in range(N) for j in range(N) if i != j]
    mats = [E(N, i, j) for (i, j) in off] + [E(N, k, k) - E(N, k + 1, k + 1) for k in range(N - 1)]

    def coords(m):
        assert int(np.trace(m)) == 0
        v = [m[i, j] for (i, j) in off]
        d = [m[k, k] for k in range(N)]
        c = np.cumsum(d)[:N - 1]
        return np.array(v + list(c), dtype=np.int64)
    return LieAlg("sl%d" % N, mats, coords)
