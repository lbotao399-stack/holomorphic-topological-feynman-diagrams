"""Faithful fixed-N coordinates: trace monomials -> functions on g^k (g = gl_N or sl_N).
relation data per (N, sl, k): evaluation matrix E[perm, col] on explicit tuples; a set of columns of full rank.
Exactness of the function-space dimension r_k:
  lower bound: exact rank of E on the chosen tuples;
  upper bound: either r_k = #perms, or exact identities (all basis tuples for small dim, or the antisymmetrizer
  identity sum_{s in S_k} sgn(s) tau_s = 0 on gl_N for k > N, restricted to sl)."""
import itertools, random
from fractions import Fraction as Fr
import numpy as np
from vlib import all_perms, cycles_of, has_fixed, perm_sign, LP
import vla


def basis_mats(N, sl):
    B = []
    for i in range(N):
        for j in range(N):
            if i != j:
                M = np.zeros((N, N), dtype=np.int64); M[i, j] = 1; B.append(M)
    if sl:
        for i in range(N - 1):
            M = np.zeros((N, N), dtype=np.int64); M[i, i] = 1; M[i + 1, i + 1] = -1; B.append(M)
    else:
        for i in range(N):
            M = np.zeros((N, N), dtype=np.int64); M[i, i] = 1; B.append(M)
    return B


def tau_val(perm, Xs):
    v = 1
    for c in cycles_of(perm):
        M = Xs[c[0]]
        for s in c[1:]:
            M = M @ Xs[s]
        v *= int(np.trace(M))
    return v


class Faithful:
    def __init__(self, N, sl, k, mode='random', npts=None, seed=0):
        self.N, self.sl, self.k = N, sl, k
        self.perms = all_perms(k, sl)
        self.pidx = {p: i for i, p in enumerate(self.perms)}
        rng = random.Random(seed + 1000 * N + 17 * k)
        B = basis_mats(N, sl)
        if mode == 'all':
            tuples = [tuple(B[a] for a in idx) for idx in itertools.product(range(len(B)), repeat=k)]
        else:
            npts = npts or 3 * len(self.perms) + 10
            tuples = []
            for _ in range(npts):
                Xs = []
                for _i in range(k):
                    M = np.array([[rng.randint(-3, 3) for _ in range(N)] for _ in range(N)], dtype=np.int64)
                    if sl:
                        M[N - 1, N - 1] -= int(np.trace(M))
                    Xs.append(M)
                tuples.append(tuple(Xs))
        E = [[tau_val(p, X) for X in tuples] for p in self.perms]
        rows = [{j: Fr(x) for j, x in enumerate(r) if x} for r in E]
        self.rank = vla.rankQ(rows)
        # choose columns: greedy via rref of the transpose
        cols_rows = [{i: Fr(E[i][j]) for i in range(len(self.perms)) if E[i][j]} for j in range(len(tuples))]
        from sympy.polys.matrices.sdm import SDM
        from sympy import QQ
        data = {}
        for j, r in enumerate(cols_rows):
            if r:
                data[j] = {i: QQ(int(x)) for i, x in r.items()}
        M = SDM(data, (len(tuples), len(self.perms)), QQ).transpose()
        _, piv = M.rref()
        self.cols = list(piv)
        assert len(self.cols) == self.rank
        self.E = E
        self.ntuples = len(tuples)

    def upper_bound_ok(self):
        """exact upper bound certificate for rank: rank == #perms, or antisymmetrizer relation (k = N+1 case)."""
        if self.rank == len(self.perms):
            return True, "rank = #taus"
        return False, "needs relations"

    def antisym_relation(self):
        """coefficient vector of sum_{s in S_k} sgn(s) tau_s restricted to sl (fixed-point terms dropped)."""
        v = {}
        for p in itertools.permutations(range(self.k)):
            if self.sl and has_fixed(p):
                continue
            v[self.pidx[p]] = Fr(perm_sign(p))
        return v

    def proj(self, w, N=None):
        N = self.N if N is None else N
        out = {}
        for (perm, mono), c in w.items():
            if perm not in self.pidx:
                continue
            cv = c.at(N) if isinstance(c, LP) else Fr(c)
            if not cv:
                continue
            i = self.pidx[perm]
            for jj, j in enumerate(self.cols):
                x = self.E[i][j]
                if x:
                    key = (mono, jj)
                    out[key] = out.get(key, 0) + cv * x
        return {k: v for k, v in out.items() if v}
