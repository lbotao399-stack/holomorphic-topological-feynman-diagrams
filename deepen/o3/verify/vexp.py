"""Explicit evaluation at fixed N (explicit matrices, exact Fractions / ints)."""
import itertools, random
from fractions import Fraction as Fr
import numpy as np
from vlib import cycles_of, LP, seq_sign


def mat(n):
    return np.zeros((n, n), dtype=object)


class Algebra:
    def __init__(self, N, sl):
        self.N, self.sl = N, sl
        B = []
        for i in range(N):
            for j in range(N):
                if i != j:
                    M = mat(N) + 0
                    M[i, j] = 1
                    B.append(M)
        if sl:
            for i in range(N - 1):
                M = mat(N) + 0
                M[i, i] = 1
                M[i + 1, i + 1] = -1
                B.append(M)
        else:
            for i in range(N):
                M = mat(N) + 0
                M[i, i] = 1
                B.append(M)
        self.basis = B
        n = len(B)
        import sympy as sp
        K = sp.Matrix(n, n, lambda a, b: int(np.trace(B[a].dot(B[b]))))
        Ki = K.inv()
        self.kinv = [[Fr(int(Ki[a, b].p), int(Ki[a, b].q)) for b in range(n)] for a in range(n)]
        self.dim = n

    def random_element(self, r=2):
        N = self.N
        M = mat(N) + 0
        for i in range(N):
            for j in range(N):
                M[i, j] = random.randint(-r, r)
        if self.sl:
            t = sum(M[i, i] for i in range(N))
            M[N - 1, N - 1] -= t
        return M


def tau_eval(perm, Xs):
    val = 1
    for c in cycles_of(perm):
        M = Xs[c[0]]
        for s in c[1:]:
            M = M.dot(Xs[s])
        val *= np.trace(M)
    return val


def mono_eval(mono, k, ps):
    n = k - 1
    v = Fr(1)
    for s in range(n):
        v *= Fr(ps[s][0]) ** mono[s] * Fr(ps[s][1]) ** mono[n + s]
    return v


def evaluate(w, k, Xs, ps, N):
    tot = Fr(0)
    cache = {}
    for (perm, mono), c in w.items():
        if perm not in cache:
            cache[perm] = tau_eval(perm, Xs)
        cv = c.at(N) if isinstance(c, LP) else Fr(c)
        tot += cv * cache[perm] * mono_eval(mono, k, ps)
    return tot


def d_explicit(w, k, Xs, ps, N):
    tot = Fr(0)
    K = k + 1
    for i in range(K):
        for j in range(i + 1, K):
            br = Xs[i].dot(Xs[j]) - Xs[j].dot(Xs[i])
            pij = (ps[i][0] + ps[j][0], ps[i][1] + ps[j][1])
            rest = [l for l in range(K) if l not in (i, j)]
            args = [br] + [Xs[l] for l in rest]
            pa = [pij] + [ps[l] for l in rest]
            tot += (-1) ** (i + j) * evaluate(w, k, args, pa, N)
    return tot


def bracket_explicit(a, m, b, n, Xs, ps, alg):
    """shuffle formula with explicit dual basis."""
    K = m + n - 2
    tot = Fr(0)
    B = alg.basis
    for S in itertools.combinations(range(K), m - 1):
        Sc = [i for i in range(K) if i not in S]
        sg = seq_sign(list(S) + Sc)
        qe = (-sum(ps[i][0] for i in S), -sum(ps[i][1] for i in S))
        qe2 = (-qe[0], -qe[1])
        av = [evaluate(a, m, [Xs[i] for i in S] + [B[e]], [ps[i] for i in S] + [qe], alg.N) for e in range(alg.dim)]
        bv = [evaluate(b, n, [B[e]] + [Xs[i] for i in Sc], [qe2] + [ps[i] for i in Sc], alg.N) for e in range(alg.dim)]
        s = Fr(0)
        for e in range(alg.dim):
            if av[e] == 0:
                continue
            for e2 in range(alg.dim):
                if alg.kinv[e][e2] and bv[e2]:
                    s += av[e] * alg.kinv[e][e2] * bv[e2]
        tot += sg * s
    return tot
