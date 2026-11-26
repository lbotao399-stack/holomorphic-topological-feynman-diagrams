"""Float cross-check of periods by tropical sampling (Hepp sectors; Borinsky arXiv:2008.12310):
  I_q(G; a) = int prod x^a U^-q prod dx/x  (projective),  omega(g) = sum_{e in g} a_e - q L(g) > 0.
Sector sigma (x_{s_1} < ... < x_{s_N} = 1, x_{s_k} = prod_{j>=k} y_j): tropical integrand prod_j y_j^{omega(g_j)-1},
g_j = {s_1..s_j};  I^tr = T(G),  T(g) = sum_{e in g} T(g - e) / omega(g - e),  T(single edge) = 1.
Samples x ~ f^tr / I^tr;  I = I^tr E[(U^tr / U)^q],  U^tr = max spanning-tree-complement monomial
(Kruskal on increasing x), ratio in [#trees^-q, 1] (finite variance).
usage: python3 tropmc.py NAME Q NSAMPLES [NVEC]"""
import math
import random
import sys
from functools import lru_cache

import numpy as np

import per5
import gutil
from numcheck import fval


class Trop:
    def __init__(self, E, a, q):
        self.E, self.a, self.q = list(E), list(a), q
        N = len(E)
        self.N = N
        full = (1 << N) - 1
        self.L = np.zeros(1 << N, dtype=np.int64)
        self.om = np.zeros(1 << N)
        for m in range(1, 1 << N):
            sub = [E[i] for i in range(N) if m >> i & 1]
            self.L[m] = gutil.cycle_rank(sub)
            self.om[m] = sum(a[i] for i in range(N) if m >> i & 1) - q * self.L[m]
        for m in range(1, full):
            if self.om[m] <= 0:
                raise ValueError('divergent subgraph')
        T = np.zeros(1 << N)
        for m in range(1, 1 << N):
            bits = [i for i in range(N) if m >> i & 1]
            if len(bits) == 1:
                T[m] = 1.0
                continue
            T[m] = sum(T[m ^ (1 << i)] / self.om[m ^ (1 << i)] for i in bits)
        self.T = T
        self.Itr = T[full]
        # cumulative choice table: from mask m, remove edge i with prob T[m-i]/om[m-i] / T[m]
        self.cum = np.zeros((1 << N, N))
        for m in range(1, 1 << N):
            p = np.zeros(N)
            for i in range(N):
                if m >> i & 1 and m != (1 << i):
                    p[i] = T[m ^ (1 << i)] / self.om[m ^ (1 << i)]
            s = p.sum()
            if s > 0:
                self.cum[m] = np.cumsum(p / s)
        cot = gutil.cotrees(tuple(E))
        self.Cmat = np.zeros((N, len(cot)))
        for k, C in enumerate(cot):
            for e in C:
                self.Cmat[e, k] = 1.0

    def sample(self, B, rng):
        N, q = self.N, self.q
        mask = np.full(B, (1 << N) - 1, dtype=np.int64)
        order = np.zeros((B, N), dtype=np.int64)          # order[:, k] = s_{k+1}
        logy = np.zeros((B, N))                           # logy[:, j-1] = log y_j, j = 1..N-1
        for k in range(N - 1, 0, -1):                     # choose s_{k+1} (largest remaining), remaining g_k
            u = rng.random(B)
            c = self.cum[mask]
            e = np.argmax(c > u[:, None], axis=1)
            order[:, k] = e
            mask = mask ^ (1 << e)
            om = self.om[mask]
            logy[:, k - 1] = np.log(rng.random(B)) / om   # y ~ om y^(om-1)
        order[:, 0] = np.log2(mask).astype(np.int64)
        # log x_{s_k} = sum_{j >= k} log y_j  (k = 1..N, j up to N-1)
        cs = np.cumsum(logy[:, ::-1], axis=1)[:, ::-1]    # cs[:, k-1] = sum_{j>=k} logy_j
        logx = np.zeros((B, N))
        rows = np.arange(B)
        for k in range(N):
            logx[rows, order[:, k]] = cs[:, k] if k < N - 1 else 0.0
        # tropical U: cotree edges of the Kruskal tree on increasing x = s_1, s_2, ...
        logutr = np.zeros(B)
        m = np.zeros(B, dtype=np.int64)
        prevL = np.zeros(B, dtype=np.int64)
        for k in range(N):
            e = order[:, k]
            m = m | (1 << e)
            Lk = self.L[m]
            logutr += (Lk - prevL) * logx[rows, e]
            prevL = Lk
        mono = logx @ self.Cmat - logutr[:, None]
        ratio = np.exp(mono).sum(axis=1)                  # U / U^tr >= 1
        return ratio ** (-q)

    def estimate(self, n, seed=0, batch=100000):
        rng = np.random.default_rng(seed)
        vals = np.concatenate([self.sample(min(batch, n - s), rng) for s in range(0, n, batch)])
        return self.Itr * vals.mean(), self.Itr * vals.std() / math.sqrt(len(vals))


if __name__ == '__main__':
    name, q, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nvec = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    E = list(per5.SKEL[name]) if name in per5.SKEL else list(per5.FOUR[name])
    tab = per5.table(name, q)
    rng = random.Random(4)
    reps = sorted(tab)
    rng.shuffle(reps)
    for r in reps[:nvec]:
        ex = fval(tab[r])
        m, e = Trop(E, list(r), q).estimate(n)
        print('%s q=%d %s: exact %.6f  tropical MC %.6f +- %.6f  (%.1f sigma)' % (name, q, r, ex, m, e, abs(m - ex) / e),
              flush=True)
