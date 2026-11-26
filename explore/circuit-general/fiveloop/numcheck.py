"""Float Monte Carlo cross-check of tabulated periods:  I_q(G; a) = prod Gamma(a_e) E[(sum t) / U(t)^q],
t_e ~ Gamma(a_e, 1) independent (homogeneity: the weight (sum t) e^{-sum t} integrates the scale to 1).
usage: python3 numcheck.py NAME Q NSAMPLES [NVEC]"""
import math
import random
import sys

import numpy as np
import networkx as nx

import per5
import gutil

ZETA = {3: 1.2020569031595942854, 5: 1.0369277551433699263, 7: 1.0083492773819228268}


def fval(v):
    tot = 0.0
    for k, x in v.items():
        f = float(x)
        if k != '1':
            for part in k.split('*'):
                f *= ZETA[int(part[1:])]
        tot += f
    return tot


def cycle_matrix(E):
    G = nx.Graph()
    G.add_edges_from(E)
    T = list(nx.minimum_spanning_tree(G).edges())
    Tset = {frozenset(e) for e in T}
    idx = {frozenset(e): k for k, e in enumerate(E)}
    Tg = nx.Graph(T)
    rows = []
    for k, (u, v) in enumerate(E):
        if frozenset((u, v)) in Tset:
            continue
        path = nx.shortest_path(Tg, v, u)
        r = np.zeros(len(E))
        r[k] = 1
        for a, b in zip(path, path[1:]):
            j = idx[frozenset((a, b))]
            r[j] += 1 if E[j] == (a, b) else -1
        rows.append(r)
    return np.array(rows)


def mc(E, a, q, n, seed=0, batch=200000):
    C = cycle_matrix(E)
    rng = np.random.default_rng(seed)
    vals = []
    for s in range(0, n, batch):
        t = np.stack([rng.gamma(x, 1.0, size=batch) for x in a], axis=1)
        M = np.einsum('re,be,se->brs', C, t, C)
        U = np.linalg.det(M)
        vals.append(t.sum(axis=1) / U ** q)
    v = np.concatenate(vals)
    pref = math.prod(math.gamma(x) for x in a)
    return pref * v.mean(), pref * v.std() / math.sqrt(len(v))


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
        m, e = mc(E, list(r), q, n)
        print('%s q=%d %s: exact %.6f  MC %.6f +- %.6f  (%.1f sigma)' % (name, q, r, ex, m, e, abs(m - ex) / e), flush=True)
