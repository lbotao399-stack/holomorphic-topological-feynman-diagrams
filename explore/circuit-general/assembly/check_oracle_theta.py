"""Assembled moments (exact) vs htk.oracle lattice integration (no reductions), D = 3 two-loop thetas
theta(1,3,3), theta(2,2,3) (N - 1 = 6), integer sources at the interior vertices."""
import time
import assembly
from htk.graph import theta
from htk import oracle

Z3 = 1.2020569031595942


def fl(v):
    return float(v.get('1', 0)) + float(v.get('z3', 0)) * Z3


for (p, q, r), lamint in [((1, 3, 3), [[1, 0, 2], [0, 1, -1], [2, -1, 0], [-1, 2, 1]]),
                          ((2, 2, 3), [[1, 1, 0], [0, 2, -1], [1, 0, 1], [-2, 1, 1]])]:
    G, interior = theta(p, q, r)
    lam = [[0, 0, 0] for _ in range(G.nv)]
    for v, l in zip([v for ints in interior for v in ints], lamint):
        lam[v] = l
    lam[G.nv - 1] = [-sum(x[i] for x in lam[:-1]) for i in range(3)]
    res = assembly.kernel_moments(G.nv, G.edges, lam, order=2, verbose=False)
    lamf = [list(map(float, v)) for v in lam[:-1]]
    for key, powers in [('bare', None), (('y', 0, 0), {(0, 0): 1}), (('y', 3, 2), {(3, 2): 1}),
                        (('yy', 1, 1, 5, 0), {(1, 1): 1, (5, 0): 1}), (('yy', 2, 2, 2, 2), {(2, 2): 2})]:
        t = time.time()
        val, err = oracle.kernel(G, lamf, T=0, powers=powers)
        ex = fl(res.get(key, {}))
        print('theta(%d,%d,%d) %-20s exact %+.8e  oracle %+.8e +- %.1e  (diff/err %+.2f) %.0fs' % (
            p, q, r, key, ex, val, err, (val - ex) / err if err else 0, time.time() - t), flush=True)
