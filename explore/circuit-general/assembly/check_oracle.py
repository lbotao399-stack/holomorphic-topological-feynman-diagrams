"""Assembled moments (exact) vs htk.oracle (independent numerical integration over the Schwinger simplex,
no reductions), C^3 three-loop K4 subdivision K4[1,1,1,2,2,3]a, sources MOM (cores/check_c3_three_loop.py)."""
import sys
import time
sys.path.insert(0, '../cores')
import assembly
from check_c3_three_loop import K4SUB, MOM
from htk.graph import Graph
from htk import oracle

Z3, Z5 = 1.2020569031595942, 1.0369277551433699


def fl(v):
    return float(v.get('1', 0)) + float(v.get('z3', 0)) * Z3 + float(v.get('z5', 0)) * Z5


name, edges, _ = K4SUB[int(sys.argv[1]) if len(sys.argv) > 1 else 1]
res = assembly.kernel_moments(8, edges, MOM, order=2, verbose=False)
G = Graph(8, edges)
lam = [list(map(float, v)) for v in MOM[:-1]]
for key, powers in [('bare', None), (('y', 3, 1), {(3, 1): 1}), (('y', 6, 2), {(6, 2): 1}),
                    (('yy', 3, 1, 5, 2), {(3, 1): 1, (5, 2): 1}), (('yy', 0, 0, 0, 0), {(0, 0): 2})]:
    t = time.time()
    val, err = oracle.kernel(G, lam, T=0, powers=powers, method='qmc', npoints=2 ** 20, nrep=16)
    ex = fl(res.get(key, {}))
    print('%s %-22s exact %.6e  oracle %.6e +- %.1e  (diff/err %.2f) %.0fs' % (name, key, ex, val, err, (val - ex) / err if err else 0, time.time() - t), flush=True)
