"""fk5 (per5 periods) against the C^2 four-loop FORCER fixtures and htk.kernel_exact on a
five-loop series-parallel D=2 fan.  No untracked census file is required."""
import random
import time
from fractions import Fraction as F

import fk5
import htk
from htk.graph import Graph, fan

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
K33 = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5)]
cases = [('C2 prism', 6, PRISM, [(-1, 0), (1, 1), (-1, -1), (0, -1), (1, 0), (0, 1)], {'1': F(-577, 40320), 'z3': F(-2, 35), 'z5': F(5, 63)}),
         ('C2 K33', 6, K33, [(-3, -1), (-2, 2), (0, -1), (1, -1), (3, 1), (1, 0)], {'1': F(-789, 640), 'z3': F(59, 80)})]
for name, nv, edges, lam, ref in cases:
    t = time.time()
    v = fk5.kernel(nv, edges, lam, verbose=False)
    print('%-20s %s %s %.1fs' % (name, {k: str(x) for k, x in v.items()}, 'MATCH' if v == ref else 'MISMATCH', time.time() - t), flush=True)
    if v != ref:
        raise AssertionError('%s period regression failed' % name)
G = fan(5)
rng = random.Random(11)
nv = G.nv
edges = G.edges
lam = [[rng.randint(-3, 3) for _ in range(2)] for _ in range(nv - 1)]
t = time.time()
ref = htk.kernel_exact(Graph(nv, edges), lam, check=1)
t1 = time.time()
v = fk5.kernel(nv, edges, lam, verbose=True)
print('n=5 fan: kernel_exact %s (%.1fs), fk5 %s (%.1fs) %s' % (ref, t1 - t, v, time.time() - t1,
      'MATCH' if v == ({'1': F(ref)} if ref else {}) else 'MISMATCH'), flush=True)
if v != ({'1': F(ref)} if ref else {}):
    raise AssertionError('five-loop fan period regression failed')
