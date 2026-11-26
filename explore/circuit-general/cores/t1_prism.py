"""T = 1 example: the four-loop prism at (H,T) = (1,1) (3d N=2 holomorphic-topological twist),
zero topological source.  Exact kernel from the half-integer-q prism tables, plus the
cross-level dimension relation between q = 7/2 and q = 9/2 as an independent check."""
import time
from fractions import Fraction as F

import cores
import gutil

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]

# cross-level check (R1) for the half-integer tables
edges = cores.FAMILIES['prism']
lo, hi = cores.table('prism', F(7, 2)), cores.table('prism', F(9, 2))
autos = gutil.edge_automorphisms(edges)
cot = gutil.cotrees(edges)
for a, va in lo.items():
    tot = {}
    for C in cot:
        b = list(a)
        for e in C:
            b[e] += 1
        for k, x in hi[gutil.orbit_rep(tuple(b), autos)].items():
            tot[k] = tot.get(k, 0) + x
    assert {k: x for k, x in tot.items() if x} == va
print('R1 between q=7/2 and q=9/2 holds on all %d orbits' % len(lo), flush=True)

for lam in ([1, 2, -1, 3, -2, -3], [1, 0, -1, 0, 0, 0], [1, 1, 0, -1, 0, -1]):
    t = time.time()
    v = cores.kernel(6, PRISM, [[x] for x in lam], (), verbose=True, T=1)
    print('prism (H,T)=(1,1), sources', lam, ':', {cores.enc(k): str(x) for k, x in v.items()},
          '%.1fs' % (time.time() - t), flush=True)
