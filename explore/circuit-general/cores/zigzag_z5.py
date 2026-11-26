"""Five-loop zigzag Z5 (C7(1,2) minus a vertex): regulator-free relation system at level q,
with boundary values from W_n (wheel.py) and the four-loop prism/K33 tables; then the
dimension relation  P(Z5; 1^10; q=2) = sum_C P(Z5; 1+1_C; 3)  against 441 zeta(7)/8."""
import sys
import time
from fractions import Fraction as F

import cores
import gutil
import msys

Z5 = cores.zigzag(5)
q = int(sys.argv[1]) if len(sys.argv) > 1 else 3
t = time.time()
S = msys.MultiSystem([Z5], q)
S.lower = cores.lower_core_value
S.build()
solved, piv = S.solve()
print('Z5 q=%d: %d unknowns, %d determined (%.1fs)' % (q, len(S.reps), len(solved), time.time() - t), flush=True)
if q == 3 and len(solved) == len(S.reps):
    autos = gutil.edge_automorphisms(Z5)
    tot = {}
    for C in gutil.cotrees(Z5):
        b = [1] * len(Z5)
        for e in C:
            b[e] += 1
        j = S.index[(0, gutil.orbit_rep(tuple(b), autos))]
        for k, x in solved[j].items():
            tot[k] = tot.get(k, 0) + x
    tot = {k: x for k, x in tot.items() if x}
    print('P(Z5) at q=2 from level 3:', tot, 'expected {z7: 441/8}', tot == {'z7': F(441, 8)}, flush=True)
