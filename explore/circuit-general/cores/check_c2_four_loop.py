"""C^2 four-loop prism and K33 kernels at the nonSP fixtures (FORCER values in
explore/nonSP/verify/prism_forcer/README.md and k33_forcer/README.md)."""
import sys
import time
from fractions import Fraction as F

import cores

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
K33 = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5)]
CASES = [
    ('prism', PRISM, [(-1, 0), (1, 1), (-1, -1), (0, -1), (1, 0), (0, 1)], (),
     {'1': F(-577, 40320), 'z3': F(-2, 35), 'z5': F(5, 63)}),
    ('prism (1,0),(0,1),(-1,0),(0,-1),0,0', PRISM, [(1, 0), (0, 1), (-1, 0), (0, -1), (0, 0), (0, 0)], (),
     {'1': F(169, 151200), 'z3': F(-3, 2800)}),
    ('K33', K33, [(-3, -1), (-2, 2), (0, -1), (1, -1), (3, 1), (1, 0)], (),
     {'1': F(-789, 640), 'z3': F(59, 80)}),
]
which = sys.argv[1:] or ['0', '1', '2']
for w in which:
    name, edges, lam, ins, ref = CASES[int(w)]
    t = time.time()
    v = cores.kernel(6, edges, lam, ins, verbose=True)
    print(name, v, 'MATCH' if v == ref else 'MISMATCH ref %s' % ref, '%.1fs' % (time.time() - t), flush=True)
