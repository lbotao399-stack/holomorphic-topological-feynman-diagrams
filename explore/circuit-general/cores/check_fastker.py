"""fastker (vectorized modular numerator) against the exact fixtures: C^2 four-loop prism and
K33 (FORCER), C^2 three-loop K4 subdivision, C^3 three-loop K4 subdivisions (k4_periods)."""
from fractions import Fraction as F
import time

import fastker
from check_c3_three_loop import K4SUB, MOM

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
K33 = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5)]
cases = [('C2 prism', 6, PRISM, [(-1, 0), (1, 1), (-1, -1), (0, -1), (1, 0), (0, 1)], {'1': F(-577, 40320), 'z3': F(-2, 35), 'z5': F(5, 63)}),
         ('C2 K33', 6, K33, [(-3, -1), (-2, 2), (0, -1), (1, -1), (3, 1), (1, 0)], {'1': F(-789, 640), 'z3': F(59, 80)}),
         ('C2 K4sub', 5, [(0, 1), (0, 2), (0, 3), (1, 4), (4, 2), (1, 3), (2, 3)], [(1, 2), (2, -1), (-1, 3), (-3, -4), (1, 0)], {'1': F(-71, 12)})]
for name, edges, (r, z) in K4SUB:
    cases.append(('C3 ' + name, 8, edges, MOM, {k: x for k, x in (('1', r), ('z3', z)) if x}))
for name, nv, edges, lam, ref in cases:
    t = time.time()
    v = fastker.kernel(nv, edges, lam, verbose=False)
    print('%-20s %s %s %.1fs' % (name, {k: str(x) for k, x in v.items()}, 'MATCH' if v == ref else 'MISMATCH', time.time() - t), flush=True)
