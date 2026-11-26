"""Moments of graph 4.11 (non-tight plan): block-triangular XP kernel versus htk spred generic, exact."""
import random
import time
from itertools import cycle

import numpy as np

from kern import GraphKernel, modq, P
from rooted_trees import graphs
from htk.conventions import S2_kernel

g = GraphKernel(4, 11, list(graphs(4))[10])
rng = random.Random(3)
basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
probes = [dict(zip(g.outputs, cycle(basis))), {v: tuple(rng.randint(-2, 2) for _ in range(3)) for v in g.outputs}]
L = len(probes)
lam = {v: np.array([[probes[r][v][a] % P for r in range(L)] for a in range(3)], dtype=np.int64) for v in g.outputs}
m = len(g.spokes)
for (j1, a1), (j2, a2) in [((0, 0), (m - 1, 2)), ((1, 1), (2, 1))]:
    t0 = time.time()
    mx = g.moment_xp(lam, L)
    val = g.kernel_numeric(lam, L, extra=[mx[j1][a1], mx[j2][a2]])
    res = []
    for r, mom in enumerate(probes):
        ref = S2_kernel(g.red.integrate(mom, powers={(g.spokes[j1], a1): 1, (g.spokes[j2], a2): 1}, check=0), 4)
        res.append((str(ref), int(val[r]) == modq(ref)))
    print('4.11 moment (p_%d)_%d (p_%d)_%d: %s [%.0fs]' % (j1, a1 + 1, j2, a2 + 1, res, time.time() - t0), flush=True)
