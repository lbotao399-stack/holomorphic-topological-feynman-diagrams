"""Float Monte Carlo cross-check (not exact) of the T=1 prism kernel: Dirichlet(a) importance
sampling of the simplex with the htk float density (oracle.density_batch)."""
import math
import sys
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../../code'))
from htk import oracle
from htk.graph import Graph

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
G = Graph(6, PRISM)
lam = [[1.0], [2.0], [-1.0], [3.0], [-2.0]]
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
a = 0.5
N = 9
vals = []
for b in range(40):
    x = rng.dirichlet([a] * N, size=200000)
    D = oracle.density_batch(G, x, lam, T=1)
    # Dirichlet(a) density on the simplex (w.r.t. dt_0..dt_{N-2})
    logq = math.lgamma(a * N) - N * math.lgamma(a) + (a - 1) * np.sum(np.log(x), axis=1)
    w = D * np.exp(-logq)
    w = w[np.isfinite(w)]
    vals.append(w.mean())
vals = np.array(vals)
exact = -808 / 63 * math.pi ** 2 + 2041 / 1680 * math.pi ** 4
print('MC %.4f +- %.4f   exact %.4f' % (vals.mean(), vals.std(ddof=1) / math.sqrt(len(vals)), exact))
