"""Block-triangular density: (i) tight graphs, against the spred local factors; (ii) graph 4.11 against the
spred generic kernel (exact Fractions) at the rotating-basis probe and one random probe."""
import sys
import time
import random
from itertools import cycle

import numpy as np

import kern
from kern import GraphKernel, modq, P, integrate_factors
from blocktri import BlockTriangular
from rooted_trees import graphs
from htk.conventions import S2_kernel

rng = random.Random(11)
for n, number in [(3, 1), (3, 4), (4, 1), (4, 10), (4, 11)]:
    edges = list(graphs(n))[number - 1]
    t0 = time.time()
    g = GraphKernel(n, number, edges)
    bt = BlockTriangular(g)
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    probes = [dict(zip(g.outputs, cycle(basis)))]
    probes.append({v: tuple(rng.randint(-3, 3) for _ in range(3)) for v in g.outputs})
    L = len(probes)
    lam = {v: np.array([[probes[r][v][a] % P for r in range(L)] for a in range(3)], dtype=np.int64)
           for v in g.outputs}
    arr = np.zeros((len(g.lam_names), L), dtype=np.int64)
    for v in g.outputs:
        for a in range(3):
            arr[g.lam_index['l%d_%d' % (v, a)]] = lam[v][a]
    val_bt = integrate_factors(bt.evaluate(arr, L), g.blocks, L) * (g.sign_n % P) % P
    t1 = time.time()
    if g.red.tight_plan:
        ref = g.kernel_numeric(lam, L)
        res = [int(a) == int(b) for a, b in zip(val_bt, ref)]
        print('%d.%d tight: Laplace terms %d, sign %+d, block det == local factors: %s [%.1fs]'
              % (n, number, len(bt.terms), bt.sign, res, t1 - t0), flush=True)
    else:
        refs = [S2_kernel(g.red.kernel(mom, check=0), n) for mom in probes]
        res = [int(val_bt[r]) == modq(refs[r]) for r in range(L)]
        print('%d.%d generic: Laplace terms %d, sign %+d, I = %s; block det == spred generic: %s [%.1fs + %.1fs]'
              % (n, number, len(bt.terms), bt.sign, [str(x) for x in refs], res, t1 - t0, time.time() - t1),
              flush=True)
