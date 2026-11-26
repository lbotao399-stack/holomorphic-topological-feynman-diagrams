"""D = 3 (H = 3, T = 0) four-loop kernels at zero displacement: graphs 191 (prism core) and 219
(K33 core) of explore/circuit-general/atlas/census.json (data['3,4']), integer sources at all
ten vertices.  fastker: modular Newton interpolation of U^7 D on 4384432 grid points, series
push-forward to the core, exact q = 7 core tables, CRT + rational reconstruction."""
import json
import os
import sys
import time

import fastker

HERE = os.path.dirname(os.path.abspath(__file__))
gs = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))['data']['3,4']['graphs']
SRC = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1), (2, 0, -1), (0, -2, 1), (-1, 2, 0), (-3, -2, -2)]
which = [int(a) for a in sys.argv[1:]] or [191, 219]
for i in which:
    g = gs[i]
    edges = [tuple(e) for e in g['edges']]
    t = time.time()
    val = fastker.kernel(g['nv'], edges, SRC)
    print('census[3,4][%d] %s  K = %s  (%.0fs)' % (i, edges, {k: str(v) for k, v in val.items()}, time.time() - t), flush=True)
