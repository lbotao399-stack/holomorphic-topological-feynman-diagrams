"""D = 4 (H = 4, T = 0) three-loop kernels with K4 core (census data['4,3']): fastker with the
wheel engine (K4 = W3) at q = 7; integer sources at all vertices."""
import json
import os
import sys
import time
from collections import Counter

import fastker
import msys
import cores

HERE = os.path.dirname(os.path.abspath(__file__))
gs = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))['data']['4,3']['graphs']
K4key = msys.struct_key([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
found = []
for i, g in enumerate(gs):
    es = [tuple(e) for e in g['edges']]
    while True:
        deg = Counter(x for e in es for x in e)
        v2 = [v for v, c in deg.items() if c == 2]
        if not v2:
            break
        v = v2[0]
        inc = [e for e in es if v in e]
        a = inc[0][0] if inc[0][1] == v else inc[0][1]
        b = inc[1][0] if inc[1][1] == v else inc[1][1]
        es = [e for e in es if v not in e] + [(a, b)]
    if len(es) == 6 and msys.struct_key(es) == K4key:
        found.append(i)
print('D=4 n=3 graphs with K4 series core:', len(found), flush=True)
import random
rng = random.Random(2)
for i in found[:int(sys.argv[1]) if len(sys.argv) > 1 else 2]:
    g = gs[i]
    nv = g['nv']
    src = [tuple(rng.randint(-2, 2) for _ in range(4)) for _ in range(nv - 1)]
    src.append(tuple(-sum(v[k] for v in src) for k in range(4)))
    t = time.time()
    val = fastker.kernel(nv, [tuple(e) for e in g['edges']], src, verbose=True)
    print('census[4,3][%d] nv=%d K = %s (%.1fs)' % (i, nv, {k: str(v) for k, v in val.items()}, time.time() - t), flush=True)
