"""D=2 n=5 graphs with a source pattern S0 for which the graph is B-class (atlas reduce_search: S, P, Y0
reduce to one edge): the kernel is rational [atlas REPORT §1.4: B => K in Q[lambda]]; check that the
zeta parts of fk5 + per5 cancel.  S0 = bivalent vertices plus EXTRA.
usage: python3 check_gen_rational.py INDEX[:EXTRA,EXTRA,...] ..."""
import json
import os
import random
import sys
import time

import fk5
import census_cores

HERE = os.path.dirname(os.path.abspath(__file__))
cen = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))
graphs = cen['data']['2,5']['graphs']
for arg in sys.argv[1:]:
    i = int(arg.split(':')[0])
    extra = [int(x) for x in arg.split(':')[1].split(',')] if ':' in arg else []
    g = graphs[i]
    nv, edges = g['nv'], [tuple(e) for e in g['edges']]
    deg = [0] * nv
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    rng = random.Random(2000 + i)
    biv = [v for v in range(nv) if deg[v] == 2] + extra
    full = {v: [0, 0] for v in range(nv)}
    for v in biv[:-1]:
        full[v] = [rng.randint(-3, 3) for _ in range(2)]
    full[biv[-1]] = [-sum(full[v][j] for v in biv[:-1]) for j in range(2)]
    lam = [full[v] for v in range(nv - 1)]
    t = time.time()
    ce, q, num = fk5.core_numerator(nv, edges, lam, verbose=False)
    K = fk5.kernel_from_core(ce, q, num)
    core, _ = census_cores.core_name(edges)
    print(json.dumps({'index': i, 'core': core, 'B': g['B'], 'sources': lam, 'S0': biv,
                      'K': {k: str(v) for k, v in K.items()}, 'rational': set(K) <= {'1'},
                      'seconds': round(time.time() - t, 1)}), flush=True)
