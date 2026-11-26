"""D=2, n=5: kernels of the 70 census graphs (atlas census.json, key '2,5') at generic integer sources
(all vertices, lambda_v in [-3,3]^2, seed 1000 + index, root = last vertex).
SP: htk.kernel_exact.  Others: exact series-skeleton numerator (fk5.core_numerator, saved in
cnum/<index>.json) contracted with per5 periods at q = 7.  usage: python3 run_d2.py [indices]"""
import json
import os
import random
import sys
import time
from fractions import Fraction as F

import fk5
import per5
import census_cores
import htk
from htk.graph import Graph

HERE = os.path.dirname(os.path.abspath(__file__))
cen = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))
graphs = cen['data']['2,5']['graphs']


def sources(g):
    rng = random.Random(1000 + g['index'])
    return [[rng.randint(-3, 3) for _ in range(2)] for _ in range(g['nv'] - 1)]


def run(i):
    g = graphs[i]
    nv, edges = g['nv'], [tuple(e) for e in g['edges']]
    lam = sources(g)
    core, _ = census_cores.core_name(edges)
    t = time.time()
    out = dict(index=i, core=core, sp=g['sp'], nv=nv, edges=edges, sources=lam)
    if g['sp']:
        K = htk.kernel_exact(Graph(nv, edges), lam, check=1)
        out['K'] = {'1': str(K)} if K else {}
        out['method'] = 'htk.kernel_exact'
    else:
        path = os.path.join(HERE, 'cnum', '%d.json' % i)
        if os.path.exists(path):
            d = json.load(open(path))
            ce, q = [tuple(e) for e in d['core_edges']], d['q']
            num = {tuple(json.loads(k)): F(v) for k, v in d['num'].items()}
        else:
            ce, q, num = fk5.core_numerator(nv, edges, lam, verbose=False)
            json.dump({'index': i, 'core_edges': ce, 'q': q, 'sources': lam,
                       'num': {json.dumps(list(a)): str(c) for a, c in num.items()}}, open(path, 'w'))
        out['core_vectors'] = len(num)
        out['tnum'] = round(time.time() - t, 1)
        try:
            K = fk5.kernel_from_core(ce, q, num)
            out['K'] = {k: str(v) for k, v in K.items()}
            out['method'] = 'fk5 + per5'
        except (ValueError, per5.UnknownCore) as e:
            out['K'] = None
            out['pending'] = str(e)[:120]
    out['seconds'] = round(time.time() - t, 1)
    return out


if __name__ == '__main__':
    idx = [int(x) for x in sys.argv[1:]] or list(range(len(graphs)))
    for i in idx:
        r = run(i)
        print(json.dumps(r), flush=True)
