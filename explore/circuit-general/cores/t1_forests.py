"""Exact (H,T) = (1,1) kernels of the n = 4 cubic product forests and tadpole(5,0) of the
ht-d2 run (graphs and sources from explore/circuit-general/ht-d2/forests.py, read only):
root sources (2,-1) for two-tree forests, 2 for single trees, lambda_A = -2 (nonsp_numeric.py /
mc_chain.py conventions).  Prints K/pi^2."""
import importlib.util
import os
import sys
import time
from fractions import Fraction as F

import cores

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('htd2_forests', os.path.join(HERE, '..', 'ht-d2', 'forests.py'))
fm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fm)


def show(v):
    """value as a polynomial in pi^2 (labels sqrtpi^(2k))"""
    out = {}
    for k, x in v.items():
        d = dict(k)
        s = d.pop('sqrtpi', 0)
        assert not d and s % 2 == 0, v
        out[s // 2] = out.get(s // 2, 0) + x
    return ' + '.join('(%s) pi^%d' % (x, p) for p, x in sorted(out.items()) if x) or '0'


names = sys.argv[1:]
n = 4
src = [2, -1, 3, 1]
for Fo in fm.forests(n):
    G, roots = fm.forest_graph(Fo)
    if names and G.name not in names:
        continue
    lam = [[0] for _ in range(n)] + [[-2]]
    for i, v in enumerate(roots):
        lam[v] = [src[i]]
    t = time.time()
    try:
        val = cores.kernel(G.nv, G.edges, lam, (), verbose=False, T=1)
        print('%-28s K = %s   K/pi^2 ~ %s  (%.1fs)' % (G.name, show(val), _f if (_f := None) else
              '%.5f' % (sum(float(x) * 3.141592653589793 ** (dict(k).get('sqrtpi', 0) / 2) for k, x in val.items()) / 9.869604401089358),
              time.time() - t), flush=True)
    except Exception as e:
        print('%-28s FAILED %s: %s' % (G.name, type(e).__name__, str(e)[:200]), flush=True)
