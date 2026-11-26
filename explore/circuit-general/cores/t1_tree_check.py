"""((A,(A,B)),(A,B)) at (1,1): exact kernels at (l0, lA) = (1,0), (1,1), (2,-2) for comparison
with the ht-d2 Monte Carlo source polynomial (run_mc_poly / mc_summary.log)."""
import importlib.util
import os
import sys
import time

import cores

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('htd2_forests', os.path.join(HERE, '..', 'ht-d2', 'forests.py'))
fm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fm)
name = sys.argv[1] if len(sys.argv) > 1 else '((A,(A,B)),(A,B))'
for Fo in fm.forests(4):
    G, roots = fm.forest_graph(Fo)
    if G.name != name:
        continue
    print(G.name, G.edges, 'roots', roots, flush=True)
    for l0, lA in ((1, 0), (1, 1), (2, -2)):
        lam = [[0] for _ in range(4)] + [[lA]]
        lam[roots[0]] = [l0]
        t = time.time()
        v = cores.kernel(G.nv, G.edges, lam, (), T=1)
        f = sum(float(x) * 3.141592653589793 ** (dict(k).get('sqrtpi', 0) / 2) for k, x in v.items()) / 9.869604401089358
        print('  (l0,lA)=(%d,%d): K = %s  K/pi^2 = %.6f  (%.0fs)' % (l0, lA, {cores.enc(k): str(x) for k, x in v.items()}, f, time.time() - t), flush=True)
