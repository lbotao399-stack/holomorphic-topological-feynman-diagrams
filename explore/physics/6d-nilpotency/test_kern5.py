"""Five-loop kernels (49 rooted graphs) by the XP machinery (4-bit exponent packing, 15 coordinates) versus the
exact witnesses of explore/six-d-Q3/five_loop_periods.jsonl (htk spred).  Usage: python3 test_kern5.py"""
import json
import os
import time
from fractions import Fraction as Fr

import numpy as np

import kern
kern.BITS = 4
kern.MASK = 15
from kern import GraphKernel, modq, P
from rooted_trees import graphs

recs = [json.loads(l) for l in open(os.path.join(kern.PROOT, 'explore', 'six-d-Q3', 'five_loop_periods.jsonl'))]
ok = 0
for rec, edges in zip(recs, graphs(5)):
    t0 = time.time()
    try:
        g = GraphKernel(5, rec['number'], edges)
        mom = {int(k): v for k, v in rec['momenta'].items()}
        lam = {v: np.array([[mom.get(v, [0, 0, 0])[a] % P] for a in range(3)], dtype=np.int64) for v in g.outputs}
        val = int(g.kernel_numeric(lam, 1)[0])
        same = val == modq(Fr(rec['value']))
        kap = str(g.kappa) == rec['wick_weight']
        ok += same and kap
        msg = 'I == witness %s, kappa %s == wick_weight %s' % (same, g.kappa, kap)
    except Exception as e:                      # report and continue
        msg = 'FAILED: %s: %s' % (type(e).__name__, e)
    print('5.%d %s: %s [%.1fs]' % (rec['number'], rec['plan'], msg, time.time() - t0), flush=True)
print('agreeing graphs: %d of %d' % (ok, len(recs)))
