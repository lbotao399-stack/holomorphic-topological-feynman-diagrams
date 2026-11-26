"""Exact Q_n kernels of the cone sets not yet evaluated: cone(D=4, n=3) at (H,T)=(4,0) (8d SYM Q_3)
and cone(D=3, n=4) at (H,T)=(2,1) (5d MSYM Q_4).  htk.kernel_exact (series-parallel reduction,
self-checked against the brute-force density), random integer sources at the bivalent vertices,
root = O.  usage: python3 e1_cone_kernels.py"""
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../../../code")))

import htk
from htk.graph import from_networkx
from atlas import set_cone

rng = random.Random(20261123)
for (D, n, H, T) in [(4, 3, 4, 0), (3, 4, 2, 1)]:
    items = set_cone(D, n)
    print("== cone D=%d n=%d, (H,T)=(%d,%d): %d (Gamma, O)" % (D, n, H, T, len(items)), flush=True)
    for idx, (G, O, S0) in enumerate(items):
        Gh, lab = from_networkx(G, root=O)
        lam = {lab[v]: [rng.randint(-3, 3) for _ in range(H)] for v in S0 if v != O}
        t0 = time.time()
        try:
            K = htk.kernel_exact(Gh, lam, T=T, check=1)
            msg = str(K)
        except Exception as e:
            msg = "FAILED: %r" % (e,)
        print("  #%d |V|=%d |E|=%d sources=%d  K = %s  (%.1fs)" % (idx, Gh.nv, len(Gh.edges), len(lam), msg, time.time() - t0), flush=True)
