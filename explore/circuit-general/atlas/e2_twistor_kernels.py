"""Exact twistor BF Q_n kernels (D=3, T=0): n=3 all graphs, n=4 series-parallel graphs.
SP: htk.kernel_exact.  n=3 non-SP: explore/nonSP/verify/k4_periods/kernel.py evaluate (Q + Q zeta(3)),
cross-checked against htk.kernel_exact on the SP graphs of the same set.
Random integer sources in C^3 at the bivalent vertices, root = O (source = minus the sum).
usage: python3 e2_twistor_kernels.py"""
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../../../code")))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../../nonSP/verify/k4_periods")))

import networkx as nx
import htk
from htk.graph import from_networkx
from atlas import set_twistor
try:
    from kernel import evaluate
except ImportError as _e:      # python-flint missing
    evaluate = None
    print('k4_periods evaluator unavailable: %r' % (_e,))

H = 3
rng = random.Random(20261123)
for n in (3, 4):
    sets = set_twistor(n)
    items = [it for key in ('tree', 'unicyclic') for it in sets[key]]
    print("== twistor n=%d: %d (Gamma, O)" % (n, len(items)), flush=True)
    for idx, item in enumerate(items):
        G, O, S0 = item[0], item[1], item[2]
        Gh, lab = from_networkx(G, root=O)
        sp = htk.catalog.is_sp(nx.Graph(G))
        lam = {lab[v]: [rng.randint(-3, 3) for _ in range(H)] for v in S0 if v != O}
        full = [lam.get(v, [0] * H) for v in range(Gh.nv - 1)]
        full.append([-sum(r[i] for r in full) for i in range(H)])
        t0 = time.time()
        out = []
        if sp:
            try:
                out.append("spred K = %s" % htk.kernel_exact(Gh, lam, T=0, check=1))
            except Exception as e:
                out.append("spred FAILED %r" % (e,))
        if n == 3 and evaluate is not None:
            try:
                a, b = evaluate(list(Gh.edges), full)
                out.append("k4_periods K = %s + %s zeta(3)" % (a, b))
            except Exception as e:
                out.append("k4_periods FAILED %r" % (e,))
        if not sp and (n == 4 or evaluate is None):
            out.append("non-SP: not evaluated")
        print("  #%d SP=%s |V|=%d |E|=%d  %s  (%.1fs)" % (idx, sp, Gh.nv, len(Gh.edges), " | ".join(out), time.time() - t0), flush=True)
