"""Kernel validation: XP-modular kernels and moments versus htk spred exact Fractions.
Usage: python3 test_kern.py N [graph numbers]"""
import sys
import time
import random
from fractions import Fraction as Fr
from itertools import cycle

import numpy as np

import kern
from kern import GraphKernel, modq, P, xp_const
from rooted_trees import graphs
from htk.conventions import S2_kernel

n = int(sys.argv[1])
only = [int(x) for x in sys.argv[2:]]
rng = random.Random(7)
for number, edges in enumerate(graphs(n), 1):
    if only and number not in only:
        continue
    t0 = time.time()
    g = GraphKernel(n, number, edges)
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    probes = [dict(zip(g.outputs, cycle(basis)))]
    for _ in range(2):
        probes.append({v: tuple(rng.randint(-3, 3) for _ in range(3)) for v in g.outputs})
    L = len(probes)
    lam = {v: np.array([[probes[r][v][a] % P for r in range(L)] for a in range(3)], dtype=np.int64)
           for v in g.outputs}
    val = g.kernel_numeric(lam, L)
    ok = []
    exact = []
    for r, mom in enumerate(probes):
        ref = S2_kernel(g.red.kernel(mom, check=0), n) if g.red.tight_plan or r == 0 else None
        exact.append(ref)
        if ref is not None:
            ok.append(int(val[r]) == modq(ref))
    # one moment: (p_0)_0 (p_last)_2 at the random probes
    m = len(g.spokes)
    mx = g.moment_xp(lam, L)
    valm = g.kernel_numeric(lam, L, extra=[mx[0][0], mx[m - 1][2]])
    okm = []
    for r, mom in enumerate(probes):
        if r == 0 or not g.red.tight_plan:
            continue
        e0, e1 = g.spokes[0], g.spokes[m - 1]
        # p_j = y_{e_j}
        ref = S2_kernel(g.red.integrate(mom, powers={(e0, 0): 1, (e1, 2): 1}, check=0), n)
        okm.append(int(valm[r]) == modq(ref))
    oks = None
    if g.red.tight_plan and n <= 2:
        lexp, coef = g.kernel_symbolic(((0,),) + ((),) * (m - 2) + ((2,),))
        args = []
        for nm in g.lam_names:
            v, a = map(int, nm[1:].split('_'))
            args.append(lam[v][a])
        vs = kern.xp_value(kern.eval_symbolic(lexp, coef, args, L))
        oks = [int(vs[r]) == int(valm[r]) for r in range(L)]
    print('graph %d.%d plan %s kappa %s: I(probe) = %s; XP==spred %s; moment XP==spred %s; symbolic==XP %s [%.1fs]'
          % (n, number, 'tight' if g.red.tight_plan else 'generic', g.kappa, exact[0], ok, okm, oks,
             time.time() - t0), flush=True)
