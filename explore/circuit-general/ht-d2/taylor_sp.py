"""Exact displacement Taylor coefficients at (H,T) = (1,1) (htk.spred.taylor, self-checked) for the
n = 4 spiders with Q_n sources and the fan-tree forests:
    K(lam; z) = sum_m K_m,  K_m = int g (-sum_e y_e z_e)^m / m!.
Symbolic sources and symbolic displacements on the spokes (v, B) and on (A, B).  Output: taylor_sp.log"""
import time

import sympy as sp

from forests import spider_graph
from htk import sp_reduce

for arms in ((1, 1, 2), (1, 1, 1, 1)):
    G, outer = spider_graph(arms)
    n = sum(arms)
    a = len(arms)
    ls = sp.symbols("l1:%d" % (a + 1))
    lA = sp.Symbol("lA")
    lam = [[0] for _ in range(n)] + [[lA]]
    for i, v in enumerate(outer):
        lam[v] = [ls[i]]
    B = G.root
    z = {}
    for e, (p, q) in enumerate(G.edges):
        if q == B:
            z[e] = [sp.Symbol("z%d" % p if p < n else "zA")]
    red = sp_reduce(G, prefer="large")
    t0 = time.time()
    Ks = red.taylor(lam, z, 2, T=1, check=1)
    print("spider %s (edges %s, root B=%d):" % (arms, G.edges, B), flush=True)
    for m, Km in enumerate(Ks):
        print("  K_%d = %s" % (m, sp.factor(Km)), flush=True)
    print("  (%.1fs)" % (time.time() - t0), flush=True)
