"""Probe: T = 1, H = 1 kernels of the cubic product forests, n = 2, 4: SP ones exactly (spred),
non-SP ones flagged.  Output: probe_forests.log"""
import sys, time
import sympy as sp
from forests import forests, forest_graph
from htk import catalog, sp_reduce, brute
from fractions import Fraction as Fr

for n in (2, 4):
    print("n =", n, flush=True)
    for F in forests(n):
        G, roots = forest_graph(F)
        sp_ok = catalog.is_series_parallel(G)
        tight = G.tight(2) and catalog.subgraph_bound_ok(G, 2)
        syms = sp.symbols("l0:%d" % n)
        lamA = sp.Symbol("lA")
        lam = [[syms[v]] if v in roots else [0] for v in range(n)] + [[lamA]]
        line = "  %-40s roots %s tight %s SP %s" % (G.name, roots, tight, sp_ok)
        if sp_ok:
            t0 = time.time()
            red = sp_reduce(G, prefer="large")
            K = red.kernel(lam, T=1)
            line += "  K = %s  (%.1fs)" % (sp.factor(K), time.time() - t0)
        print(line, flush=True)
