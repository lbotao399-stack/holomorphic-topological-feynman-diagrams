"""Exact pointwise T = 1, H = 1 densities (htk.brute) at random rational times for graphs whose
Monte Carlo value was 0 +- 0: tadpole(3,2), tadpole(4,1), tadpole(5,0), (A,(B,(B,(A,B)))) and the
single trees at generic sources.  Output: zero_density_check.log"""
import os, sys, random
from fractions import Fraction as Fr
from forests import tadpole_graph, forests, forest_graph
from htk import brute
rng = random.Random(5)
cases = []
for (c, L) in ((3, 2), (4, 1), (5, 0)):
    G, q = tadpole_graph(c, L)
    lam = [[Fr(0)] for _ in range(G.nv - 1)]
    lam[q] = [Fr(1)]
    cases.append((G, lam, "lambda_q = 1"))
for F in forests(4):
    G, roots = forest_graph(F)
    if len(F) == 1:
        for (l0, lA) in ((2, -2), (2, 1)):
            lam = [[Fr(0)] for _ in range(G.nv - 1)]
            lam[roots[0]] = [Fr(l0)]
            lam[G.nv - 2] = [Fr(lA)]
            cases.append((G, lam, "l0=%d lA=%d" % (l0, lA)))
for G, lam, tag in cases:
    vals = []
    for _ in range(4):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 9)) for _ in range(G.ne)]
        vals.append(brute.density(G, t, lam, T=1).coeff)
    print("%-24s %-12s D at 4 random rational t: %s" % (G.name, tag, "all zero" if all(v == 0 for v in vals) else "nonzero"), flush=True)
