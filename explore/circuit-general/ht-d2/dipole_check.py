"""Exact test (htk.brute, Fractions): at (H,T) = (1,1), sources lambda at e(0) and -lambda at e(1)
for an edge e (all other sources zero) give D_Gamma(t) = 0 at random rational t, for every edge of
every tight D = 2 graph with n = 2, 4 (catalog), and for comparison a non-edge pair (p, q).
Output: dipole_check.log"""
import itertools, random
from fractions import Fraction as Fr
from forests import Graph  # noqa: F401  (sets sys.path)
from htk import catalog, brute

rng = random.Random(9)
for n in (2, 4):
    gs = catalog.tight_graphs(2, n)
    tot_edge = zero_edge = tot_non = zero_non = 0
    for G in gs:
        adj = {frozenset(e) for e in G.edges}
        pairs = list(itertools.combinations(range(G.nv), 2))
        for (p, q) in pairs:
            lam = [[Fr(0)] for _ in range(G.nv)]
            lam[p] = [Fr(1)]
            lam[q] = [Fr(-1)]
            vals = []
            for _ in range(2):
                t = [Fr(rng.randint(1, 30), rng.randint(1, 9)) for _ in range(G.ne)]
                vals.append(brute.density(G, t, lam, T=1).coeff)
            z = all(v == 0 for v in vals)
            if frozenset((p, q)) in adj:
                tot_edge += 1
                zero_edge += z
            else:
                tot_non += 1
                zero_non += z
    print("n=%d: %d graphs; dipole across an edge: D = 0 at 2 points for %d of %d (graph, edge); "
          "dipole on a non-adjacent pair: D = 0 for %d of %d" % (n, len(gs), zero_edge, tot_edge, zero_non, tot_non), flush=True)
