"""Exact T = 1, H = 1 spider kernels at n = 6 (and n = 4 for comparison) at integer Q_n sources
(htk.spred, self-checked against brute force): values q pi^3.  Tests the banana value
K_(1^n) = (-1)^(n/2) pi^(n/2) prod lambda_i / (n+1).  Output: spiders_n6_numeric.log"""
import sys, time
from fractions import Fraction as Fr
from forests import spider_graph
from spiders_exact import arm_lists
from htk import sp_reduce

srcs = [(2, -1, 3, 1, -2, 1, 3), (1, 1, 2, -3, 1, 2, -1)]
for n in (4, 6):
    for arms in arm_lists(n):
        arms = tuple(sorted(arms))
        G, outer = spider_graph(arms)
        red = sp_reduce(G, prefer="large")
        for s in srcs:
            lam = [[0] for _ in range(n)] + [[s[-1]]]
            for i, v in enumerate(outer):
                lam[v] = [s[i]]
            t0 = time.time()
            K = red.kernel(lam, T=1, check=1)
            extra = ""
            if all(l == 1 for l in arms):
                p = Fr(1)
                for i in range(len(arms)):
                    p *= s[i]
                pred = (-1) ** (n // 2) * p / (n + 1)
                extra = "  banana prediction %s pi^%d: %s" % (pred, n // 2, K.coeff(n // 2) == pred)
            print("arms %s sources %s lamA %d: K = %s  (%.1fs)%s" % (arms, s[:len(arms)], s[-1], K, time.time() - t0, extra), flush=True)
