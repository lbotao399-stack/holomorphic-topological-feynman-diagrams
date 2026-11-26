"""Exact T = 1, H = 1 spider kernels (LLC Definition b-00380) by htk.spred, self-checked against the
brute-force layer.  Q_n sources: lambda at the outer arm vertices v_{i,1} and at the hub A; inner arm
vertices source-free; B = root.  Also: general sources for small cases.  Output: spiders_exact.log"""
import itertools
import sys
import time

import sympy as sp

from forests import spider_graph
from htk import sp_reduce, catalog


def arm_lists(n, amin=3):
    out = []
    def rec(rem, mx, acc):
        if rem == 0:
            if len(acc) >= amin:
                out.append(tuple(acc))
            return
        for l in range(min(rem, mx), 0, -1):
            rec(rem - l, l, acc + [l])
    rec(n, n, [])
    return out


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(3, nmax + 1):
        for arms in arm_lists(n):
            arms_sorted = tuple(sorted(arms))
            G, outer = spider_graph(arms_sorted)
            a = len(arms_sorted)
            ls = sp.symbols("l1:%d" % (a + 1))
            lA = sp.Symbol("lA")
            lam = [[0] for _ in range(n)] + [[lA]]
            for i, v in enumerate(outer):
                lam[v] = [ls[i]]
            t0 = time.time()
            assert catalog.is_series_parallel(G) and G.tight(2)
            red = sp_reduce(G, prefer="large")
            K = red.kernel(lam, T=1, check=1)
            print("arms %s n=%d: K = %s   (%.1fs)" % (arms_sorted, n, sp.factor(K), time.time() - t0), flush=True)
