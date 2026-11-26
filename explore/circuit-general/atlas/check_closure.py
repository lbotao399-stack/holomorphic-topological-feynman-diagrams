"""Exhaustive search (no safe closure) versus the closure search: same k(Gamma, S) for every graph
(patterns gen and all) in the ranges given on the command line."""
import sys
import time
import reduce_search as rs
from graphs_enum import admissible


def ks(G, S):
    net = rs.Network(G, S)
    r = rs.search(net, ("S", "Y0"))
    if r["reducible"]:
        return 0
    return rs.search(net, ("S", "Y0", "Y1"))["min_Y1"]


rng = [tuple(map(int, a.split(","))) for a in sys.argv[1:]]
for D, n in rng:
    t0 = time.time()
    diff = 0
    tot = 0
    for G in admissible(D, n):
        biv = [v for v in G if G.degree(v) == 2]
        for S in (biv, list(G.nodes())):
            rs.CLOSURE = True
            a = ks(G, S)
            rs.CLOSURE = False
            b = ks(G, S)
            tot += 1
            if a != b:
                diff += 1
    print("D=%d n=%d: %d (graph, pattern) pairs, closure vs exhaustive differ: %d | %.1fs" % (D, n, tot, diff, time.time() - t0), flush=True)
