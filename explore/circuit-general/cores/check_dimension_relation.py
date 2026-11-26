"""Independent check of the period tables: dimension relation (R1)
   I(a; q-1) = sum_{C cotree} I(a + 1_C; q)
for every convergent orbit a at level q-1 whose raised vectors are all in the level-q table
(R1 alone is not among the relations used to solve a single level)."""
import random
import sys
import time

import cores
import gutil


def check(name, q, sample=None):
    edges = cores.FAMILIES[name]
    lo, hi = cores.table(name, q - 1), cores.table(name, q)
    autos = gutil.edge_automorphisms(edges)
    cot = gutil.cotrees(edges)
    keys = list(lo)
    if sample and len(keys) > sample:
        keys = random.Random(5).sample(keys, sample)
    n_ok = 0
    for a in keys:
        tot = {}
        for C in cot:
            b = list(a)
            for e in C:
                b[e] += 1
            v = hi[gutil.orbit_rep(tuple(b), autos)]
            for k, x in v.items():
                tot[k] = tot.get(k, 0) + x
        tot = {k: x for k, x in tot.items() if x}
        assert tot == lo[a], (name, q, a, tot, lo[a])
        n_ok += 1
    return n_ok


if __name__ == '__main__':
    for name, q in [('prism', 7), ('K33', 7)]:
        t = time.time()
        n = check(name, q)
        print('%s: R1 between q=%d and q=%d holds on all %d orbits of level %d  (%.1fs)' % (name, q - 1, q, n, q - 1, time.time() - t), flush=True)
