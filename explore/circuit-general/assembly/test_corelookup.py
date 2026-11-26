"""CoreLookup (vectorized) against periods.period (msys.classify) at random power vectors."""
import random
import time
import numpy as np
import periods
from corelookup import CoreLookup
from modlin import fmod

P = 2147483629
cases = [('prism', [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)], 7),
         ('K33', [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)], 6),
         ('prism relabelled', [(2, 0), (5, 1), (3, 4), (0, 1), (3, 5), (1, 2), (2, 4), (0, 3), (4, 5)], 7),
         ('K4', [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], 7),
         ('W4', [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (1, 4)], 6),
         ('cubic+double', [(0, 1), (0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 5)], 6),
         ('W4 double', [(0, 1), (0, 1), (0, 2), (2, 3), (0, 3), (1, 3), (1, 2), (2, 3)], 6),
         ('prism q8', [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)], 8),
         ('K33 q8', [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)], 8),
         ('W4 q9', [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (1, 4)], 9),
         ('prism+double q8', [(0, 1), (0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (2, 5)], 8)]
rng = random.Random(3)
for name, edges, q in cases:
    n = len(edges) - len({x for e in edges for x in e}) + 1
    rows = []
    while len(rows) < 300:
        a = [rng.randint(1, 2 * q) for _ in edges]
        s = sum(a)
        if s == n * q:
            rows.append(a)
    t = time.time()
    cl = CoreLookup(edges, q)
    prep = cl.prepare(np.array(rows))
    vals, ok = cl.evaluate(np.array(rows), prep, P)
    t1 = time.time() - t
    bad = nconv = 0
    for k, a in enumerate(rows):
        try:
            v = periods.period(edges, a, q)
            conv = True
        except ValueError:
            conv = False
        if conv != bool(ok[k]):
            bad += 1
            continue
        if conv:
            nconv += 1
            for lab in set(v) | set(vals):
                if fmod(v.get(lab, 0), P) != int(vals.get(lab, np.zeros(len(rows), dtype=np.int64))[k]):
                    bad += 1
    print('%-18s kind=%-6s q=%d rows %d convergent %d mismatches %d  (%.2fs)' % (name, cl.kind, q, len(rows), nconv, bad, t1), flush=True)
