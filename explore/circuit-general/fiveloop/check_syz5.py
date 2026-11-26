"""Exact check of the (L) syzygy relations (syz5.py) on cores with known periods."""
import random
import time

import per5
import gutil
import syz5

K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
W4 = [(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (1, 2), (2, 3), (0, 3)]
CASES = [('K4', K4, (1, 2), (4, 5)), ('W4', W4, (1, 2), (4, 5)), ('prism', list(per5.PRISM), (1, 2), (4, 5)),
         ('K33', list(per5.K33), (1, 2), (4, 5)), ('W5', list(per5.SKEL['W5']), (1, 2), (3, 4)),
         ('K5-e', list(per5.SKEL['K5-e']), (1, 2), (4,))]
rng = random.Random(3)
t0 = time.time()
for name, E, ks, qs in CASES:
    n = gutil.cycle_rank(E)
    N = len(E)
    for k in ks:
        sy = syz5.syzygies(E, k)
        for q in qs:
            bases = list(gutil.compositions(n * q - k + 1, N, 1, q))
            rng.shuffle(bases)
            nrel = bad = skip = 0
            for beta in bases[:60]:
                for s in sy:
                    items = syz5.relation(E, s, list(beta), q)
                    if items is None:
                        skip += 1
                        continue
                    tot = {}
                    ok = True
                    for e, g, c in items:
                        if min(g) < 1:
                            ok = False
                            break
                        v = per5.value(e, g, q)
                        if v is None:
                            ok = False
                            break
                        tot = per5.vadd(tot, v, c)
                    if not ok:
                        skip += 1
                        continue
                    nrel += 1
                    bad += (tot != {})
            print('%-6s k=%d (%d syzygies) q=%d: %d relations checked, %d bad, %d skipped' % (
                name, k, len(sy), q, nrel, bad, skip), flush=True)
print('%.1fs' % (time.time() - t0))
