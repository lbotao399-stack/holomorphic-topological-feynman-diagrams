"""Exact check of the conformal-completion relation (C) of rel5.py on cores with known periods."""
import random
import time

import per5
import gutil
import rel5

K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
W4 = [(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (1, 2), (2, 3), (0, 3)]
CASES = [('K4', K4, (3, 4, 5, 6)), ('W4', W4, (3, 4, 5, 6)), ('prism', list(per5.PRISM), (4, 5, 6)),
         ('K33', list(per5.K33), (4, 5, 6)), ('W5', list(per5.SKEL['W5']), (3, 4, 5)),
         ('K5-e', list(per5.SKEL['K5-e']), (4, 5)), ('K4+K4', list(per5.SKEL['K4+K4']), (4, 5))]
rng = random.Random(7)
t0 = time.time()
for name, E, qs in CASES:
    n = gutil.cycle_rank(E)
    for q in qs:
        reps = [list(c) for c in gutil.compositions(n * q, len(E), 1, q - 1) if gutil.convergent(tuple(E), list(c), q)]
        rng.shuffle(reps)
        nrel = bad = skipped = 0
        targets = {}
        for r in reps[:400]:
            for items in rel5.completion_items(E, r, q):
                vals = [per5.value(e, b, q) for e, b, c in items]
                if any(v is None for v in vals):
                    skipped += 1
                    continue
                tot = {}
                for (e, b, c), v in zip(items, vals):
                    tot = per5.vadd(tot, v, c)
                nrel += 1
                bad += (tot != {})
                ne = items[1][0]
                nl = gutil.cycle_rank(ne)
                targets[nl] = targets.get(nl, 0) + 1
        print('%-6s q=%d: %d relations checked, %d bad, %d skipped (divergent/unknown); target loop numbers %s' % (
            name, q, nrel, bad, skipped, dict(sorted(targets.items()))), flush=True)
print('%.1fs' % (time.time() - t0))
