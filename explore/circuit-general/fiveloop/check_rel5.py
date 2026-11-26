"""Exact check of the (Y), (X), (P) relations of rel5.py on cores with known periods:
K4, W4, W5 (wheel.py), prism (tables), K5-e (dual prism), K4+K4 (2-sum)."""
import random
import time
from fractions import Fraction as F

import per5
import gutil
import rel5

K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
W4 = [(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (1, 2), (2, 3), (0, 3)]
CASES = [('K4', K4, (3, 4, 5, 6, 7)), ('W4', W4, (3, 4, 5, 6, 7)), ('prism', list(per5.PRISM), (4, 5)),
         ('W5', list(per5.SKEL['W5']), (3, 4, 5, 6)), ('K5-e', list(per5.SKEL['K5-e']), (4, 5)),
         ('K4+K4', list(per5.SKEL['K4+K4']), (4, 5))]
rng = random.Random(5)


def check_items(items):
    tot = {}
    for e, b, c in items:
        v = per5.value(e, b, q)
        if v is None:
            return None
        tot = per5.vadd(tot, v, c)
    return tot == {}


t0 = time.time()
for name, E, qs in CASES:
    n = gutil.cycle_rank(E)
    tris = rel5.triangles(E)
    xr = rel5.x_rules(E)
    sdp = rel5.self_dual_perms(E)
    for q in qs:
        reps = [list(c) for c in gutil.compositions(n * q, len(E), 1, q - 1) if gutil.convergent(tuple(E), list(c), q)]
        rng.shuffle(reps)
        reps = reps[:150]
        cnt = {'Y': [0, 0], 'X': [0, 0], 'P': [0, 0]}
        for r in reps:
            for tri in tris:
                if r[tri[3]] + r[tri[4]] + r[tri[5]] == 2 * q:
                    ne, npw = rel5.tri_to_star(E, r, tri, q)
                    ok = check_items([(E, r, F(1)), (ne, npw, F(-1))])
                    if ok is not None:
                        cnt['Y'][0] += 1
                        cnt['Y'][1] += (not ok)
            for rule in xr:
                items = rel5.x_items(E, r, rule, q)
                if items is None:
                    continue
                ok = check_items(items)
                if ok is not None:
                    cnt['X'][0] += 1
                    cnt['X'][1] += (not ok)
            for perm in sdp:
                b = [q - r[perm[k]] for k in range(len(r))]
                if min(b) >= 1:
                    ok = check_items([(E, r, F(1)), (E, b, F(-1))])
                    if ok is not None:
                        cnt['P'][0] += 1
                        cnt['P'][1] += (not ok)
        print('%-6s q=%d: %d vectors; (Y) %d checked %d bad; (X) %d checked %d bad; (P) %d checked %d bad' % (
            name, q, len(reps), cnt['Y'][0], cnt['Y'][1], cnt['X'][0], cnt['X'][1], cnt['P'][0], cnt['P'][1]), flush=True)
print('%.1fs' % (time.time() - t0))
