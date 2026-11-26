"""Checks of per5.lower: W5, K4+K4 at q=2; 2-sum rule vs series-parallel rules; dimension relation
I_{q-1}(a) = sum_C I_q(a+1_C) for K4+K4 (q=3), K5-e via prism duality (q=3->4, 4->5), W5 (q=3)."""
import random
import time
from fractions import Fraction as F

import per5
import gutil
import per5
import gutil

t0 = time.time()
SK = per5.SKEL
W5 = SK['W5']
KK = SK['K4+K4']
K5e = SK['K5-e']
print('W5  q=2 ones:', per5.fmt(per5.period(W5, [1] * 10, 2)))
print('K4+K4 q=2 ones:', per5.fmt(per5.period(KK, [1] * 10, 2)))

# 2-sum rule on a graph that series-parallel rules also reduce: K4 minus (0,1) plus paths 0-4-1, 0-5-1
G = [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (0, 4), (4, 1), (0, 5), (5, 1)]
rng = random.Random(3)
nt = bad = 0
for q in (3, 4, 5):
    for _ in range(300):
        pw = [rng.randint(1, q - 1) for _ in G]
        if sum(pw) != 4 * q or not gutil.convergent(tuple(G), pw, q):
            continue
        a = per5.value(G, pw, q)
        b = per5.lower(G, pw, q)
        nt += 1
        bad += (a != b)
print('2-sum rule vs SP rules: %d vectors, %d mismatches' % (nt, bad))


def dimrel(E, q, nmax, label):
    cot = gutil.cotrees(tuple(E))
    N = len(E)
    n = gutil.cycle_rank(E)
    tested = bad = 0
    for c in gutil.compositions(n * (q - 1), N, 1, q - 2):
        if not gutil.convergent(tuple(E), list(c), q - 1):
            continue
        lhs = per5.value(E, list(c), q - 1)
        rhs = {}
        ok = True
        for C in cot:
            b = list(c)
            for e in C:
                b[e] += 1
            v = per5.value(E, b, q)
            if v is None:
                ok = False
                break
            rhs = per5.vadd(rhs, v)
        if not ok or lhs is None:
            continue
        tested += 1
        bad += (lhs != rhs)
        if tested >= nmax:
            break
    print('dimension relation %s q=%d->%d: %d vectors, %d mismatches' % (label, q - 1, q, tested, bad), flush=True)


dimrel(KK, 3, 50, 'K4+K4')
dimrel(K5e, 4, 200, 'K5-e')
dimrel(K5e, 5, 200, 'K5-e')
dimrel(W5, 3, 5, 'W5')
dimrel(W5, 4, 100, 'W5')
dimrel(SK['V7_1'], 4, 100, 'V7_1')
dimrel(SK['V8_0'], 4, 100, 'V8_0')
print('%.1fs' % (time.time() - t0))
