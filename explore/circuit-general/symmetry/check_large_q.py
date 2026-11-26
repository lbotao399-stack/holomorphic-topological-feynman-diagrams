"""Reach at large q (hyperplane conditions of measure O(1/q)): coverage samples at q = 20, 40;
targets W_n and SP only (prism/K33 not tabulated there)."""
import random
import time
from collections import Counter

import engine as en
import gutil
from coverage import SK, sample

en.NUMER[0] = 3
for q in (20, 40):
    for name in ['K5-e', 'W5', 'K33+e', 'Z5', 'K4+K4', 'V7_0', 'V7_2', 'V8_3', 'V8_4', 'K5']:
        G = SK[name]
        n = gutil.cycle_rank(G)
        rng = random.Random(99)
        acc = hit = 0
        t = time.time()
        while acc < 60:
            a = sample(n * q, len(G), q - 1, rng)
            if not gutil.convergent(tuple(G), a, q):
                continue
            acc += 1
            hit += en.period(G, a, q, budget=60) is not None
        print('%-6s q=%d reached %d/%d  (%.1fs)' % (name, q, hit, acc, time.time() - t), flush=True)
