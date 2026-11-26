"""Route independence: for sampled convergent weight vectors on the five-loop skeletons and K5,
collect every chain the search finds (all_routes, budget 60, two-point terminal k <= 3) and compare the
values; distinct routes use different identities and different tabulated data (W_n, prism, K33, SP)."""
import random
import sys
import time
from collections import Counter

import engine as en
import gutil
from coverage import SK, ORDER, sample

en.NUMER[0] = 3
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
for q in [int(x) for x in sys.argv[2:]] or [5, 7]:
    for name in ORDER:
        G = SK[name]
        n = gutil.cycle_rank(G)
        rng = random.Random(777 + q)
        t = time.time()
        acc = multi = bad = 0
        terms = Counter()
        while acc < N:
            a = sample(n * q, len(G), q - 1, rng)
            if not gutil.convergent(tuple(G), a, q):
                continue
            acc += 1
            V, E, inf = en.complete(G, a, q)
            rs = en.solve(V, E, q, budget=60, all_routes=True)
            if len(rs) < 2:
                continue
            vals = {tuple(sorted(r[0].items())) for r in rs}
            kinds = {r[1][-1].split(':')[-1] for r in rs}
            multi += 1
            terms[len(kinds)] += 1
            if len(vals) > 1:
                bad += 1
                print('MISMATCH', name, q, a, [(en.vstr(r[0]), r[1]) for r in rs[:3]], flush=True)
        print('%-6s q=%d samples %d  with >=2 routes %d (distinct terminal types: %s)  mismatches %d  (%.1fs)'
              % (name, q, acc, multi, dict(terms), bad, time.time() - t), flush=True)
