"""q = Q prism / K33 tables (usage: python3 check_tables_dim.py Q [names]) against the dimension relation  I_{q-1}(a) = sum_{C cotree} I_q(a + 1_C)
(not among the relations used to solve the systems [src cores/REPORT.md §1.3]), all q = Q-1 orbits."""
import time
import periods
import gutil

import sys
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 8
names = sys.argv[2:] or ["prism", "K33"]
for name in names:
    t = time.time()
    fe = periods.FAMILIES[name]
    t7 = periods.table(name, Q - 1)
    t8 = periods.table(name, Q)
    cot = gutil.cotrees(tuple(fe))
    bad = n = 0
    for rep, val in t7.items():
        tot = {}
        for C in cot:
            a = list(rep)
            for e in C:
                a[e] += 1
            v = periods.period(fe, a, Q)
            for k, x in v.items():
                tot[k] = tot.get(k, 0) + x
        tot = {k: x for k, x in tot.items() if x}
        n += 1
        if tot != val:
            bad += 1
    print('%s: q=%d orbits %d; dimension relation q=%d -> %d on %d orbits: %d failures (%.1fs)' % (name, Q, len(t8), Q - 1, Q, n, bad, time.time() - t), flush=True)
