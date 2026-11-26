"""q = 8 prism / K33 tables against the dimension relation  I_{q-1}(a) = sum_{C cotree} I_q(a + 1_C)
(not among the relations used to solve the systems [src cores/REPORT.md §1.3]), all q = 7 orbits."""
import time
import periods
import gutil

for name in ('prism', 'K33'):
    t = time.time()
    fe = periods.FAMILIES[name]
    t7 = periods.table(name, 7)
    t8 = periods.table(name, 8)
    cot = gutil.cotrees(tuple(fe))
    bad = n = 0
    for rep, val in t7.items():
        tot = {}
        for C in cot:
            a = list(rep)
            for e in C:
                a[e] += 1
            v = periods.period(fe, a, 8)
            for k, x in v.items():
                tot[k] = tot.get(k, 0) + x
        tot = {k: x for k, x in tot.items() if x}
        n += 1
        if tot != val:
            bad += 1
    print('%s: q=8 orbits %d; dimension relation q=7 -> 8 on %d orbits: %d failures (%.1fs)' % (name, len(t8), n, bad, time.time() - t), flush=True)
