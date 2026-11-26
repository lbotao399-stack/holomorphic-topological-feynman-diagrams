"""count S6 classes of convergent balanced K6 weights at level q; kinds T1/T2/G; K5-decompletable"""
import sys, time
import k6
for q in map(int, sys.argv[1:]):
    t = time.time()
    n = 0
    cl = k6.classes(q)
    tot = sum(cl.values())
    kinds = {}
    for c in cl:
        kd = k6.kind(c, q)
        dec = bool(k6.k5_decompletions(c, q))
        pos = min(c) >= 1
        kinds[(kd, dec, pos)] = kinds.get((kd, dec, pos), 0) + 1
    print('q=%d labelled %d classes %d  (%.1fs)' % (q, tot, len(cl), time.time() - t))
    for k in sorted(kinds):
        print('   kind %s  K5-decompletable %s  all beta>=1 %s : %d' % (k[0], k[1], k[2], kinds[k]))
