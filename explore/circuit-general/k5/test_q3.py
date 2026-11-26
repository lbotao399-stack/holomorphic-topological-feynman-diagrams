"""level-3 (and level-2) K6 classes: values by engine / per5 / graphical functions; consistency"""
import pickle, time, sys
import k6, bval
for q in (2, 3):
    from solve_level import classes; cl = classes(q)
    for b in sorted(cl):
        out = []
        for m in ('E', 'P', 'S'):
            t = time.time()
            r = bval.value(b, q, methods=(m,))
            out.append((m, r, time.time() - t))
        vals = [r[0] for m, r, dt in out if r is not None]
        agree = all(v == vals[0] for v in vals)
        print(q, b, k6.kind(b, q), '|', ' ; '.join('%s:%s (%.1fs)' % (m, (r[0] if r else None), dt) for m, r, dt in out),
              '| agree' if agree else '| DISAGREE', flush=True)
