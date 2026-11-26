"""Product identity at a 3-vertex cut {a, b, inf} = 2-vertex cut {a, b} of the decompleted graph.
Independent check on K4+K4 (2-sum of two K4 without the common edge, five loops):
values from the regulator-free relation system of cores/msys.py (triangle rules, IBP, boundaries,
star-triangle; boundary data from wheels and SP) versus the engine (product split -> two completed K4)."""
import sys
import time

import engine as en
import msys
import cores

K4K4 = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5))

for q in [int(x) for x in sys.argv[1:]] or [3, 4]:
    t = time.time()
    S = msys.MultiSystem([K4K4], q)
    S.lower = cores.lower_core_value
    S.build(verbose=True)
    solved, piv = S.solve(verbose=True)
    ok = bad = und = 0
    for j, (f, rep) in enumerate(S.reps):
        r = en.period(K4K4, rep, q, use_m2=False)
        if r is None:
            und += 1
            continue
        if j not in solved:
            continue
        ref = en.from_cores(solved[j])
        if r[0] == ref:
            ok += 1
        else:
            bad += 1
            print('MISMATCH', rep, en.vstr(r[0]), '|', en.vstr(ref), flush=True)
    ex = S.reps[sorted(solved)[0]][1] if solved else None
    print('K4+K4 q=%d: orbits %d, msys-determined %d, engine-evaluated %d, compared %d, mismatches %d  (%.1fs)'
          % (q, len(S.reps), len(solved), len(S.reps) - und, ok + bad, bad, time.time() - t), flush=True)
    if ex:
        r = en.period(K4K4, ex, q, use_m2=False)
        print('  example', ex, en.vstr(r[0]), r[1], flush=True)
