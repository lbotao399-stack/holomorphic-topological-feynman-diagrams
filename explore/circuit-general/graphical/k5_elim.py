"""K5 (six-loop core) by vertex elimination with function-valued response (graphical functions, d = 2q).
Part 1: the three engine values of symmetry/REPORT.md 3.3 (edge order 01,02,03,04,12,13,14,23,24,34).
Part 2: K5 weights with beta_{v,oo} = 0 for some v (T1) that the symmetry engine does not reach.
usage: python3 k5_elim.py ref | python3 k5_elim.py missed Q N SEED"""
import os
import random
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import gutil
import seq
import mzv
from validate import to_labels

E = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
REF = [
    (8, (5, 5, 6, 1, 6, 4, 5, 6, 5, 5), {'1': F(105943, 1569209241600), 'z3': F(-17, 302702400)}),
    (8, (4, 4, 6, 7, 6, 1, 5, 5, 5, 5), {'1': F(-62311, 224172748800), 'z3': F(1, 4324320)}),
    (7, (4, 2, 5, 5, 5, 6, 5, 2, 5, 3), {'1': F(3251321, 64665216000), 'z3': F(52643, 449064000),
                                         'z5': F(-23, 124740)}),
]


def beta_inf(alpha, q):
    deg = [0] * 5
    for (u, v), a in zip(E, alpha):
        deg[u] += q - a
        deg[v] += q - a
    return [2 * q - d for d in deg]


def run(q, a):
    t = time.time()
    val, pl = seq.period(E, list(a), q)
    return val, pl, time.time() - t


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'ref':
        for q, a, ref in REF:
            val, pl, dt = run(q, a)
            ok = to_labels(val) == ref
            print('I_%d(K5; %s) = %s   beta_oo=%s  [%s, %d appends, %.0fs]  engine value: %s' % (
                q, a, mzv.tostr(val), beta_inf(a, q), pl['route'], pl['appends'], dt, 'AGREE' if ok else 'DIFFER'),
                flush=True)
    else:
        q, N, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        sys.path.insert(2, os.path.join(HERE, '..', 'symmetry'))
        import engine as en
        rng = random.Random(seed)
        found = 0
        tried = 0
        while found < N and tried < 4000:
            tried += 1
            a = [rng.randint(1, q - 1) for _ in range(10)]
            if sum(a) != 6 * q or not gutil.convergent(tuple(E), a, q):
                continue
            if 0 not in beta_inf(a, q):
                continue
            if not seq.plans_for(E, a, q):
                continue
            r = en.period(E, a, q, budget=60)
            if r is not None:
                continue
            found += 1
            val, pl, dt = run(q, a)
            print('I_%d(K5; %s) = %s   beta_oo=%s  [%s, %d appends, %.0fs]  engine: not reached' % (
                q, tuple(a), mzv.tostr(val), beta_inf(a, q), pl['route'], pl['appends'], dt), flush=True)
        print('tried %d random vectors' % tried)
