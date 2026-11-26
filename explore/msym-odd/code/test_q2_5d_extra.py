"""Further operators for {Q_0, Qhat_A} and {Q_0, Qhat_B} (5d): 3 and 4 letters, up to two input
derivatives.  Usage: python3 test_q2_5d_extra.py so3|so5 [max_tests] [seed]"""
import itertools
import random
import sys
import time

import sympy as sp

from salg import structure_constants, add, mul
from ops5 import Theory, L, KernelData, TwoLoop

name = sys.argv[1] if len(sys.argv) > 1 else 'so3'
maxt = int(sys.argv[2]) if len(sys.argv) > 2 else 40
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
n, f = structure_constants(name)
T = Theory(n, f)
K = {k: KernelData(k) for k in ('A', 'B')}
T2 = TwoLoop(T, K)
rng = random.Random(seed)
types = ['c', 'g1', 'g2', 'b']
jets = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0), (1, 1, 0), (0, 2, 0)]
tests = []
seen = set()
while len(tests) < maxt:
    k = rng.choice([3, 3, 4])
    lets = []
    nder = 0
    for _ in range(k):
        j = rng.choice(jets[:3] if nder >= 1 else jets)
        nder += sum(j)
        lets.append((rng.choice(types), rng.randrange(n), j))
    if nder > 2:
        continue
    key = tuple(sorted(lets))
    if key in seen:
        continue
    seen.add(key)
    p = {((), ()): 1}
    for (tt, A, j) in lets:
        p = mul(p, L(tt, A, j))
    if not p:
        continue
    tests.append((' '.join('%s%s_%d' % ('d' * sum(j) + ''.join(str(i + 1) * j[i] for i in range(2)) if sum(j) else '', tt, A) for tt, A, j in lets), p))
tot_eq, nonzero = 0, 0
rows = []
t0 = time.time()
for nm, p in tests:
    t1 = time.time()
    a = add(T.Q0(T2.Q2('A', p)), T2.Q2('A', T.Q0(p)))
    b = add(T.Q0(T2.Q2('B', p)), T2.Q2('B', T.Q0(p)))
    qa = T2.Q2('A', p)
    eqs = [(sp.Rational(a.get(m, 0)), sp.Rational(b.get(m, 0))) for m in set(a) | set(b)]
    tot_eq += len(eqs)
    nonzero += bool(eqs)
    rows += eqs
    print('%-28s |QA O| = %5d  |{Q0,QA}O| = %5d  |{Q0,QB}O| = %5d  (%.1fs)' % (nm, len(qa), len(a), len(b), time.time() - t1), flush=True)
rk = sp.Matrix(rows).rank() if rows else 0
print('%s: %d operators, %d with nonzero anticommutators, %d equations, rank %d, null space %s  (%.1fs)'
      % (name, len(tests), nonzero, tot_eq, rk, [list(v) for v in sp.Matrix(rows).nullspace()] if rows else 'all', time.time() - t0))
