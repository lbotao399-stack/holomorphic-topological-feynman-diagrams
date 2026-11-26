"""<Phi, hatQ_Gamma O> for the four-loop graph operators on su(2) monomials with jets (cochain level).
Usage: python3 single_op_jets.py algebra O_types jets F seed     (jets as in run_nilp: k:a.b;k2:c)"""
import sys
import random

import numpy as np

from comp import Op, Functionals, evaluate1, so_n, su2
from kern import all_graphs

alg, Os, js, F, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
dim, f = su2() if alg == 'su2' else so_n(int(alg[2:]))
Ot = Os.split(',')
jets = {}
for item in js.split(';'):
    if item:
        k, comps = item.split(':')
        jets[('O', int(k))] = tuple(int(c) for c in comps.split('.'))
O = tuple((('O', k), t) for k, t in enumerate(Ot))
O_cols = {('O', k): np.array([rng.randint(-2, 2) for _ in range(dim)], dtype=np.int64) for k in range(len(Ot))}
DEG = {'c': (0, 0), 'g1': (1, 0), 'g2': (0, 1), 'b': (1, 1)}
d = (sum(DEG[t][0] for t in Ot) - 4, sum(DEG[t][1] for t in Ot) - 4)
ext = {(0, 0): [], (1, 0): ['g1'], (0, 1): ['g2'], (1, 1): ['b']}[d]
ft = ['c'] * (len(Ot) + 1 - len(ext)) + ext
func = Functionals(ft, dim, F, rng)
gs = all_graphs(4)
vals = [evaluate1(Op(g, 'o'), O, O_cols, func, f, jets=jets) for g in gs]
print('%s O=%s jets %s -> %s: nonzero graphs %s' % (alg, ''.join(Ot), jets, ''.join(ft),
      [g.name for g, v in zip(gs, vals) if any(int(x) for x in v)]), flush=True)
