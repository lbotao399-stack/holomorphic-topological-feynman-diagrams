"""<Phi, hatQ_Gamma O> for the 15 four-loop graph operators on jet-free so(5) / su(2) monomials (cochain level).
Usage: python3 single_op_n4.py algebra F seed"""
import sys
import random

import numpy as np

from comp import Op, Functionals, evaluate1, so_n, su2
from kern import all_graphs, ratrec

alg, F, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed)
dim, f = su2() if alg == 'su2' else so_n(int(alg[2:]))
gs = all_graphs(4)
cases = [('b,b,b,b,b', 'c,c,c,c,c,b'), ('b,b,b,b,b', 'c,c,c,c,g1,g2'), ('g1,g2,b,b,b', 'c,c,c,c,c,c'),
         ('g1,b,b,b,b', 'c,c,c,c,c,g1'), ('b,b,b,b,c', 'c,c,c,c,c,c')]
for Os, fs in cases:
    Ot, ft = Os.split(','), fs.split(',')
    O = tuple((('O', k), t) for k, t in enumerate(Ot))
    O_cols = {('O', k): np.array([rng.randint(-2, 2) for _ in range(dim)], dtype=np.int64) for k in range(len(Ot))}
    func = Functionals(ft, dim, F, rng)
    vals = [evaluate1(Op(g, 'o'), O, O_cols, func, f) for g in gs]
    print('%s O=%s -> %s: nonzero graphs %s' % (alg, ''.join(Ot), ''.join(ft),
          [g.name for g, v in zip(gs, vals) if any(int(x) for x in v)]), flush=True)
