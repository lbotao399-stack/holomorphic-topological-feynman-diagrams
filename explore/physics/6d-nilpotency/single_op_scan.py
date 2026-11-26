"""<Phi, hatQ_Gamma O> for single graph operators on jet-free (and one-jet) monomials.
Usage: python3 single_op_scan.py algebra nmax F seed"""
import sys
import itertools
import random

import numpy as np

import kern
from comp import Op, Functionals, evaluate1, so_n, su2
from kern import all_graphs

alg, nmax, F, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = random.Random(seed)
dim, f = su2() if alg == 'su2' else so_n(int(alg[2:]))
DEG = {'c': (0, 0), 'g1': (1, 0), 'g2': (0, 1), 'b': (1, 1)}


def finals(Ot, n):
    d1 = sum(DEG[t][0] for t in Ot) - n
    d2 = sum(DEG[t][1] for t in Ot) - n
    N = len(Ot) + 1
    out = []
    for ext in [('b',), ('g1', 'g2'), ('g1',), ('g2',), ()]:
        e1 = sum(DEG[t][0] for t in ext)
        e2 = sum(DEG[t][1] for t in ext)
        if (e1, e2) == (d1, d2) and len(ext) <= N:
            out.append(['c'] * (N - len(ext)) + list(ext))
    return out


for n in range(1, nmax + 1):
    gs = all_graphs(n)
    for Ot in itertools.combinations_with_replacement(['g1', 'g2', 'b'], n + 1):
        for jet in (None, 0):
            for ft in finals(Ot, n):
                O = tuple((('O', k), t) for k, t in enumerate(Ot))
                O_cols = {('O', k): np.array([rng.randint(-2, 2) for _ in range(dim)], dtype=np.int64)
                          for k in range(len(Ot))}
                func = Functionals(ft, dim, F, rng)
                jets = {('O', 0): (0,)} if jet is not None else {}
                vals = [evaluate1(Op(g, 'o'), O, O_cols, func, f, jets=jets) for g in gs]
                nz = sum(1 for v in vals for x in v if x)
                print('n=%d %s O=%s%s -> %s: nonzero <Phi, hatQ_Gamma O> %d of %d'
                      % (n, alg, ''.join(Ot), ' (d_1 on letter 0)' if jets else '', ''.join(ft), nz,
                         len(gs) * F), flush=True)
