"""Validation of the own operator engine (vops.py) against S2's statements:
  Q_0^2 = 0, {Q_0, Q_1} = 0 (S2 eq-res-nilp), and the order-hbar^2 relation
  Q_1^2 + {Q_0, Q_2} = 0 with Q_2 = (i/64pi^6) sum_X kappa_X hatQ^I_X, i.e.
  kappa_A <{hQ0,hQ^I_A}O> + kappa_B <{hQ0,hQ^I_B}O> = -(1/4) <hQ_box^2 O>, S2 eq-kAkB: kappa = (1/2, 1/2)."""
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from vops import Op, pair, so_n

dim, f = so_n(5)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)


def funcs(types):
    N = len(types)
    A = np.array([[rng.randint(-2, 2) for _ in range(dim)] for _ in range(N)], dtype=float)
    u = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(N)]
    return (list(types), A, u)


MODE = sys.argv[2] if len(sys.argv) > 2 else 'vec'


def mkO(types):
    O = [(('O', k), t) for k, t in enumerate(types)]
    if MODE == 'onehot':
        cols = {('O', k): c for k, c in zip(range(len(types)), rng.sample(range(dim), len(types)))}
    else:
        cols = {('O', k): [rng.randint(-2, 2) for _ in range(dim)] for k in range(len(types))}
    return O, cols


def anti(a, b, O, cols, jets, F, cache):
    """<F, (a b + b a) O> where 'a b' means b acts first."""
    return (pair([b, a], O, cols, jets, F, dim, f, cache) + pair([a, b], O, cols, jets, F, dim, f, cache))


cache = {}
t0 = time.time()
q0, p0 = Op('q0', 'q'), Op('q0', 'p')
box1, box2 = Op('box', 'u'), Op('box', 'd')
# Q0^2 = 0
for Ot, Ft in [(['b'], ['c', 'c', 'b']), (['b'], ['c', 'g1', 'g2']), (['g1', 'b'], ['c', 'c', 'g1', 'b']),
               (['c', 'g2'], ['c', 'c', 'c', 'g2'])]:
    O, cols = mkO(Ot)
    F = funcs(Ft)
    v = pair([q0, p0], O, cols, {}, F, dim, f, cache)
    print('Q0^2 %s, functionals %s: %s' % (Ot, Ft, v))
# {Q0, Q1} = 0
for Ot, Ft, jets in [(['g1', 'g2'], ['c'] * 4, {}), (['b', 'b'], ['c', 'c', 'c', 'b'], {}),
                     (['b', 'b'], ['c', 'c', 'g1', 'g2'], {}), (['g1', 'b'], ['c', 'c', 'c', 'g1'], {}),
                     (['g1', 'g2'], ['c'] * 4, {('O', 0): [0]}), (['b', 'g2'], ['c', 'c', 'c', 'g2'], {('O', 1): [2]})]:
    O, cols = mkO(Ot)
    F = funcs(Ft)
    v1 = pair([box1, q0], O, cols, jets, F, dim, f, cache)
    v2 = pair([q0, box1], O, cols, jets, F, dim, f, cache)
    print('{Q0,Q_box} %s jets %s, functionals %s: Q_box Q0 = %s, Q0 Q_box = %s, sum = %s' % (Ot, jets, Ft, v1, v2, v1 + v2))
print('[%.0fs]' % (time.time() - t0), flush=True)
# order hbar^2
A_, B_ = Op('A', 'u'), Op('B', 'u')
for Ot, Ft, jets, nf in [(['g1', 'g2', 'b'], ['c'] * 5, {}, 4), (['b', 'b', 'b'], ['c', 'c', 'c', 'g1', 'g2'], {}, 3),
                         (['b', 'b', 'g1'], ['c'] * 4 + ['g1'], {}, 3),
                         (['g1', 'g2', 'b'], ['c'] * 5, {('O', 0): [0]}, 3),
                         (['b', 'b', 'g2'], ['c'] * 4 + ['g2'], {('O', 2): [1]}, 2)]:
    rows = []
    t1 = time.time()
    O, cols = mkO(Ot)
    for _ in range(nf):
        F = funcs(Ft)
        eA = anti(q0, Op('A', 'd'), O, cols, jets, F, cache)
        eB = anti(q0, Op('B', 'd'), O, cols, jets, F, cache)
        bb = pair([box1, box2], O, cols, jets, F, dim, f, cache)
        rows.append((eA, eB, -Fr(1, 4) * bb))
        print('   %s jets %s: <{Q0,A}> = %s, <{Q0,B}> = %s, -(1/4)<box^2> = %s' % (Ot, jets, eA, eB, -Fr(1, 4) * bb), flush=True)
    M = sp.Matrix([[sp.Rational(r[0].numerator, r[0].denominator), sp.Rational(r[1].numerator, r[1].denominator)] for r in rows])
    b = sp.Matrix([sp.Rational(r[2].numerator, r[2].denominator) for r in rows])
    aug = M.row_join(b)
    sol = sp.linsolve((M, b))
    print('n=2, O = %s, jets %s: %d equations, rank %d, augmented rank %d, solution (kappa_A, kappa_B) = %s  [%.0fs]'
          % (Ot, jets, len(rows), M.rank(), aug.rank(), sol, time.time() - t1), flush=True)
print('total [%.0fs]' % (time.time() - t0))
