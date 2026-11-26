"""C7, C8, C10 with the own engine (vops.py): order hbar^3 of Q^2 = 0,
   sum_Gamma kappa_Gamma <{hQ0, hQ^I_Gamma} O> = -(1/4) <{hQ^I_box, hQ^I_A + hQ^I_B} O>,
on test operators different from the explorer's (so(5), colour vectors in the operator letters).
Usage: python3 v_nilp3.py MODE SEED
  MODE = main   : O = bbbb, bbbg1, bbbg2 (no jets) and d_2 b bbb, b bb d_3 g1 (one jet)
         c8     : O = bbb, every term separately
         c10    : O = g1g2bb and bbbc with the same colour vectors and the same functionals
         so6    : O = bbbb in so(6)
"""
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from vops import Op, pair, so_n

MODE = sys.argv[1] if len(sys.argv) > 1 else 'main'
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 3
rng = random.Random(SEED)
dim, f = so_n(6 if MODE == 'so6' else 5)
GS = ['G1', 'G2', 'G3', 'G4', 'G5']
KAPPA_EXPL = [Fr(-1, 2), Fr(1, 2), Fr(-1), Fr(1, 8), Fr(-1, 2)]
cache = {}
q0 = Op('q0', 'q')


def funcs(types):
    N = len(types)
    A = np.array([[rng.randint(-2, 2) for _ in range(dim)] for _ in range(N)], dtype=float)
    u = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(N)]
    return (list(types), A, u)


def mkO(types, cols=None):
    O = [(('O', k), t) for k, t in enumerate(types)]
    if cols is None:
        cols = [[rng.randint(-2, 2) for _ in range(dim)] for _ in types]
    return O, {('O', k): c for k, c in enumerate(cols)}


def terms(O, cols, jets, F):
    """returns dict of the separate terms."""
    out = {}
    for g in GS:
        out[(g, 'G.Q0')] = pair([q0, Op(g, 'd')], O, cols, jets, F, dim, f, cache)     # hQ_G hQ0 O
        out[(g, 'Q0.G')] = pair([Op(g, 'u'), q0], O, cols, jets, F, dim, f, cache)     # hQ0 hQ_G O
    for X in ['A', 'B']:
        out[('box.' + X)] = pair([Op(X, 'u'), Op('box', 'd')], O, cols, jets, F, dim, f, cache)
        out[(X + '.box')] = pair([Op('box', 'u'), Op(X, 'd')], O, cols, jets, F, dim, f, cache)
    return out


def row(T):
    E = [T[(g, 'G.Q0')] + T[(g, 'Q0.G')] for g in GS]
    rhs = -Fr(1, 4) * (T['box.A'] + T['A.box'] + T['box.B'] + T['B.box'])
    return E, rhs


def solve(rows, label):
    M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in E] for E, _ in rows])
    b = sp.Matrix([sp.Rational(r.numerator, r.denominator) for _, r in rows])
    sol = sp.linsolve((M, b), *sp.symbols('k1:6'))
    res = [sum(E[i] * KAPPA_EXPL[i] for i in range(5)) - r for E, r in rows]
    print('%s: %d equations, rank %d, augmented rank %d, solution (kappa_G1..G5) = %s; '
          'residual of the explorer kappa = (-1/2, 1/2, -1, 1/8, -1/2): %s'
          % (label, len(rows), M.rank(), M.row_join(b).rank(), sol, [str(x) for x in res]), flush=True)
    return M, b


T00 = time.time()
if MODE in ('main', 'so6'):
    cases = [(['b', 'b', 'b', 'b'], {}, [['c', 'c', 'c', 'c', 'g1', 'g2'], ['c', 'c', 'c', 'c', 'c', 'b']], 3),
             (['b', 'b', 'b', 'g1'], {}, [['c'] * 5 + ['g1']], 4),
             (['b', 'b', 'b', 'g2'], {}, [['c'] * 5 + ['g2']], 3),
             (['b', 'b', 'b', 'b'], {('O', 0): [1]}, [['c', 'c', 'c', 'c', 'g1', 'g2']], 4),
             (['b', 'b', 'b', 'g1'], {('O', 3): [2]}, [['c'] * 5 + ['g1']], 3)]
    if MODE == 'so6':
        cases = cases[:1]
    allrows = []
    for Ot, jets, Fts, nf in cases:
        t0 = time.time()
        O, cols = mkO(Ot)
        rows = []
        for Ft in Fts:
            for _ in range(nf):
                F = funcs(Ft)
                T = terms(O, cols, jets, F)
                E, rhs = row(T)
                rows.append((E, rhs))
                print('   O=%s jets=%s functionals %s:  E = %s | -(1/4)rhs = %s'
                      % (''.join(Ot), jets, ''.join(Ft), [str(x) for x in E], rhs), flush=True)
        solve(rows, 'so(%d) O = %s, jets %s  [%.0fs]' % (6 if MODE == 'so6' else 5, ''.join(Ot), jets, time.time() - t0))
        allrows += rows
    solve(allrows, 'all operators combined')
elif MODE == 'c8':
    O, cols = mkO(['b', 'b', 'b'])
    for k in range(3):
        F = funcs(['c'] * 5)
        T = terms(O, cols, {}, F)
        print('O = bbb, functionals ccccc #%d: ' % k + ', '.join('%s: %s' % (key, v) for key, v in T.items()), flush=True)
    for Ft in (['c', 'c', 'c', 'c', 'b'], ['c', 'c', 'c', 'g1', 'g2']):
        F = funcs(Ft)
        T = terms(O, cols, {}, F)
        print('O = bbb, functionals %s (type mismatch control): ' % ''.join(Ft) + ', '.join('%s' % v for v in T.values()))
elif MODE == 'c10':
    w = [[rng.randint(-2, 2) for _ in range(dim)] for _ in range(4)]
    O1, c1 = mkO(['g1', 'g2', 'b', 'b'], w)
    O2, c2 = mkO(['b', 'b', 'b', 'c'], w)
    rows1, rows2 = [], []
    for k in range(int(sys.argv[3]) if len(sys.argv) > 3 else 3):
        t0 = time.time()
        F = funcs(['c'] * 6)
        T1 = terms(O1, c1, {}, F)
        T2 = terms(O2, c2, {}, F)
        r1, r2 = row(T1), row(T2)
        rows1.append(r1)
        rows2.append(r2)
        print('functional set %d: g1g2bb E = %s | -(1/4)rhs = %s' % (k, [str(x) for x in r1[0]], r1[1]))
        print('                 bbbc   E = %s | -(1/4)rhs = %s' % ([str(x) for x in r2[0]], r2[1]))
        print('                 term by term T(bbbc) == -T(g1g2bb): %s;  == +T: %s  [%.0fs]'
              % (all(T2[key] == -T1[key] for key in T1), all(T2[key] == T1[key] for key in T1), time.time() - t0),
              flush=True)
    solve(rows1, 'O = g1g2bb')
    solve(rows2, 'O = bbbc')
print('total %.0fs' % (time.time() - T00))
