"""Relations with one input derivative (displacement maps at the level of the operators):
  n = 2:  O = d_1 gamma^1 gamma^2 b:      kappa_A <{Q0,Q^I_A}O> + kappa_B <{Q0,Q^I_B}O> = -(1/4) <Q_box^2 O>
  n = 3:  O = d_1 gamma^1 gamma^2 b b:    sum kappa_Gamma <{Q0,Q^I_Gamma}O> = -(1/4) <{Q_box, Q^I_A + Q^I_B}O>
(normalizations as in q3check.py), so(5), random colour vectors, plane-wave functionals."""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import opcheck as oc

part = sys.argv[1] if len(sys.argv) > 1 else 'both'
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 8
n, f, fd = oc.fdense('so5')
rng = np.random.RandomState(13)
d1 = (1, 0, 0)
z0 = (0, 0, 0)


def solve(cols, rhs_vals, names):
    ks = sp.symbols('k1:%d' % (len(names) + 1))
    M = sp.Matrix([[sp.Rational(c[i].numerator, c[i].denominator) for c in cols] for i in range(len(rhs_vals))])
    rhs = sp.Matrix([-sp.Rational(x.numerator, x.denominator) / 4 for x in rhs_vals])
    print('   rank = %d, rank augmented = %d, solution %s = %s'
          % (M.rank(), M.row_join(rhs).rank(), tuple(names), sp.linsolve((M, rhs), *ks)), flush=True)


if part in ('both', '2'):
    t0 = time.time()
    Ov = [rng.randint(-2, 3, size=n).astype(np.int64) for _ in range(3)]
    tds = [oc.TestDataV(n, 900 + s, n=5) for s in range(6)]
    types, jets = [3, 4, 2], [d1, z0, z0]
    E11 = oc.evaluate4(['box', 'box'], types, jets, f, Ov, tds)[0]
    EA = [a + b for a, b in zip(oc.evaluate5(['A', 'Q0'], types, jets, f, Ov, tds)[0], oc.evaluate5(['Q0', 'A'], types, jets, f, Ov, tds)[0])]
    EB = [a + b for a, b in zip(oc.evaluate5(['B', 'Q0'], types, jets, f, Ov, tds)[0], oc.evaluate5(['Q0', 'B'], types, jets, f, Ov, tds)[0])]
    print('n = 2, O = d_1 gamma^1 gamma^2 b, %d equations  [%.0fs]' % (len(tds), time.time() - t0))
    for a, b, c in zip(EA, EB, E11):
        print('   %s  %s  | Q_box^2: %s' % (a, b, c))
    solve([EA, EB], E11, ['kappa_A', 'kappa_B'])

if part in ('both', '3'):
    T0 = time.time()
    Ov = [rng.randint(-2, 3, size=n).astype(np.int64) for _ in range(4)]
    tds = [oc.TestDataV(n, 950 + s, n=6) for s in range(NT)]
    types, jets = [3, 4, 2, 2], [d1, z0, z0, z0]
    E = {}
    for g in ['G1', 'G2', 'G3', 'G4', 'G5']:
        t0 = time.time()
        a, h1 = oc.evaluate5([g, 'Q0'], types, jets, f, Ov, tds)
        b, h2 = oc.evaluate5(['Q0', g], types, jets, f, Ov, tds)
        E[g] = [x + y for x, y in zip(a, b)]
        print('%s: %d + %d histories [%.0fs]' % (g, h1, h2, time.time() - t0), flush=True)
    Rhs = [Fr(0)] * len(tds)
    for a, b in [('box', 'A'), ('A', 'box'), ('box', 'B'), ('B', 'box')]:
        t0 = time.time()
        v, h = oc.evaluate4([a, b], types, jets, f, Ov, tds)
        Rhs = [x + y for x, y in zip(Rhs, v)]
        print('%s*%s: %d histories [%.0fs]' % (a, b, h, time.time() - t0), flush=True)
    print('n = 3, O = d_1 gamma^1 gamma^2 b b, %d equations' % len(tds))
    for i in range(len(tds)):
        print('   ' + '  '.join(str(E[g][i]) for g in ['G1', 'G2', 'G3', 'G4', 'G5']) + '  | -(1/4) rhs = %s' % (-Rhs[i] / 4))
    solve([E[g] for g in ['G1', 'G2', 'G3', 'G4', 'G5']], Rhs, ['G1', 'G2', 'G3', 'G4', 'G5'])
    print('total %.0fs' % (time.time() - T0))
