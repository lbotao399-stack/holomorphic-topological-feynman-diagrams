"""Order hbar^3 of Q^2 = 0:  {Q_0, Q_3} + {Q_1, Q_2} = 0 on four-letter operators O (six final c letters).
On three-letter operators every term vanishes separately (an operator consuming two letters
produced at one theta' gives [(theta_1 - theta_1')(theta_2 - theta_2')]^2 = 0).

Normalizations (S2 (5.15), (7.15)):  Q_0 = i hat Q_0,  Q_1 = (i/8pi^3)(1/2) hat Q^I_box,
Q_2 = (i/64pi^6)(1/2)(hat Q^I_A + hat Q^I_B),  ansatz  Q_3 = (i/512pi^9) sum_Gamma kappa_Gamma hat Q^I_Gamma,
hat Q^I_Gamma: operator of opcheck.py with the signed kernel I_Gamma.  Then
   sum_Gamma kappa_Gamma <{hat Q_0, hat Q^I_Gamma} O, Phi> = -(1/4) <{hat Q^I_box, hat Q^I_A + hat Q^I_B} O, Phi>.
Prediction of S2 (5.15) with W_Gamma of the leg pattern of S2 (1.4), (1.6), (1.7):  kappa_Gamma = 1/|Aut Gamma|.
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import opcheck as oc

alg = sys.argv[1] if len(sys.argv) > 1 else 'so5'
NO = int(sys.argv[2]) if len(sys.argv) > 2 else 2
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 6
n, f, fd = oc.fdense(alg)
rng = np.random.RandomState(11)
TYPES = {'g1g2bb': [3, 4, 2, 2], 'bbbc': [2, 2, 2, 1]}
oname = sys.argv[4] if len(sys.argv) > 4 else 'g1g2bb'
O_types = TYPES[oname]
Ovs = [[rng.randint(-2, 3, size=n).astype(np.int64) for _ in range(len(O_types))] for _ in range(NO)]
tds = [oc.TestDataV(n, 500 + s, n=6) for s in range(NT)]

T0 = time.time()
E = {}
for g in ['G1', 'G2', 'G3', 'G4', 'G5']:
    t0 = time.time()
    v1 = []; v2 = []
    for Ov in Ovs:
        x, h1 = oc.evaluate3([g, 'Q0'], O_types, f, Ov, tds); v1 += x
        x, h2 = oc.evaluate3(['Q0', g], O_types, f, Ov, tds); v2 += x
    E[g] = [a + b for a, b in zip(v1, v2)]
    print('%s: histories %d (Q_%s Q_0), %d (Q_0 Q_%s)  [%.0fs]' % (g, h1, g, h2, g, time.time() - t0), flush=True)
Rhs = None
for a, b in [('box', 'A'), ('A', 'box'), ('box', 'B'), ('B', 'box')]:
    t0 = time.time()
    flat = []
    for Ov in Ovs:
        x, h = oc.evaluate3([a, b], O_types, f, Ov, tds); flat += x
    Rhs = flat if Rhs is None else [x + y for x, y in zip(Rhs, flat)]
    print('%s*%s: histories %d  [%.0fs]' % (a, b, h, time.time() - t0), flush=True)

ks = sp.symbols('k1:6')
M = sp.Matrix([[sp.Rational(E[g][i].numerator, E[g][i].denominator) for g in ['G1', 'G2', 'G3', 'G4', 'G5']]
               for i in range(len(Rhs))])
rhs = sp.Matrix([-sp.Rational(x.numerator, x.denominator) / 4 for x in Rhs])
print('\nalgebra %s, O = %s with %d random colour assignments, %d plane-wave functionals each: %d equations' % (alg, oname, NO, NT, len(Rhs)))
for i in range(len(Rhs)):
    print('   ' + '  '.join('%s' % E[g][i] for g in ['G1', 'G2', 'G3', 'G4', 'G5']) + '   | -(1/4) rhs = %s' % (-Rhs[i] / 4))
print('rank of [E_Gamma] = %d, rank augmented = %d' % (M.rank(), M.row_join(rhs).rank()))
sol = sp.linsolve((M, rhs), *ks)
print('solution (kappa_G1, ..., kappa_G5) =', sol)
print('prediction 1/|Aut| =', tuple(sp.Rational(1, oc.AUT[g]) for g in ['G1', 'G2', 'G3', 'G4', 'G5']))
pred = sp.Matrix([sp.Rational(1, oc.AUT[g]) for g in ['G1', 'G2', 'G3', 'G4', 'G5']])
print('residual of the prediction:', list(M * pred - rhs))
print('total time %.0fs' % (time.time() - T0))
