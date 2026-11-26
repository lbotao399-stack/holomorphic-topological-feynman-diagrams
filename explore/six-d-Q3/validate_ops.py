"""Validation of the slot evaluation of opcheck.py against the operators of S2's code
(s2code/ops.py, ops2.py, copied verbatim), and the two-loop relation of S2 (7.16)-(7.17).

Cases (O with fixed colours = unit vectors):
  V1  <hat Q_box (gamma^1_0 gamma^2_1), Phi_3>            su(2), so(5)
  V2  <hat Q_A (gamma^1_0 gamma^2_1 b_2), Phi_4>, <hat Q_B ...>   so(5)   (hat Q_B of S2 has J_B = -I_B)
  V3  <hat Q_box hat Q_box O>, <hat Q_A hat Q_0 O>, <hat Q_0 hat Q_A O>, same for B, O = gamma^1 gamma^2 b, su(2)
  V4  sum_Gamma kappa_Gamma <{hat Q_0, hat Q^I_Gamma} O> = -(1/4) <hat Q_box^2 O>, so(5), solve for kappa_A, kappa_B
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import opcheck as oc
from salg import structure_constants, mul, add
from ops import Theory, L
from ops2 import TwoLoop


def pair(P, td):
    tot = Fr(0)
    for (odd, even), c in P.items():
        assert not even and all(l[0] == 1 for l in odd) and len(odd) == td.n, (odd, even)
        for sig in itertools.permutations(range(td.n)):
            val = 1
            for l, i in zip(odd, sig):
                A, k = l[1], l[2:]
                val *= int(td.a[i][A]) * td.u[i][0] ** k[0] * td.u[i][1] ** k[1] * td.u[i][2] ** k[2]
            if val:
                tot += oc.perm_sign(sig) * val * c
    return tot


def unit(n, A):
    e = np.zeros(n, dtype=np.int64)
    e[A] = 1
    return e


def _main():
  pass

for alg in ([] if __name__ != "__main__" else ["su2", "so5"]):
    n, f, fd = oc.fdense(alg)
    T = Theory(n, fd)
    T2 = TwoLoop(T)
    # V1
    for seed in range(3):
        td = oc.TestData(n, seed, n=3)
        p = mul(L('g1', 0), L('g2', 1))
        s2 = pair(T.Q1(p), td)
        sk, nh = oc.evaluate(['box'], [3, 4], f, [[unit(n, 0), unit(n, 1)]], [td])
        print('V1 %s seed %d: S2 code %s, slot evaluation %s, equal %s' % (alg, seed, s2, sk[0], s2 == sk[0]), flush=True)
    if alg == 'so5':
        p = mul(mul(L('g1', 0), L('g2', 1)), L('b', 2))
        for kind in ['A', 'B']:
            t0 = time.time()
            P2 = T2.Q2(kind, p)
            t1 = time.time() - t0
            for seed in range(2):
                td = oc.TestData(n, 10 + seed, n=4)
                s2 = pair(P2, td)
                sk, nh = oc.evaluate([kind], [3, 4, 2], f, [[unit(n, 0), unit(n, 1), unit(n, 2)]], [td])
                rel = 1 if kind == 'A' else -1
                print('V2 so5 %s seed %d: S2 code %s, slot evaluation (kernel I) %s, I/J = %+d: %s  [S2 code %.1fs]'
                      % (kind, seed, s2, sk[0], rel, sk[0] == rel * s2, t1), flush=True)
    if alg == 'su2':
        p = mul(mul(L('g1', 0), L('g2', 1)), L('b', 2))
        Ov = [[unit(n, 0), unit(n, 1), unit(n, 2)]]
        t0 = time.time()
        Q11 = T.Q1(T.Q1(p))
        QA0 = T2.Q2('A', T.Q0(p)); Q0A = T.Q0(T2.Q2('A', p))
        QB0 = T2.Q2('B', T.Q0(p)); Q0B = T.Q0(T2.Q2('B', p))
        print('   S2 code compositions (su2): %.1fs; sizes %d %d %d %d %d' % (time.time() - t0, len(Q11), len(QA0), len(Q0A), len(QB0), len(Q0B)))
        for seed in range(2):
            td = oc.TestData(n, 20 + seed, n=5)
            rows = []
            for seq, P, rel in [(['box', 'box'], Q11, 1), (['A', 'Q0'], QA0, 1), (['Q0', 'A'], Q0A, 1),
                                (['B', 'Q0'], QB0, -1), (['Q0', 'B'], Q0B, -1)]:
                s2 = pair(P, td)
                sk, nh = oc.evaluate(seq, [3, 4, 2], f, Ov, [td])
                print('V3 su2 seed %d %-12s: S2 code %s, slot evaluation %s, equal (times %+d): %s  [%d histories]'
                      % (seed, '*'.join(seq), s2, sk[0], rel, sk[0] == rel * s2, nh), flush=True)

if __name__ == "__main__":
    # V4: two-loop relation, so(5), O = gamma^1 gamma^2 b with random colour vectors, several test data
    n, f, fd = oc.fdense('so5')
    rng = np.random.RandomState(3)
    tds, Ovs = [], []
    for seed in range(8):
        tds.append(oc.TestData(n, 100 + seed, n=5))
        Ovs.append([rng.randint(-2, 3, size=n).astype(np.int64) for _ in range(3)])
    t0 = time.time()
    E11, h1 = oc.evaluate(['box', 'box'], [3, 4, 2], f, Ovs, tds)
    EA = [a + b for a, b in zip(oc.evaluate(['A', 'Q0'], [3, 4, 2], f, Ovs, tds)[0], oc.evaluate(['Q0', 'A'], [3, 4, 2], f, Ovs, tds)[0])]
    EB = [a + b for a, b in zip(oc.evaluate(['B', 'Q0'], [3, 4, 2], f, Ovs, tds)[0], oc.evaluate(['Q0', 'B'], [3, 4, 2], f, Ovs, tds)[0])]
    kA, kB = sp.symbols('kA kB')
    M = sp.Matrix([[sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)] for a, b in zip(EA, EB)])
    rhs = sp.Matrix([-sp.Rational(e.numerator, e.denominator) / 4 for e in E11])
    print('V4 so5, O = gamma^1 gamma^2 b (random colour vectors), 8 test functionals:')
    for a, b, c in zip(EA, EB, E11):
        print('    <{Q0,Q^I_A}O> = %s, <{Q0,Q^I_B}O> = %s, <Q_box^2 O> = %s' % (a, b, c))
    print('    rank [A|B] = %d, rank augmented = %d, solution (kappa_A, kappa_B) = %s  [%.0fs]'
          % (M.rank(), M.row_join(rhs).rank(), sp.linsolve((M, rhs), kA, kB), time.time() - t0))
