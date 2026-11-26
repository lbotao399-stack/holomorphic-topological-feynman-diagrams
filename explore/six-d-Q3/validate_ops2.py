"""V5: nonzero compositions in so(5), slot evaluation against S2's code:
   hat Q_0 hat Q_box (g1_0 g2_1), hat Q_box hat Q_0 (g1_0 g2_1)                 (4 final letters)
   hat Q_0 hat Q_A (g1_0 g2_1 b_2), hat Q_box hat Q_box (g1_0 g2_1 b_2)         (5 final letters)
   hat Q_A hat Q_0 (g1_0 g2_1 b_2)                                             (5 final letters)"""
import time
import numpy as np
import opcheck as oc
from salg import mul
from ops import Theory, L
from ops2 import TwoLoop
from validate_ops import pair, unit

n, f, fd = oc.fdense('so5')
T = Theory(n, fd)
T2 = TwoLoop(T)
p2 = mul(L('g1', 0), L('g2', 1))
p3 = mul(p2, L('b', 2))
O2 = [[unit(n, 0), unit(n, 1)]]
O3 = [[unit(n, 0), unit(n, 1), unit(n, 2)]]
cases = [('Q0*box', lambda: T.Q0(T.Q1(p2)), ['Q0', 'box'], [3, 4], O2, 4),
         ('box*Q0', lambda: T.Q1(T.Q0(p2)), ['box', 'Q0'], [3, 4], O2, 4),
         ('Q0*A', lambda: T.Q0(T2.Q2('A', p3)), ['Q0', 'A'], [3, 4, 2], O3, 5),
         ('box*box', lambda: T.Q1(T.Q1(p3)), ['box', 'box'], [3, 4, 2], O3, 5),
         ('A*Q0', lambda: T2.Q2('A', T.Q0(p3)), ['A', 'Q0'], [3, 4, 2], O3, 5)]
for name, fn, seq, types, Ov, nf in cases:
    t0 = time.time()
    P = fn()
    ts = time.time() - t0
    for seed in range(2):
        td = oc.TestData(n, 40 + seed, n=nf)
        s2 = pair(P, td)
        t1 = time.time()
        sk, nh = oc.evaluate(seq, types, f, Ov, [td])
        print('V5 so5 %-8s seed %d: S2 code %s, slot evaluation %s, equal %s  [S2 code %.0fs, %d monomials; slot %.1fs]'
              % (name, seed, s2, sk[0], s2 == sk[0], ts, len(P), time.time() - t1), flush=True)
