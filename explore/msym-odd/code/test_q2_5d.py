"""{Q_0, k_A Qhat_A + k_B Qhat_B} = 0 for the HT twist of 5d MSYM on C^{2|2} x R (adapted from
6d-holomorphic-twist-note/code/test_q2.py).  Q_1 = 0 (rho_1 = 0), so the order-hbar^2 relation is
{Q_0, Q_2} = 0; the null space of the linear system gives k_A : k_B.
Usage: python3 test_q2_5d.py so3|so5 [test names]"""
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

from salg import structure_constants, add, mul
from ops5 import Theory, L, KernelData, TwoLoop

name = sys.argv[1] if len(sys.argv) > 1 else 'so3'
n, f = structure_constants(name)
T = Theory(n, f)
t0 = time.time()
K = {k: KernelData(k) for k in ('A', 'B')}
print('kernels: s_A = %+d, s_B = %+d, D_A = %s, D_B = %s  (%.1fs)' % (K['A'].s, K['B'].s, sp.factor(K['A'].D), sp.factor(K['B'].D), time.time() - t0))
print('weights at z = 0 (K/pi):  A:', K['A'].weights(((0, 0, 0),) * 3), '  B:', K['B'].weights(((0, 0, 0),) * 3))
T2 = TwoLoop(T, K)

# Q_0^2 = 0
tests0 = [L('c', 0), L('g1', 1), L('b', 2, (1, 0, 0)), mul(L('g1', 0), L('g2', 1)), mul(L('b', 0), L('c', 1, (0, 1, 0))),
          mul(mul(L('g1', 0), L('g2', 1)), L('b', 2))]
t1 = time.time()
print('Q0^2 = 0 on %d operators:' % len(tests0), all(T.Q0(T.Q0(p)) == {} for p in tests0), '(%.1fs)' % (time.time() - t1), flush=True)


def anti(kind, p):
    return add(T.Q0(T2.Q2(kind, p)), T2.Q2(kind, T.Q0(p)))


tests = {
    'g1_0 g2_1 b_2': mul(mul(L('g1', 0), L('g2', 1)), L('b', 2)),
    'b_0 b_1 b_2': mul(mul(L('b', 0), L('b', 1)), L('b', 2)),
    'g1_0 b_1 b_2': mul(mul(L('g1', 0), L('b', 1)), L('b', 2)),
    'g1_0 g2_0 b_1': mul(mul(L('g1', 0), L('g2', 0)), L('b', 1)),
    'c_0 b_1 b_2': mul(mul(L('c', 0), L('b', 1)), L('b', 2)),
    'g1_0 g2_1': mul(L('g1', 0), L('g2', 1)),
    'b_0 b_1': mul(L('b', 0), L('b', 1)),
    'g1_0 g1_1 g2_2': mul(mul(L('g1', 0), L('g1', 1)), L('g2', 2)),
    'dg1_0 g2_1 b_2': mul(mul(L('g1', 0, (1, 0, 0)), L('g2', 1)), L('b', 2)),
    'g1_0 g2_1 db_2': mul(mul(L('g1', 0), L('g2', 1)), L('b', 2, (0, 1, 0))),
}
if len(sys.argv) > 2:
    tests = {k: v for k, v in tests.items() if k in sys.argv[2:]}
rows = []
for nm, p in tests.items():
    t1 = time.time()
    qa, qb = T2.Q2('A', p), T2.Q2('B', p)
    a = anti('A', p)
    b = anti('B', p)
    eqs = [(sp.Rational(a.get(m, 0)), sp.Rational(b.get(m, 0))) for m in set(a) | set(b)]
    if eqs:
        Mx = sp.Matrix(eqs)
        rk = Mx.rank()
        ns = Mx.nullspace()
        nsv = [list(v) for v in ns]
    else:
        rk, nsv = 0, 'all'
    print('%s: |QA O| = %d, |QB O| = %d, |{Q0,QA} O| = %d, |{Q0,QB} O| = %d, equations %d, rank %d, null space %s  (%.1fs)'
          % (nm, len(qa), len(qb), len(a), len(b), len(eqs), rk, nsv, time.time() - t1), flush=True)
    rows += eqs
if rows:
    Mx = sp.Matrix(rows)
    print('all tests: equations %d, rank %d, null space (k_A, k_B): %s' % (len(rows), Mx.rank(), [list(v) for v in Mx.nullspace()]))
print('total %.1fs' % (time.time() - t0))
