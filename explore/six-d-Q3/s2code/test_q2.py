import sys, time
from fractions import Fraction as Fr
from salg import structure_constants, add, mul, scale
from ops import Theory, L
from ops2 import TwoLoop

name = sys.argv[1] if len(sys.argv) > 1 else 'su2'
n, f = structure_constants(name)
T = Theory(n, f); T2 = TwoLoop(T)

def anti(kind, p):
    return add(T.Q0(T2.Q2(kind, p)), T2.Q2(kind, T.Q0(p)))

tests = {
    'g1_0 g2_1 b_2': mul(mul(L('g1', 0), L('g2', 1)), L('b', 2)),
    'b_0 b_1 b_2':   mul(mul(L('b', 0), L('b', 1)), L('b', 2)),
    'g1_0 b_1 b_2':  mul(mul(L('g1', 0), L('b', 1)), L('b', 2)),
    'g1_0 g2_0 b_1': mul(mul(L('g1', 0), L('g2', 0)), L('b', 1)),
    'c_0 b_1 b_2':   mul(mul(L('c', 0), L('b', 1)), L('b', 2)),
    'g1_0 g2_1':     mul(L('g1', 0), L('g2', 1)),
    'b_0 b_1':       mul(L('b', 0), L('b', 1)),
    'dg1_0 g2_1 b_2': mul(mul(L('g1', 0, (1, 0, 0)), L('g2', 1)), L('b', 2)),
    'g1_0 g2_1 db_2': mul(mul(L('g1', 0), L('g2', 1)), L('b', 2, (0, 1, 0))),
}
if len(sys.argv) > 2:
    tests = {k: v for k, v in tests.items() if k in sys.argv[2:]}
rows = []
for nm, p in tests.items():
    t0 = time.time()
    lhs = T.Q1(T.Q1(p)); a = anti('A', p); b = anti('B', p)
    print(f'{nm}: |Q1^2| = {len(lhs)}, |{{Q0,QA}}| = {len(a)}, |{{Q0,QB}}| = {len(b)}  ({time.time()-t0:.1f}s)', flush=True)
    rows.append((nm, lhs, a, b))

# solve  Q1hat^2 = kA {Q0hat,QAhat} + kB {Q0hat,QBhat}
import sympy as sp
kA, kB = sp.symbols('kA kB')
eqs = []
for nm, lhs, a, b in rows:
    for m in set(lhs) | set(a) | set(b):
        eqs.append((sp.Rational(lhs.get(m, 0)), sp.Rational(a.get(m, 0)), sp.Rational(b.get(m, 0))))
M = sp.Matrix([[e[1], e[2]] for e in eqs]); rhs = sp.Matrix([e[0] for e in eqs])
print('equations:', len(eqs), 'rank of [A|B]:', M.rank(), 'rank augmented:', M.row_join(rhs).rank())
sol = sp.linsolve((M, rhs), kA, kB)
print('solution (kA, kB):', sol)
