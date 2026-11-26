"""Control: {hat Q_0, hat Q_1} with a wrong pentagon kernel (moments of the cube (0,1)^4 instead of the ordered
simplex 1>r1>..>r4>0, i.e. wrong displacement weights) must fail on operators with input derivatives."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
from fractions import Fraction as Fr
import ops_c41 as O
O.pent_moment = lambda a, b, c, d: Fr(1, (a + 1) * (b + 1) * (c + 1) * (d + 1))
n, f = O.structure_constants('su2')
T = O.Theory(n, f)
e = lambda a: tuple(1 if i == a else 0 for i in range(4))
for nm, p in [('g0 g1', O.mul(O.L('g', 0), O.L('g', 1))), ('d1g0 g1', O.mul(O.L('g', 0, e(0)), O.L('g', 1))),
              ('d1g0 d2g1', O.mul(O.L('g', 0, e(0)), O.L('g', 1, e(1)))), ('d1d2g0 g1', O.mul(O.L('g', 0, (1, 1, 0, 0)), O.L('g', 1)))]:
    s = O.add(T.Q0(T.Q1(p)), T.Q1(T.Q0(p)))
    print('control (cube moments): {Q0,Q1} on %-10s -> %s' % (nm, 'ZERO' if s == {} else 'NONZERO (%d terms)' % len(s)))
