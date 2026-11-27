"""2d: hatQ1(d^a phi^A d^b phi^B) = (a! b!/(a+b+1)!) * [hatQ1(phi^A phi^B) with d -> d^(a+b+1)] for 0 <= a, b <= 3,
all letter pairs, theories (0,2) sl2+fund*4, (2,2) sl2+fund, (4,4) sl2+fund.  Log: c_q1jets.log"""
import math
from fractions import Fraction as Fr
import polyops as po
import theories
from wick import Algebra
from explicit import GraphOp
import rgraphs
for th, gname, matter in (('2d02', 'sl2', 'fund*4'), ('2d22', 'sl2', 'fund'), ('2d44', 'sl2', 'fund')):
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    g = dict(rgraphs.rooted_trees(1, 1)[0]); g['name'] = 'B1'
    Q1 = GraphOp(alg, g)
    ok, tested = True, 0
    for A in range(alg.N):
        for B in range(alg.N):
            s0, m0 = po.canon([(A, (0,)), (B, (0,))], alg.pi)
            if not s0:
                continue
            base = Q1.apply({m0: s0})
            if not base:
                continue
            for a in range(4):
                for b in range(4):
                    s, m = po.canon([(A, (a,)), (B, (b,))], alg.pi)
                    if not s:
                        continue
                    img = Q1.apply({m: s})
                    c = po.modq(Fr(math.factorial(a) * math.factorial(b), math.factorial(a + b + 1)))
                    ref = {((mm[0][0], (a + b + 1,)),): v * c % po.P for mm, v in base.items()}
                    tested += 1
                    if img != ref:
                        ok = False
    print('%s %s %s: jet rule holds on %d (pair, a, b) cases: %s' % (th, gname, matter, tested, ok), flush=True)
