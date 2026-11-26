"""d = 4 periods (q = 2, unit weights) and the d = 6 conformal star."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import seq, gf, mzv, sv
from rat import R
from cores import zigzag
def wheel(n): return [(0, i) for i in range(1, n + 1)] + [(i, i % n + 1) for i in range(1, n + 1)]
for name, E, ref in [('W3=K4', wheel(3), '6*z3'), ('W4', wheel(4), '20*z5'), ('W5', wheel(5), '70*z7'),
                     ('W6', wheel(6), '252*z9'), ('Z5', list(zigzag(5)), '441/8*z7')]:
    t = time.time()
    val, pl = seq.period(E, [1] * len(E), 2)
    print('%-6s q=2: %s   expected %s   %s (%.1fs)' % (name, mzv.tostr(val), ref, mzv.tostr(val) == ref, time.time() - t), flush=True)
lam = 2
f = gf.VD.rational(lam, gf.edge_factor(2, 2))
b = gf.bounds_for([('u1', '0', 2), ('u1', '1', 2), ('u1', 'z', 2)], lam)
f = gf.append_particular(f).symmetrize()
f, free = gf.fix_kernel(f, b)
print('d=6 star (2,2,2) == 1/(z zb (z-1)(zb-1)):', (f.explicit() - sv.SVF.rat(R({(0, 0): 1}, (1, 1, 1, 1, 0)))).iszero())
