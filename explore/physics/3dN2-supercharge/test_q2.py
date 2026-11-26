"""{Q_0, Q_2} = 0 on all monomials of bounded degree/derivative order; separate fan/tadpole parts;
nontriviality (Q0Q2, Q2Q0 individually nonzero).  Usage: python3 test_q2.py <theory> <maxdeg> <maxd>
Log: test_q2_<theory>.log"""
import sys, time, itertools
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels
from theories import *

name = sys.argv[1]; maxdeg = int(sys.argv[2]); maxd = int(sys.argv[3])
TH = {
    "su2adj": lambda: su2_adj(1),
    "su2adj_mass_cs": lambda: su2_adj(1, W=W_mass_adj(1), kappa=Fr(1, 3)),
    "su2adj2": lambda: su2_adj(2, W={((0, a), (1, a)): 1 for a in range(3)}, kappa=Fr(2, 5)),
    "su2N8": lambda: su2_adj(3, W=W_N8()),
    "u1": lambda: u1([1, -1], W={((0, 0), (1, 0)): 1}, kappa=2),
    "u1q": lambda: u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}, kappa=Fr(1, 2)),
    "su2pure": lambda: su2_adj(0, kappa=Fr(3, 7)),
}
th = TH[name]()
K = load_kernels(4)
E = Engine(th, K)
A = E.A
pool = [('c', a, d) for a in range(th.G) for d in range(maxd + 1)] + [('B', a, d) for a in range(th.G) for d in range(maxd + 1)]
for k, mlt in enumerate(th.multiplets):
    pool += [('p%d' % k, i, d) for i in range(mlt[0]) for d in range(maxd + 1)]
    pool += [('s%d' % k, i, d) for i in range(mlt[0]) for d in range(maxd + 1)]
print("theory:", th.name, "W:", th.W, "kappa:", th.kappa, "pool size", len(pool))
t0 = time.time()
n_mono = n_nontriv = n_bad = n_tad = n_q2 = 0
for deg in range(1, maxdeg + 1):
    for combo in itertools.combinations_with_replacement(pool, deg):
        if sum(L[2] for L in combo) > maxd:
            continue
        O = mono(A, *combo)
        if not O:
            continue
        n_mono += 1
        q2 = E.Q2(O)
        tad = E.Q2(O, 'tad')
        a = E.Q0(q2)
        b = E.Q2(E.Q0(O))
        C = psum(a, b)
        if q2:
            n_q2 += 1
        if tad:
            n_tad += 1
            if n_tad <= 3:
                print("  tadpole part nonzero on", pstr(O), ":", pstr(tad))
        if a and b:
            n_nontriv += 1
        if C:
            n_bad += 1
            if n_bad <= 5:
                print("  FAIL on", pstr(O), ": {Q0,Q2} =", pstr(C))
print("monomials %d (deg<=%d, total derivatives<=%d): Q2 != 0 on %d, tadpole part != 0 on %d, "
      "Q0Q2 and Q2Q0 both != 0 on %d, {Q0,Q2} != 0 on %d   [%.0fs]"
      % (n_mono, maxdeg, maxd, n_q2, n_tad, n_nontriv, n_bad, time.time() - t0))
