"""Q_0 = {I,.} from the one-edge graph: Q_0^2 = 0 and comparison with CDG (3.34).  Log: test_q0.log"""
import sys, itertools, random
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, psum, pscale
from theories import *

def check_sq(E, letters_pool, nmax=3, ntest=200, seed=1):
    A = E.A
    rng = random.Random(seed)
    bad = 0
    tested = 0
    for k in range(1, nmax + 1):
        for combo in itertools.combinations(letters_pool, k):
            P = mono(A, *combo)
            if not P:
                continue
            R = E.Q0(E.Q0(P))
            tested += 1
            if R:
                bad += 1
                if bad < 4:
                    print("   Q0^2 !=0 on", pstr(P), "->", pstr(R))
    return tested, bad

for sg in (-1, 1):
    for th in [su2_adj(1, s_g=sg), su2_adj(1, W=W_mass_adj(1), kappa=Fr(1, 3), s_g=sg),
               u1([1, -1], W={((0, 0), (1, 0)): 1}, kappa=2, s_g=sg),
               su2_adj(3, W=W_N8(), s_g=sg)]:
        E = Engine(th, kernels={})
        pool = [('c', a, d) for a in range(th.G) for d in (0, 1)] + [('B', a, 0) for a in range(th.G)]
        for k, mlt in enumerate(th.multiplets):
            pool += [('p%d' % k, i, 0) for i in range(mlt[0])] + [('s%d' % k, i, 0) for i in range(mlt[0])]
        nmax = 3 if len(pool) < 20 else 2
        t, b = check_sq(E, pool, nmax=nmax)
        print("s_g=%d %-28s kappa=%s: Q0^2=0 tested %d monomials (deg<=%d), failures %d" % (sg, th.name, th.kappa, t, nmax, b))
    print()

th = su2_adj(1, W=W_mass_adj(1), kappa=Fr(1, 3), s_g=1)
E = Engine(th, kernels={})
A = E.A
print("letter Q_0 (s_g=+1, SU(2)+adj, W=(1/2)phi.phi, kappa=1/3):")
for L in [('c', 0, 0), ('B', 0, 0), ('p0', 0, 0), ('s0', 0, 0), ('c', 0, 1), ('B', 0, 1)]:
    print("  Q0", pstr(mono(A, L)), "=", pstr(E.Q0(mono(A, L))))
