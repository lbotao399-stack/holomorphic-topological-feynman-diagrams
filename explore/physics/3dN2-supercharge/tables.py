"""Explicit Qh2 (= Q_2 / N_2) on undifferentiated quadratic and cubic monomials; log tables.log"""
import itertools
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, load_kernels
from theories import *
K = load_kernels(4)
for th in [u1([1]), su2_adj(1)]:
    E = Engine(th, K); A = E.A
    pool = [('c', a, 0) for a in range(th.G)] + [('B', a, 0) for a in range(th.G)]
    for k, m in enumerate(th.multiplets):
        pool += [('p%d' % k, i, 0) for i in range(m[0])] + [('s%d' % k, i, 0) for i in range(m[0])]
    print("==", th.name, "(I = (1/2) f B c c + psi c rho phi, eta(c,B) = eta(phi,psi) = 1, kappa = W = 0)")
    cnt = 0
    for deg in (2, 3):
        for combo in itertools.combinations_with_replacement(pool, deg):
            O = mono(A, *combo)
            if not O:
                continue
            f, t = E.Q2(O, 'fan'), E.Q2(O, 'tad')
            if f or t:
                cnt += 1
                if th.G == 1 or cnt <= 12:
                    print("  Qh2 %s = fan: %s ; tad: %s" % (pstr(O), pstr(f), pstr(t)))
    print("  (%d nonzero)" % cnt)
th = u1([1])
E = Engine(th, K); A = E.A
print("== U(1), one chiral q=1: first-derivative examples")
for lets in [[('B', 0, 0), ('B', 0, 0), ('B', 0, 1)], [('B', 0, 1), ('B', 0, 1), ('B', 0, 0)], [('B', 0, 0), ('B', 0, 0), ('s0', 0, 1)], [('B', 0, 0), ('B', 0, 1)]]:
    O = mono(A, *lets)
    print("  Qh2 %s = %s" % (pstr(O), pstr(E.Q2(O))))
