"""Are fan and tadpole parts of Qh2 separately Q0-closed?  {Q0, Qh2^fan}, {Q0, Qh2^tad} on monomials;
negative control: perturbed kernel (fan K_1 moments scaled by 2).  Log: check_split.log"""
import itertools
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels
from theories import *
K = load_kernels(4)
# perturbed tables: multiply all moments of total insertion order 1 by 2 (fan-type and tadpole-type)
Kp = {}
for edges, (mom, M) in K.items():
    Kp[edges] = ({p: ({k: v * (2 if sum(p) == 1 else 1) for k, v in d.items()}) for p, d in mom.items()}, M)
for th in [u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}), gl_theory(2, 1, 0), su2_adj(1)]:
    E = Engine(th, K); Ep = Engine(th, Kp); A = E.A
    pool = [('c', a, 0) for a in range(th.G)] + [('B', a, d) for a in range(th.G) for d in (0, 1)]
    for k, m in enumerate(th.multiplets):
        pool += [('p%d' % k, i, d) for i in range(m[0]) for d in (0, 1)] + [('s%d' % k, i, 0) for i in range(m[0])]
    nf = nt = n = npert = 0
    for deg in (2, 3, 4):
        for combo in itertools.combinations_with_replacement(pool, deg):
            if sum(L[2] for L in combo) > (0 if th.G == 4 else 1) or (deg == 4 and th.G > 1):
                continue
            O = mono(A, *combo)
            if not O:
                continue
            n += 1
            cf = psum(E.Q0(E.Q2(O, 'fan')), E.Q2(E.Q0(O), 'fan'))
            ct = psum(E.Q0(E.Q2(O, 'tad')), E.Q2(E.Q0(O), 'tad'))
            cp = psum(E.Q0(Ep.Q2(O)), Ep.Q2(E.Q0(O)))
            nf += bool(cf); nt += bool(ct); npert += bool(cp)
    print("%s: %d monomials (D <= 1): {Q0,fan} != 0 on %d, {Q0,tad} != 0 on %d ; perturbed K_1 (x2): {Q0,Q2'} != 0 on %d"
          % (th.name, n, nf, nt, npert))
