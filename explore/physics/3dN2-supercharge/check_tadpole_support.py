"""Tadpole part of Qh2 on all undifferentiated quadratic monomials: support only on B_a B_b, value (2/3) d_abc d^2 c^c.
Log: check_tadpole_support.log"""
import itertools
from ht3d import Engine, mono, pstr, load_kernels
from theories import *
K = load_kernels(4)
for th in [gl_theory(2, 1, 0), gl_theory(2, 2, 1), u1([1, 2, -3]), u1([1, 1])]:
    E = Engine(th, K); A = E.A
    pool = [('c', a, 0) for a in range(th.G)] + [('B', a, 0) for a in range(th.G)]
    for k, m in enumerate(th.multiplets):
        pool += [('p%d' % k, i, 0) for i in range(m[0])] + [('s%d' % k, i, 0) for i in range(m[0])]
    nonBB = []
    nBB = 0
    for combo in itertools.combinations_with_replacement(pool, 2):
        O = mono(A, *combo)
        if not O:
            continue
        t = E.Q2(O, 'tad')
        if t:
            if all(L[0] == 'B' for L in combo):
                nBB += 1
            else:
                nonBB.append((O, t))
    print("%s: quadratic monomials with nonzero tadpole part: B B: %d ; other: %d %s" % (th.name, nBB, len(nonBB), [pstr(o) for o, _ in nonBB[:3]]))
