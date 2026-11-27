"""[fork A] Explicit Q_1 = kappa_1 hatQ_1, Q_2 = kappa_2 hatQ_2 (kappa_n = 1/2, c = 1, Q_0 = hatQ_0 as in n4_llc_check)
on gauge-invariant c_0-free bilinears (Q_1) and trilinears (Q_2) without derivatives, 4d N=2 sl2 + one fundamental hyper.
Letters (theories.build('4dN2','sl2','fund')): X = c (ghost, odd), u = hyper phi (even, U = R + R^*), Xv = xi (Koszul
generator, odd), X* = b (even), u* = varsigma (odd), Xv* = a (adjoint scalar, even); sl2 basis (e, f, h) = (0, 1, 2).
Fan kernels (no input jets): I_n = (-1)^{n(n-1)/2} det(lam_1, lam_2)^n / ((n!)^3 (n+1)) [chk validate_kernels.log (a)].
usage: python3 a_qn_examples.py [matter]   log: a_qn_examples.log"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np

import polyops as po
import theories
import lifts as lf
from wick import Algebra
from explicit import graph_ops

matter = sys.argv[1] if len(sys.argv) > 1 else 'fund'
dq2 = tuple(int(x) for x in sys.argv[2].split('.')) if len(sys.argv) > 2 else (0, 0)
only2 = len(sys.argv) > 3
H, sl = theories.build('4dN2', 'sl2', matter)
alg = Algebra(sl, H)
g = theories.gauge('sl2')
gidx = sl.blocks['g']
th = lf.LiftTheory(alg, cartan=[gidx[2]], posroots=[gidx[0]])
G1 = graph_ops(alg, 1)[0]
G2 = graph_ops(alg, 2)[0]
k1 = k2 = po.modq(Fr(1, 2))
names = sl.names


def show(p, m=14):
    return ' + '.join(po.show(p, names, H, maxterms=m).split('\n'))


def block_content(mono):
    return tuple(sorted(th.blk_of[A] for A, _ in mono))


print('4d N=2, gauge sl2, matter %s: blocks %s' % (matter, {b: len(v) for b, v in th.part.items()}))
print('charges per block:', th.charges)
for k, (Gop, kap, lab) in (((3, (G2, k2, 'Q_2')),) if only2 else ((2, (G1, k1, 'Q_1')), (3, (G2, k2, 'Q_2')))):
    # all charge sectors with k letters, no derivatives
    qs = {}
    for m in lf.sector_monomials(th, None, dq2 if k == 3 else (0, 0), k):
        q = tuple(sum(th.wt[A][i] for A, _ in m) for i in range(len(th.charges)))
        qs.setdefault(q, []).append(m)
    for q in sorted(qs):
        monos = qs[q]
        inv = lf.invariants(th, monos)
        for v in inv:
            O = lf.vec_to_poly(v, monos)
            img = lf.drop_c0(th, po.pscale(Gop.apply(O), kap))
            cl = lf.drop_c0(th, alg.q0(O))
            contents = sorted({block_content(m) for m in O})
            if only2 and not img:
                continue
            print('%s  q=%s  O ~ %s  [Q0 O %s 0]' % (lab, q, contents, '=' if not cl else '!='))
            print('    O     = %s' % show(O))
            print('    %s O = %s' % (lab, show(img) if img else '0'))
