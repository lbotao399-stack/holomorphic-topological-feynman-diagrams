"""Symbolic det[t; dt/dr] of the S1 chart (own t_of_r), via log-derivatives and a DomainMatrix
determinant over QQ(u, v): det[t; dt/dr] = (prod_e t_e) det[1; d log t_e / d r_k]."""
import sys, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from vlib import t_of_r

for n in [int(x) for x in sys.argv[1:]]:
    t0 = time.time()
    us = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    vs = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    te = [sp.factor(x) for x in t_of_r(n, us, vs)]
    order = []
    for k in range(1, n + 1):
        order += [us[k], vs[k]]
    rows = [[sp.Integer(1)] * len(te)] + [[sp.cancel(sp.diff(sp.log(x), s)) for x in te] for s in order]
    gens = order
    K = sp.QQ.frac_field(*gens)
    dm = DomainMatrix.from_list_sympy(len(rows), len(te), rows).convert_to(K)
    d = dm.det()
    d = K.to_sympy(d)
    pt = sp.Integer(1)
    for x in te:
        pt *= x
    J = sp.factor(sp.cancel(d * pt))
    print("n=%d: det[t; dt/du_1; dt/dv_1; ...; dt/du_n; dt/dv_n] = %s   (%.1f s)" % (n, J, time.time() - t0))
