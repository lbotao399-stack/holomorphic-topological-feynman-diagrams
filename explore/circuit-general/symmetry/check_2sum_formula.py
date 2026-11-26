"""Decompleted product identity (2-vertex cut {a,b}, no edge ab):
  I_q(H; a) = Gamma(q) I_q(H1 + e_{w1}) I_q(H2 + e_{w2}) / (Gamma(q - w1) Gamma(q - w2)),
  w_i = q (n_i + 1) - sum_{H_i} a  (alpha-weight of the closing edge ab),  w1 + w2 = q.
Check on K4+K4 (H_i = K4 - ab) at q = 3..6 against msys-independent engine values is circular, so
the right side is evaluated by cores.period (K4 = W3) and the left side by msys where determined
(q = 3, 4) and by cores-free engine routes otherwise (product split = completed form)."""
import sys
from fractions import Fraction as F
from math import factorial

import engine as en
import cores
import gutil
import msys

K4K4 = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5))
S1 = [0, 1, 4, 5, 6]          # side {1,2}: 01, 02, 12, 15, 25
S2 = [2, 3, 7, 8, 9]          # side {3,4}: 03, 04, 34, 35, 45
for q in [3, 4, 5, 6]:
    autos = gutil.edge_automorphisms(K4K4)
    reps = set()
    for c in gutil.compositions(5 * q, 10, 1, q - 1):
        if gutil.convergent(K4K4, list(c), q):
            reps.add(gutil.orbit_rep(c, autos))
    ref = {}
    if q <= 4:
        S = msys.MultiSystem([K4K4], q)
        S.lower = cores.lower_core_value
        S.build(verbose=False)
        solved, piv = S.solve(verbose=False)
        ref = {S.reps[j][1]: en.from_cores(v) for j, v in solved.items()}
    ok = bad = okm = 0
    for al in sorted(reps):
        val = {(): F(1)}
        fac = F(factorial(q - 1))
        for side in (S1, S2):
            w = 3 * q - sum(al[i] for i in side)          # n_i = 2
            edges = [K4K4[i] for i in side] + [(0, 5)]
            pw = [al[i] for i in side] + [w]
            I = en.from_cores(cores.period(edges, pw, q))
            val = en.vmul(val, I)
            fac /= factorial(q - w - 1)
        rhs = en.vscale(val, fac)
        lhs = en.period(K4K4, al, q, use_m2=False)[0]
        if lhs == rhs:
            ok += 1
        else:
            bad += 1
            print('MISMATCH', q, al, flush=True)
        if al in ref:
            okm += ref[al] == rhs
    print('K4+K4 q=%d orbits %d: formula = engine %d, mismatches %d; formula = msys-determined %d/%d'
          % (q, len(reps), ok, bad, okm, len(ref)), flush=True)
