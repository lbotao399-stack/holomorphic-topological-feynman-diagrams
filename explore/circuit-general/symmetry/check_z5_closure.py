"""Z5 at level q: msys relation system (triangle rules, IBP, boundaries, star-triangle) closed with
symmetry data: (i) engine values (boundary rows), (ii) symmetry equalities between Z5 orbits
I(Z5; a) = f I(Z5; a') from completion / Fourier / twist / uniqueness orbits (decompletions isomorphic to Z5)."""
import sys
import time
from fractions import Fraction as F
from math import factorial

import engine as en
import cores
import gutil
import msys

q = int(sys.argv[1]) if len(sys.argv) > 1 else 4
Z5 = cores.zigzag(5)
t = time.time()
S = msys.MultiSystem([Z5], q)
S.lower = cores.lower_core_value
S.build(verbose=True)
base_rows = len(S.rows)
solved0, piv0 = S.solve(verbose=True)
nrep = len(S.reps)
# (i) engine values
vals = {}
for j, (f, rep) in enumerate(list(S.reps)):
    r = en.period(Z5, rep, q)
    if r is not None:
        vals[j] = r[0]
for j, v in vals.items():
    S.rows.append(({j: F(1)}, {('1' if m == () else m[0]): x for m, x in v.items()}))
print('engine-evaluated orbits %d of %d' % (len(vals), nrep), flush=True)
solved1, piv1 = S.solve(verbose=True)
# (ii) symmetry equalities
fam_key = msys.struct_key(Z5)
neq = 0
for j, (f, rep) in enumerate(list(S.reps)[:nrep]):
    if j in vals:
        continue
    V, E, inf = en.complete(Z5, rep, q)
    root = F(1) / factorial(q - 1)
    for a in rep:
        root *= factorial(q - a - 1)           # I(root) = prod Gamma(b)/Gamma(q) Pi
    for V2, E2, fac in en.orbit_states(V, E, q, budget=60):
        for c in range(V2):
            H = en.decomplete(V2, E2, c, q)
            if H is None or len(H[0]) != 10:
                continue
            edges, al, be = H
            if msys.struct_key(edges) != fam_key:
                continue
            mapped = S.map_to_family(edges, al, 0)
            if mapped is None or not gutil.convergent(tuple(Z5), mapped, q):
                continue
            j2 = S.var(0, mapped)
            if j2 == j:
                continue
            # Pi(state) = Gamma(q) I(H)/prod Gamma(b_H)
            g = F(factorial(q - 1))
            for b in be:
                g /= factorial(b - 1)
            coef = root * fac * g             # I(root) = coef * I(Z5; mapped)
            S.rows.append(({j: F(1), j2: -coef} if j2 != j else {}, {}))
            neq += 1
print('symmetry equalities %d' % neq, flush=True)
solved2, piv2 = S.solve(verbose=True)
print('Z5 q=%d: orbits %d | msys alone: determined %d | + engine values (%d): determined %d | + symmetry equalities: determined %d  (%.1fs)'
      % (q, nrep, len(solved0), len(vals), len([j for j in solved1 if j < nrep]),
         len([j for j in solved2 if j < nrep]), time.time() - t), flush=True)
