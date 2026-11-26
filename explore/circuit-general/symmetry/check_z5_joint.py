"""Z5, joint levels q = 3, 4: msys relation systems at both levels + engine values + symmetry equalities
+ dimension links I_3(a) = sum_C I_4(a + 1_C) (rows used only if every term is convergent).
Afterwards: P(Z5) at q = 2, unit weights, = sum_C I_3(1 + 1_C), compared with 441/8 zeta(7)
[src arXiv:1208.1890]; not used as input."""
import time
from fractions import Fraction as F
from math import factorial

import engine as en
import cores
import gutil
import msys

Z5 = cores.zigzag(5)
FK = msys.struct_key(Z5)
t0 = time.time()


def sym_rows(S, f, q, reps, vals):
    rows = []
    for rep in reps:
        j = S.var(f, rep)
        if j in vals:
            continue
        V, E, inf = en.complete(Z5, rep, q)
        root = F(1) / factorial(q - 1)
        for a in rep:
            root *= factorial(q - a - 1)
        for V2, E2, fac in en.orbit_states(V, E, q, budget=60):
            for c in range(V2):
                H = en.decomplete(V2, E2, c, q)
                if H is None or len(H[0]) != 10:
                    continue
                edges, al, be = H
                if msys.struct_key(edges) != FK:
                    continue
                mapped = S.map_to_family(edges, al, f)
                if mapped is None or not gutil.convergent(tuple(Z5), mapped, q):
                    continue
                j2 = S.var(f, mapped)
                if j2 == j:
                    continue
                g = F(factorial(q - 1))
                for b in be:
                    g /= factorial(b - 1)
                rows.append(({j: F(1), j2: -root * fac * g}, {}))
    return rows


S3 = msys.MultiSystem([Z5], 3)
S3.lower = cores.lower_core_value
S3.build(verbose=False)
S = msys.MultiSystem([Z5], 4)
S.lower = cores.lower_core_value
S.build(verbose=False)
n4 = len(S.reps)
# second family slot = level 3
d = dict(S.fam[0])
S.fam.append(d)
m3 = {}
for j3, (f, rep) in enumerate(S3.reps):
    m3[j3] = S.var(1, rep)
for row, const in S3.rows:
    S.rows.append(({m3[j]: v for j, v in row.items()}, dict(const)))
reps3 = [rep for f, rep in S3.reps]
reps4 = [rep for f, rep in S.reps[:n4]]
vals = {}
for f, q, reps in ((1, 3, reps3), (0, 4, reps4)):
    for rep in reps:
        r = en.period(Z5, rep, q)
        if r is not None:
            j = S.var(f, rep)
            vals[j] = r[0]
            S.rows.append(({j: F(1)}, {('1' if m == () else m[0]): x for m, x in r[0].items()}))
nsym = 0
for f, q, reps in ((1, 3, reps3), (0, 4, reps4)):
    rows = sym_rows(S, f, q, reps, vals)
    S.rows.extend(rows)
    nsym += len(rows)
cot = gutil.cotrees(tuple(Z5))
nlink = 0
for rep in reps3:
    row = {S.var(1, rep): F(1)}
    ok = True
    for C in cot:
        b = list(rep)
        for e in C:
            b[e] += 1
        if not gutil.convergent(tuple(Z5), b, 4):
            ok = False
            break
        j = S.var(0, gutil.orbit_rep(tuple(b), S.fam[0]['autos']))
        row[j] = row.get(j, 0) - 1
    if ok:
        S.rows.append(({k: v for k, v in row.items() if v}, {}))
        nlink += 1
print('unknowns %d (level 3: %d, level 4: %d), engine values %d, symmetry rows %d, dimension links %d'
      % (len(S.reps), len(reps3), n4, len(vals), nsym, nlink), flush=True)
solved, piv = S.solve(verbose=True)
d3 = sum(1 for rep in reps3 if S.var(1, rep) in solved)
d4 = sum(1 for j in range(n4) if j in solved)
print('determined: level 3 %d/%d, level 4 %d/%d' % (d3, len(reps3), d4, n4), flush=True)
tot = {}
good = True
for C in cot:
    b = [1] * 10
    for e in C:
        b[e] += 1
    if not gutil.convergent(tuple(Z5), b, 3):
        good = False
        break
    j = S.var(1, gutil.orbit_rep(tuple(b), S.fam[1]['autos']))
    if j not in solved:
        good = False
        break
    for k, v in solved[j].items():
        tot[k] = tot.get(k, 0) + v
tot = {k: v for k, v in tot.items() if v}
print('P(Z5) = sum_C I_3(1+1_C) =', tot if good else 'not determined', ' expected {z7: 441/8}:',
      good and tot == {'z7': F(441, 8)}, ' (%.1fs)' % (time.time() - t0), flush=True)
