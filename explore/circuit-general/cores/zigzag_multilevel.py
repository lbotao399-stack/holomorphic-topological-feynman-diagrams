"""Z5: joint relation system over levels q = 3..qmax, linked by the dimension relation
I(a; q-1) = sum_C I(a + 1_C; q); then P(Z5; 1^10; 2) = sum_C I(1 + 1_C; 3) vs 441 zeta(7)/8."""
import sys
import time
from fractions import Fraction as F

import cores
import gutil
import msys

Z5 = cores.zigzag(5)
qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
autos = gutil.edge_automorphisms(Z5)
cot = gutil.cotrees(Z5)
systems = {}
offset = {}
rows = []
total = 0
t = time.time()
for q in range(3, qmax + 1):
    S = msys.MultiSystem([Z5], q)
    S.lower = cores.lower_core_value
    S.build()
    systems[q] = S
    offset[q] = total
    for row, const in S.rows:
        rows.append(({offset[q] + j: v for j, v in row.items()}, const))
    total += len(S.reps)
# R1 links
nlink = 0
for q in range(4, qmax + 1):
    lo, hi = systems[q - 1], systems[q]
    for j, (f, rep) in enumerate(lo.reps):
        row = {offset[q - 1] + j: F(1)}
        ok = True
        for C in cot:
            b = list(rep)
            for e in C:
                b[e] += 1
            key = (0, gutil.orbit_rep(tuple(b), autos))
            if key not in hi.index:
                ok = False
                break
            jj = offset[q] + hi.index[key]
            row[jj] = row.get(jj, 0) - 1
        if ok:
            rows.append(({k: v for k, v in row.items() if v}, {}))
            nlink += 1
print('joint: %d unknowns, %d rows (%d R1 links)  %.1fs' % (total, len(rows), nlink, time.time() - t), flush=True)


import modsolve
# complexity order: pivot on the highest level, then largest max power, then lexicographic
level_of = {}
for q, S in systems.items():
    for j, (f, rep) in enumerate(S.reps):
        level_of[offset[q] + j] = (q, max(rep), sum(x * x for x in rep), rep)
sol, rank = modsolve.solve(rows, total, rank_key=lambda c: level_of[c])
solved = sol
print('rank %d of %d; determined per level:' % (rank, total),
      {q: sum(1 for j in range(len(S.reps)) if offset[q] + j in solved) for q, S in systems.items()},
      'of', {q: len(S.reps) for q, S in systems.items()}, flush=True)
S3 = systems[3]
tot = {}
okall = True
for C in cot:
    b = [1] * len(Z5)
    for e in C:
        b[e] += 1
    key = (0, gutil.orbit_rep(tuple(b), autos))
    j = offset[3] + S3.index[key]
    if j not in solved:
        okall = False
        break
    for k, x in solved[j].items():
        tot[k] = tot.get(k, 0) + x
if okall:
    tot = {k: x for k, x in tot.items() if x}
    print('P(Z5) at q=2 via R1:', tot, '| 441/8 zeta(7):', tot == {'z7': F(441, 8)}, flush=True)
else:
    print('level-3 vectors needed for P(Z5) not all determined', flush=True)
