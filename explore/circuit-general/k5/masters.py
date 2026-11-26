"""Master sets of the tower: level 2 (one class) + level 3 (15 classes) with (X) at level 3 and (DS) 2 -> 3.
Lists the triples of level-3 classes that, together with J_2, determine every level-3 class; prints values."""
import itertools
import os

import k6
import linsys
from solve_level import classes, load_values, read_rows, SCR, class_file
import subprocess

q = 3
cls = classes(3)
low = classes(2)
rf = os.path.join(SCR, 'rows3.txt')
with open(rf, 'w') as f:
    subprocess.check_call([os.path.join(SCR, 'xrows'), '3', class_file(3), class_file(2)], stdout=f)
rows = read_rows(rf)
v3 = load_values(3)
v2 = load_values(2)
n = len(cls)
print('level 2:', low[0], v2[low[0]])
for j, b in enumerate(cls):
    print('level 3 class %2d %s %s decompletable at %s: J = %s' % (j, b, k6.kind(b, 3), k6.k5_decompletions(b, 3), v3[b]))
good = []
for tri in itertools.combinations(range(n), 3):
    unk = [j for j in range(n) if j not in tri]
    cidx = {j: i for i, j in enumerate(unk)}
    mrows = []
    for r in rows:
        m = {cidx[k]: linsys.to_mod(v) for k, v in r.items() if not isinstance(k, tuple) and k in cidx}
        m = {k: v for k, v in m.items() if v}
        if m:
            mrows.append(m)
    piv = linsys.rref_mod(mrows, len(unk))
    if len(piv) == len(unk):
        good.append(tri)
print('master triples (with J_2): %d of %d' % (len(good), len(list(itertools.combinations(range(n), 3)))))
print('examples:', good[:12])
# rank of X alone and X + DS
for tag, use in (('X only', False), ('X + DS(2->3)', True)):
    mrows = []
    for r in rows:
        if any(isinstance(k, tuple) for k in r) and not use:
            continue
        m = {k: linsys.to_mod(v) for k, v in r.items() if not isinstance(k, tuple)}
        m = {k: v for k, v in m.items() if v}
        if m:
            mrows.append(m)
    print('%s: rank %d of %d' % (tag, len(linsys.rref_mod(mrows, n)), n))
