"""Classes of level q-1 left free by (X)+(DS from q-2) carry labels u<q-1>_n.  Their DS rows to level q,
  (q-1) J_{q-1}(b) = sum_T prod b_e J_q(...),  one per decompletion vertex c,
fix them once level q is solved without them.  All rows of one class must agree.
usage: python3 fix_u.py Q [write]     (Q = the upper level)"""
import os
import sys
from fractions import Fraction as F

from solve_level import classes, read_rows, load_values, save_values, SCR

q = int(sys.argv[1])
low = classes(q - 1)
lv = load_values(q - 1)
cls = classes(q)
vq = load_values(q)
free = [b for b in low if any(l.startswith('u') for l in lv[b])]
print('free level-%d classes: %s' % (q - 1, free))
if any(any(l.startswith('u%d_' % (q - 1)) for l in v) for v in vq.values()):
    print('WARNING: level-%d values depend on level-%d free labels' % (q, q - 1))
lowidx = {b: i for i, b in enumerate(low)}
want = {lowidx[b]: b for b in free}
vals = {}
for line in open(os.path.join(SCR, 'rows%d.txt' % q)):
    if not line.startswith('D'):
        continue
    p = line.split()
    li = int(p[2][1:])
    if li not in want:
        continue
    cl = F(int(p[1]))
    acc = {}
    items = p[3:]
    for i in range(0, len(items), 2):
        c, j = int(items[i]), int(items[i + 1])
        for lab, x in vq[cls[j]].items():
            acc[lab] = acc.get(lab, 0) + c * x
    vals.setdefault(want[li], []).append({lab: -x / cl for lab, x in acc.items() if x})
for b in free:
    vs = vals.get(b, [])
    print('class %s: %d DS rows, agree %s, J_%d = %s' % (b, len(vs), all(v == vs[0] for v in vs) if vs else None,
                                                       q - 1, vs[0] if vs else None))
if 'write' in sys.argv:
    for b in free:
        vs = vals.get(b, [])
        if vs and all(v == vs[0] for v in vs):
            lv[b] = vs[0]
    save_values(q - 1, low, {i: lv[b] for i, b in enumerate(low)})
    print('written')
