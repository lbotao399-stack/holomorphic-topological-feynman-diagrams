"""The level-5 class left free by (X)+(DS from 4): its 6 DS rows to level 6 (one per decompletion vertex c)
give (q-1) J_5 = sum_T prod b_e J_6(...); all 6 must agree.  Prints the value; with 'write' substitutes it
into the level-5 table (label u5_0 -> value) and re-checks the level-6 rows."""
import sys
from fractions import Fraction as F
import k6
from solve_level import classes, read_rows, load_values, save_values, SCR
import os
q = 6
low = classes(5)
lv = load_values(5)
cls = classes(6)
v6 = load_values(6)
ub = [b for b in low if any(l.startswith('u') for l in lv[b])]
print('free level-5 classes:', ub)
rows = read_rows(os.path.join(SCR, 'rows6.txt'))
cands = []
for r in rows:
    Ls = [k for k in r if isinstance(k, tuple)]
    if not Ls:
        continue
    b = low[Ls[0][1]]
    if b not in ub:
        continue
    acc = {}
    for k, c in r.items():
        if isinstance(k, tuple):
            cl = c
            continue
        for lab, x in v6[cls[k]].items():
            acc[lab] = acc.get(lab, 0) + c * x
    val = {lab: -x / cl for lab, x in acc.items() if x}
    cands.append(val)
print('DS rows from the free class: %d; values agree: %s' % (len(cands), all(c == cands[0] for c in cands)))
print('J_5', ub[0], '=', cands[0])
if 'write' in sys.argv:
    lv[ub[0]] = cands[0]
    save_values(5, low, {i: lv[b] for i, b in enumerate(low)})
    print('written')
