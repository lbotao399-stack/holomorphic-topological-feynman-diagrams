"""Insert the exact value of the isolated odd-level class (check_aux_odd.py) into the level-q table
(replacing its free label) and recheck all level-q rows.  usage: python3 set_aux.py Q NUM DEN"""
import sys
from fractions import Fraction as F
from solve_level import classes, load_values, save_values
import fix_rows_lib
q, num, den = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
k = (q - 3) // 2
lam = q - 1
b = (-k, -k, lam, lam, lam, -k, lam, lam, lam, lam, lam, lam, -k, -k, -k)
cls = classes(q)
v = load_values(q)
old = v[b]
lab = [l for l in old if l.startswith('u')]
assert len(old) == 1 and lab, old
val = F(num, den)
new = {}
for j, c in enumerate(cls):
    nv = {}
    for l, x in v[c].items():
        if l == lab[0]:
            nv['1'] = nv.get('1', 0) + x * val
        else:
            nv[l] = nv.get(l, 0) + x
    new[j] = {l: x for l, x in nv.items() if x}
save_values(q, cls, new)
print('class %s: %s -> %s; rows with nonzero residual after insertion: %d' % (
    b, old, new[cls.index(b)], fix_rows_lib.count_bad(q)))
