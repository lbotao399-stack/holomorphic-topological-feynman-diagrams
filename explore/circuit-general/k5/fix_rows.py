"""Close the remaining free labels u<q>_n of level q with the full row set (X + DS from q-1), exactly.
Each row residual is linear in the free labels: sum_u c_u u + (constants) = 0; the small system is solved
over Fractions (constants carried per MZV label), the solution substituted, the table saved, and all rows
rechecked.  usage: python3 fix_rows.py Q [write]"""
import os
import sys
from fractions import Fraction as F

from solve_level import classes, load_values, save_values, SCR

q = int(sys.argv[1])
cls = classes(q)
low = classes(q - 1)
vq = load_values(q)
lv = load_values(q - 1)
pre = 'u%d_' % q


def residuals(vals_by_idx, lowvals):
    for line in open(os.path.join(SCR, 'rows%d.txt' % q)):
        p = line.split()
        r = {}
        if p[0] == 'X':
            items = p[1:]
        else:
            r[('L', int(p[2][1:]))] = F(int(p[1]))
            items = p[3:]
        for i in range(0, len(items), 2):
            j = int(items[i + 1])
            r[j] = r.get(j, 0) + int(items[i])
        acc = {}
        for kk, c in r.items():
            if not c:
                continue
            vv = lowvals[low[kk[1]]] if isinstance(kk, tuple) else vals_by_idx[kk]
            for lab, x in vv.items():
                acc[lab] = acc.get(lab, 0) + c * x
        acc = {k: v for k, v in acc.items() if v}
        if acc:
            yield acc


vals = {j: vq[b] for j, b in enumerate(cls)}
free = sorted({l for v in vals.values() for l in v if l.startswith(pre)})
print('level %d free labels: %s' % (q, free))
eqs = []
for acc in residuals(vals, lv):
    eqs.append(acc)
print('rows with nonzero residual: %d' % len(eqs))
# eliminate: unknowns = free labels; constants = other labels
piv = {}
for acc in eqs:
    u = {l: x for l, x in acc.items() if l.startswith(pre)}
    cst = {l: x for l, x in acc.items() if not l.startswith(pre)}
    while u:
        l0 = min(u)
        if l0 in piv:
            pu, pc = piv[l0]
            f = u[l0]
            for l, x in pu.items():
                u[l] = u.get(l, 0) - f * x
                if not u[l]:
                    del u[l]
            for l, x in pc.items():
                cst[l] = cst.get(l, 0) - f * x
                if not cst[l]:
                    del cst[l]
        else:
            f = 1 / u[l0]
            piv[l0] = ({l: x * f for l, x in u.items()}, {l: x * f for l, x in cst.items()})
            break
    else:
        if cst:
            print('INCONSISTENT residual', cst)
sol = {}
for l0 in sorted(piv, reverse=True):
    pu, pc = piv[l0]
    val = {l: -x for l, x in pc.items()}
    ok = True
    for l, x in pu.items():
        if l == l0:
            continue
        if l not in sol:
            ok = False
            break
        for l2, y in sol[l].items():
            val[l2] = val.get(l2, 0) - x * y
    if ok:
        sol[l0] = {l: x for l, x in val.items() if x}
print('determined labels:', {l: sol[l] for l in sorted(sol)})
print('still free:', [l for l in free if l not in sol])
new = {}
for j, v in vals.items():
    nv = {}
    for l, x in v.items():
        if l in sol:
            for l2, y in sol[l].items():
                nv[l2] = nv.get(l2, 0) + x * y
        else:
            nv[l] = nv.get(l, 0) + x
    new[j] = {l: x for l, x in nv.items() if x}
nbad = sum(1 for _ in residuals(new, lv))
print('after substitution: rows with nonzero residual %d' % nbad)
for l in free:
    if l in sol:
        js = [cls[j] for j, v in vals.items() if l in v and len(v) == 1]
        print('  %s = class %s' % (l, js))
if 'write' in sys.argv:
    save_values(q, cls, new)
    print('written')
