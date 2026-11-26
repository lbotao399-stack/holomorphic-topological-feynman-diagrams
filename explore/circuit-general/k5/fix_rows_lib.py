"""row residual count for a level-q table (all X + DS rows of rows<q>.txt)"""
import os
from fractions import Fraction as F
from solve_level import classes, load_values, SCR


def count_bad(q):
    cls = classes(q)
    low = classes(q - 1)
    vq = load_values(q)
    lv = load_values(q - 1)
    vals = [vq[b] for b in cls]
    nbad = 0
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
            vv = lv[low[kk[1]]] if isinstance(kk, tuple) else vals[kk]
            for lab, x in vv.items():
                acc[lab] = acc.get(lab, 0) + c * x
        if any(acc.values()):
            nbad += 1
    return nbad
