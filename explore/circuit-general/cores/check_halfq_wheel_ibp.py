"""Half-integer q: W4 periods from the momentum-space reducer (halfq.wheel_value) satisfy the
parametric IBP relation (I) at q = 7/2 and the dimension relation (R1) between q = 5/2 and 7/2."""
from fractions import Fraction as F

import gutil
import halfq

W4 = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (1, 4))
cot = gutil.cotrees(W4)


def P(pw, q):
    return halfq.wheel_value(pw[:4], pw[4:], q)


def add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


q = F(7, 2)
n_ibp = 0
for b in [(1, 1, 1, 1, 1, 2, 1, 2), (2, 1, 1, 1, 1, 1, 2, 1), (1, 2, 1, 1, 2, 1, 1, 1)]:
    for i in range(8):
        if b[i] >= q - 1:
            continue
        lhs, rhs = {}, {}
        ok = True
        for C in cot:
            a = list(b)
            for e in C:
                a[e] += 1
            if not gutil.convergent(W4, a, q):
                ok = False
                break
            v = P(a, q)
            if i in C:
                rhs = add(rhs, v, (q - 1 - b[i]) / b[i])
            else:
                lhs = add(lhs, v)
        if ok:
            assert lhs == rhs, (b, i, lhs, rhs)
            n_ibp += 1
n_r1 = 0
for b in [(1, 1, 1, 1, 1, 2, 1, 2), (1, 1, 1, 2, 2, 1, 1, 1)]:
    if not gutil.convergent(W4, list(b), F(5, 2)):
        continue
    tot = {}
    for C in cot:
        a = list(b)
        for e in C:
            a[e] += 1
        tot = add(tot, P(a, F(7, 2)))
    assert tot == P(list(b), F(5, 2)), (b, tot)
    n_r1 += 1
print('W4 at half-integer q: %d IBP relations (q=7/2) and %d dimension relations (5/2 -> 7/2) hold exactly' % (n_ibp, n_r1))
