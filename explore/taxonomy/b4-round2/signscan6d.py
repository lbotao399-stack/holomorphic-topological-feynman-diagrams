"""6d (1,0) sl2 + 16 fund (K-interpolated), order 2: residual for kappa_G = s_G/|Aut G|, s in {0,+-1}^5 (0 = graph dropped)."""
import json
import itertools
from fractions import Fraction as Fr
names = ['G2.1', 'G2.2', 'G2.3u', 'G2.4u', 'G2.5u']
aut = dict(zip(names, [2, 2, 1, 2, 2]))
rows, rhs = [], []
for tag in ('gvgvuJ', 'gv3J'):
    d = json.load(open('o2_6d10/E_f16interp_%s.json' % tag))
    P, F, E = d['P'], d['F'], d['E']
    for f in range(F):
        rows.append([(E['Q0|' + nm][f] + E[nm + '|Q0'][f]) % P for nm in names])
        rhs.append((-pow(4, P - 2, P) * E['G1.1|G1.1'][f]) % P)
mq = lambda x: x.numerator % P * pow(x.denominator % P, P - 2, P) % P
good = []
for sg in itertools.product((0, 1, -1), repeat=5):
    res = [(sum(mq(Fr(s, aut[nm])) * r[i] for i, (s, nm) in enumerate(zip(sg, names))) - b) % P for r, b in zip(rows, rhs)]
    if not any(res):
        good.append(sg)
print('patterns s in {0,+1,-1}^5 with residual 0 on %d equations: %d of 243' % (len(rows), len(good)))
for g in good:
    print('  ', dict(zip(names, g)), ' 15k3+6k4+4k5 =', 15 * Fr(g[2], 1) + 6 * Fr(g[3], 2) + 4 * Fr(g[4], 2))
