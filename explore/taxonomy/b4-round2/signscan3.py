"""8d order 3, O = gamma^4, functionals 0..2: residual for every sign pattern kappa_G = s_G/|Aut G| (4096) and for kappa_3 = 0."""
import json
import itertools
from fractions import Fraction as Fr
d = json.load(open('o3g4/A_all.json'))
P, F, E = d['P'], d['F'], d['E']
kap = {k: Fr(v) for k, v in json.load(open('../qn-holomorphic/b_kappa_8d_so5.json')).items()}
aut = dict(zip(['G3.%d' % i for i in range(1, 13)], [2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 8, 1]))
mq = lambda x: x.numerator % P * pow(x.denominator % P, P - 2, P) % P
names = list(aut)
rows = [[(E['Q0|' + nm][f] + E[nm + '|Q0'][f]) % P for nm in names] for f in range(F)]
rhs = [sum(-mq(kap[k.split('|')[0]] * kap[k.split('|')[1]]) * E[k][f] for k in E if 'Q0' not in k) % P for f in range(F)]
print('kappa_3 = 0: residual', [(-b) % P for b in rhs])
good = []
for sg in itertools.product((1, -1), repeat=12):
    res = [(sum(mq(Fr(s, aut[nm])) * r[i] for i, (s, nm) in enumerate(zip(sg, names))) - b) % P for r, b in zip(rows, rhs)]
    if not any(res):
        good.append(sg)
print('sign patterns s in {+-1}^12 with residual 0 on all %d functionals: %d of 4096' % (F, len(good)))
for g in good:
    print('  ', dict(zip(names, g)))
