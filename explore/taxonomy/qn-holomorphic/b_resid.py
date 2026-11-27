"""Fork B: residual sum kappa kappa' E of an order-n log of run_nilp_f.py with prescribed kappa (default 1/|Aut|),
and the rank of the kappa_n system augmented by the prescribed solution.  usage: python3 b_resid.py LOG [P]"""
import re, sys
from fractions import Fraction as Fr
P = int(sys.argv[2]) if len(sys.argv) > 2 else 2147483647
txt = open(sys.argv[1]).read()
aut = {}
for m in re.finditer(r"\('(G[\d.]+u?)', (\d+), \[", txt):
    aut[m.group(1)] = int(m.group(2))
E = {}
for m in re.finditer(r"<hatQ\[(\S+)\] hatQ\[(\S+)\] O>: \[([\d, ]*)\]", txt):
    E[(m.group(1), m.group(2))] = [int(x) for x in m.group(3).split(',')]
kap = {'Q0': 1}
for nm, a in aut.items():
    kap[nm] = Fr(1, a)
def md(x):
    return x.numerator % P * pow(x.denominator, P - 2, P) % P
F = len(next(iter(E.values())))
res = [0] * F
for (a, b), v in E.items():
    k = md(kap[a]) * md(kap[b]) % P
    for i in range(F):
        res[i] = (res[i] + k * v[i]) % P
print('kappa = 1/|Aut|:', {k: str(v) for k, v in kap.items() if k != 'Q0'})
print('residual sum kappa kappa E with kappa = 1/|Aut|:', res)
import itertools
n = max(int(k[1]) for k in aut)
top = sorted(k for k in aut if k.startswith('G%d.' % n))
print('sign patterns at order %d (lower orders kappa = 1/|Aut|):' % n)
for signs in itertools.product((1, -1), repeat=len(top)):
    kk = dict(kap)
    for s, nm in zip(signs, top):
        kk[nm] = s * Fr(1, aut[nm])
    r = [0] * F
    for (a, b), v in E.items():
        k = md(kk[a]) * md(kk[b]) % P
        for i in range(F):
            r[i] = (r[i] + k * v[i]) % P
    if not any(r):
        print('  residual 0 for kappa = %s' % {nm: str(kk[nm]) for nm in top})
