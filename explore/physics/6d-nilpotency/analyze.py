"""Residuals of sum kappa kappa' E from a run_nilp log, and the same sum with one weight altered.
Usage: python3 analyze.py LOG n [prime]"""
import re
import sys
from fractions import Fraction as Fr

sys.argv += [] if len(sys.argv) > 3 else ['2147483647']
log, n, P = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
import kern
kern.set_prime(P)
from kern import all_graphs, modq

E = {}
for line in open(log):
    m = re.match(r'\s+<hatQ\[(\S+)\] hatQ\[(\S+)\] O>: \[([0-9, ]+)\]', line)
    if m:
        E[(m.group(1), m.group(2))] = [int(x) for x in m.group(3).split(',')]
F = len(next(iter(E.values())))
G = {m: all_graphs(m) for m in range(1, n + 1)}
kap = {'Q0': Fr(1)}
for m in G:
    for g in G[m]:
        kap[g.name] = g.kappa


def resid(k):
    out = []
    for fi in range(F):
        s = 0
        for (a, b), v in E.items():
            s = (s + modq(k[a]) * modq(k[b]) % P * v[fi]) % P
        out.append(s)
    return out


print('compositions read: %d, functional sets: %d' % (len(E), F))
print('kappa = wick/|Aut|: residual', resid(kap))
nonzero = sum(1 for v in E.values() for x in v if x % P)
print('nonzero composition values: %d of %d' % (nonzero, len(E) * F))
for g in G[n]:
    k2 = dict(kap)
    k2[g.name] = -kap[g.name]
    r1 = resid(k2)
    k2[g.name] = Fr(0)
    r2 = resid(k2)
    print('  %s kappa %s: flipped -> zero residual %s; dropped -> zero residual %s'
          % (g.name, kap[g.name], all(x == 0 for x in r1), all(x == 0 for x in r2)))
k2 = {a: abs(v) for a, v in kap.items()}
print('kappa = +1/|Aut| (no Wick sign): zero residual', all(x == 0 for x in resid(k2)))
