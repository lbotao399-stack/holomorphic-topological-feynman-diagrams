"""Solve the order-n nilpotency system from saved composition values (poolnilp.py --save).
usage: solve.py n kappa.json AUTS group1 [group2 ...]   group = comma-separated json files of one functional family
AUTS: 'name:aut,...' for the order-n graphs.  Lower-order kappa from kappa.json.  Prints rank, solution, residual,
and the residual of the ansatz kappa_n = sigma/|Aut| for sigma = +1, -1."""
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'qn-holomorphic'))
from nilp_explicit import rank_solve  # noqa: E402

n = int(sys.argv[1])
kap = {k: Fr(v) for k, v in json.load(open(sys.argv[2])).items()}
auts = {}
for it in sys.argv[3].split(','):
    a, b = it.split(':')
    auts[a] = int(b)
names = list(auts)
rows, rhs = [], []
P = None


def mq(x, P):
    return x.numerator % P * pow(x.denominator % P, P - 2, P) % P


def ratrec(a, P):
    a %= P
    r0, r1, s0, s1 = P, a, 0, 1
    while r1 * r1 * 2 > P:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


for grp in sys.argv[4:]:
    E = {}
    F = None
    for fn in grp.split(','):
        d = json.load(open(fn))
        P = d['P'] if P is None else P
        assert d['P'] == P
        F = d['F'] if F is None else F
        assert d['F'] == F
        E.update(d['E'])
    need = ['Q0|%s' % nm for nm in names] + ['%s|Q0' % nm for nm in names]
    missing = [k for k in need if k not in E]
    lower = [k for k in E if 'Q0' not in k.split('|')]
    print('group %s: F = %d, compositions %d, missing %s' % (grp, F, len(E), missing or 'none'))
    for fi in range(F):
        rows.append([(E.get('Q0|%s' % nm, [0] * F)[fi] + E.get('%s|Q0' % nm, [0] * F)[fi]) % P for nm in names])
        r = 0
        for k in lower:
            ka, kb = k.split('|')
            r = (r - mq(kap[ka] * kap[kb], P) * E[k][fi]) % P
        rhs.append(r)
    print('  lower compositions used: %s' % sorted(lower))
rk, aug, sol = rank_solve(rows, rhs, P)
print('order %d system: %d equations, %d unknowns, rank %d, augmented rank %d' % (n, len(rows), len(names), rk, aug))
nzcols = [nm for i, nm in enumerate(names) if any(r[i] for r in rows)]
print('nonzero columns: %s' % nzcols)
if sol is not None:
    for nm, s in zip(names, sol):
        k = ratrec(s, P)
        print('  kappa[%s] = %s   |Aut| = %d   kappa*|Aut| = %s' % (nm, k, auts[nm], k * auts[nm] if k is not None else None))
for sg in (1, -1):
    res = []
    for r, b in zip(rows, rhs):
        res.append((sum(mq(Fr(sg, auts[nm]), P) * r[i] for i, nm in enumerate(names)) - b) % P)
    print('ansatz kappa = %+d/|Aut|: residual %s' % (sg, res))
# rank of the homogeneous system restricted: which unknown directions are fixed
if sol is None and rk == aug:
    print('consistent, underdetermined: rank %d < %d' % (rk, len(names)))
elif rk != aug:
    print('INCONSISTENT')


def nullspace(rows, P):
    M = [list(r) for r in rows]
    nc = len(M[0])
    piv = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, len(M)) if M[i][c] % P), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        iv = pow(M[r][c], P - 2, P)
        M[r] = [x * iv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % P:
                f = M[i][c]
                M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * nc
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i][fc]) % P
        basis.append(v)
    return basis


ns = nullspace(rows, P)
print('null space of the coefficient matrix (%d vectors):' % len(ns))
for v in ns:
    print('  ', {nm: str(ratrec(x, P)) for nm, x in zip(names, v) if x})
