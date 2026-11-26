"""K5 period tables from the K6 class values:
   I_q(K5; alpha) = J_q(beta) prod_{e in K5} Gamma(beta_e) / Gamma(q),   beta = q - alpha on K5 = K6 - c.
All S5-orbits of convergent alpha in [1, q-1]^10 (sum 6q), K5 edge order 01,02,03,04,12,13,14,23,24,34.
Kinds: T1 (some beta_{v,oo} = 0), T2 (a triangle of K6 with sum beta = q), G (generic).
Output: tables/K5_q<q>.json  {'values': [[alpha, kind, {label: str}], ...]}
usage: python3 k5tables.py Q"""
import itertools
import json
import os
import sys
from fractions import Fraction as F
from math import factorial

import k6
from solve_level import load_values

HERE = os.path.dirname(os.path.abspath(__file__))
E5 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
P5 = []
for s in itertools.permutations(range(5)):
    idx = {frozenset(e): k for k, e in enumerate(E5)}
    P5.append([idx[frozenset((s[u], s[v]))] for (u, v) in E5])


def canon5(a):
    best = None
    for p in P5:
        b = [0] * 10
        for k in range(10):
            b[p[k]] = a[k]
        b = tuple(b)
        if best is None or b < best:
            best = b
    return best


def complete(a, q):
    """K6 weights (vertex 5 = oo) of the K5 vector alpha"""
    b = [0] * 15
    for k, (u, v) in enumerate(E5):
        b[k6.EIDX[u, v]] = q - a[k]
    for v in range(5):
        s = sum(b[k6.EIDX[v, w]] for w in range(5) if w != v)
        b[k6.EIDX[v, 5]] = 2 * q - s
    return b


def orbits(q):
    vals = load_values(q)
    out = {}
    for c6, v in vals.items():
        for c in k6.k5_decompletions(c6, q):
            verts = [x for x in range(6) if x != c]
            # map verts -> 0..4 (any order; S5-canonicalized below)
            a = [q - c6[k6.EIDX[verts[u], verts[w]]] for (u, w) in E5]
            key = canon5(a)
            if key not in out:
                out[key] = c6
    return out, vals


def I_from_J(a, Jv, q):
    n = F(1, factorial(q - 1))
    for x in a:
        n *= factorial(q - x - 1)
    return {k: x * n for k, x in Jv.items() if x}


if __name__ == '__main__':
    q = int(sys.argv[1])
    orb, vals = orbits(q)
    rows = []
    cnt = {}
    for a in sorted(orb):
        b6 = orb[a]
        kd = k6.kind(b6, q)
        cnt[kd] = cnt.get(kd, 0) + 1
        rows.append([list(a), kd, {k: str(x) for k, x in I_from_J(a, vals[b6], q).items()}])
    os.makedirs(os.path.join(HERE, 'tables'), exist_ok=True)
    json.dump({'core': 'K5', 'edges': E5, 'q': q, 'values': rows}, open(os.path.join(HERE, 'tables', 'K5_q%d.json' % q), 'w'))
    labs = set()
    for r in rows:
        labs |= set(r[2])
    print('q=%d K5 orbits %d: %s; labels %s' % (q, len(rows), cnt, sorted(labs, key=lambda s: (len(s), s))))
