"""Exact periods of any graph whose series-parallel reduction ends in a core with loop
number <= 4 (K4 = W3, W4, ..., W_n, prism, K33), any integer powers, any q.

  period(edges, powers, q) -> {label: Fraction}, labels '1','z3','z5','z7',...
  kernel(nv, edges, lam, insertions) -> exact holomorphic kernel (T = 0)

Wheels: radial-angular residue algorithm (wheel.py).  Prism and K33: tables of all
convergent orbits at level q, solved from regulator-free relations (msys.py): triangle
rules (prism), diamond rules (K33), parametric IBP, W4 boundaries, star-triangle at unique
vertices.  Tables are cached in tables/<core>_q<q>.json.
"""
import json
import os
import time
from fractions import Fraction as F

import gutil
import msys

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, 'tables')
PRISM = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))
K33 = ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5))
FAMILIES = {'prism': PRISM, 'K33': K33}

_systems = {}
_tables = {}


def system(q):
    """a MultiSystem used only for reduction / family recognition at level q"""
    if q not in _systems:
        _systems[q] = msys.MultiSystem([PRISM, K33], q)
    return _systems[q]


def enc(label):
    if isinstance(label, str):
        return label
    return '*'.join('%s^%d' % (s, p) for s, p in label) or '1'


def dec(text, tuple_labels):
    if not tuple_labels:
        return text
    if text == '1':
        return ()
    out = []
    for part in text.split('*'):
        s, p = part.split('^')
        out.append((s, int(p)))
    return tuple(out)


def qname(q):
    q = F(q)
    return '%d' % q if q.denominator == 1 else '%do%d' % (q.numerator, q.denominator)


def table(name, q, verbose=True):
    q = F(q)
    if q.denominator == 1:
        q = int(q)
    key = (name, q)
    if key in _tables:
        return _tables[key]
    os.makedirs(TABLES, exist_ok=True)
    path = os.path.join(TABLES, '%s_q%s.json' % (name, qname(q)))
    tl = F(q).denominator != 1
    if os.path.exists(path):
        raw = json.load(open(path))
        tab = {tuple(r): {dec(k, tl): F(v) for k, v in val.items()} for r, val in raw['values']}
        _tables[key] = tab
        return tab
    t = time.time()
    S = msys.MultiSystem([FAMILIES[name]], q)
    S.build(verbose=verbose)
    solved, piv = S.solve(verbose=verbose)
    if len(solved) != len(S.reps):
        und = [S.reps[j] for j in range(len(S.reps)) if j not in solved]
        raise ArithmeticError('%s q=%d: %d undetermined orbits, e.g. %s' % (name, q, len(und), und[:3]))
    tab = {S.reps[j][1]: v for j, v in solved.items()}
    json.dump({'core': name, 'edges': FAMILIES[name], 'q': str(q), 'seconds': round(time.time() - t, 1),
               'values': [[list(r), {enc(k): str(x) for k, x in v.items()}] for r, v in sorted(tab.items())]},
              open(path, 'w'))
    _tables[key] = tab
    return tab


def period(edges, powers, q):
    S = system(q)
    res = S.classify(list(edges), list(powers), q)
    if res is None:
        raise ValueError('divergent or unrecognized core for powers %s' % (powers,))
    const, unk = res
    out = dict(const)
    for j, s in unk:
        f, rep = S.reps[j]
        name = 'prism' if f == 0 else 'K33'
        val = table(name, q)[tuple(rep)]
        for k, v in val.items():
            out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def kernel(nv, edges, lam, insertions=(), verbose=False, T=0):
    """exact kernel  K = int D prod y_ins  (htk conventions), T = 0 or 1 (zero kappa)."""
    import numer
    G = numer.Graph(nv, list(edges))
    t = time.time()
    if T == 0:
        poly, q = numer.numerator(G, lam, insertions, verbose)
    else:
        poly, q = numer.numerator_T1(G, lam, insertions, verbose)
    t1 = time.time()
    # orbit sums under the edge automorphisms of the graph (periods are invariant);
    # orbits with vanishing coefficient sum are dropped before any period is evaluated
    autos = gutil.edge_automorphisms(tuple(tuple(e) for e in edges))
    orb = {}
    for ex, c in poly.items():
        rep = gutil.orbit_rep(tuple(k + 1 for k in ex), autos)
        orb[rep] = orb.get(rep, 0) + c
    out = {}
    cache = {}
    for alpha, c in orb.items():
        if not c:
            continue
        if alpha not in cache:
            cache[alpha] = period(edges, alpha, q)
        for k, v in cache[alpha].items():
            out[k] = out.get(k, 0) + c * v
    if verbose:
        print('  numerator %d terms (%.1fs), periods %d (%.1fs)' % (len(poly), t1 - t, len(cache), time.time() - t1))
    return {k: v for k, v in out.items() if v}


_LOWER_KEYS = {}


def lower_core_value(edges, pw, q):
    """value of a prism or K33 period (edges isomorphic to the core) from the tables, else None"""
    if not _LOWER_KEYS:
        for name, fe in FAMILIES.items():
            _LOWER_KEYS[msys.struct_key(fe)] = name
    sk = msys.struct_key(edges)
    name = _LOWER_KEYS.get(sk)
    if name is None:
        return None
    S = system(q)
    res = S.classify(list(edges), list(pw), q)
    if res is None:
        return None
    const, unk = res
    out = dict(const)
    for j, s in unk:
        f, rep = S.reps[j]
        val = table('prism' if f == 0 else 'K33', q)[tuple(rep)]
        for k, v in val.items():
            out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def zigzag(n):
    """decompleted zigzag Z_n: circulant C_{n+2}(1,2) minus one vertex"""
    m = n + 2
    es = set()
    for i in range(m):
        for d in (1, 2):
            a, b = i, (i + d) % m
            if m - 1 in (a, b):
                continue
            es.add((min(a, b), max(a, b)))
    return tuple(sorted(es))
