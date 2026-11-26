"""Boundary values J_q(beta) = Pi(K6; beta) for K6 classes, via
  (E) symmetry engine on the completed graph (symmetry/engine.py: solve), numerator terminal on;
  (P) per5 (fiveloop): a decompletion c with beta in [1, q-1] on K6 - c, then K5 (Y)->K33+e, K5-e (dual prism);
  (S) graphical functions (graphical/seq.py) on a decompletion with alpha in [1, q-1] (T1 sequential plans).
Values: dict label -> Fraction, labels '1', 'z3', 'z5', 'z3*z3', ...
"""
import os
import sys
from fractions import Fraction as F
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
CG = os.path.abspath(os.path.join(HERE, '..'))
for sub in ('graphical', 'fiveloop', 'symmetry', 'cores'):
    p = os.path.join(CG, sub)
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)      # cores first (its wheel.py shadows the site-packages 'wheel')

import k6  # noqa: E402

_engine = None


def engine():
    global _engine
    if _engine is None:
        import engine as en
        _engine = en
    return _engine


def lab_of_mon(mon):
    if not mon:
        return '1'
    return '*'.join(mon)


def from_engine(v):
    return {lab_of_mon(m): x for m, x in v.items() if x}


def via_engine(b, q, budget=60, numer=3):
    en = engine()
    E = {}
    for k, (i, j) in enumerate(k6.EDGES):
        if b[k]:
            E[(i, j)] = b[k]
    old = en.NUMER[0]
    en.NUMER[0] = numer
    try:
        r = en.solve(6, E, q, budget=budget)
    except (ArithmeticError, ZeroDivisionError, ValueError, KeyError, AssertionError):
        r = None
    finally:
        en.NUMER[0] = old
    if r is None:
        return None
    return from_engine(r[0]), 'E:' + ' '.join(r[1])


def decomp(b, q, c):
    """K6 - c as momentum-space K5: (edges, alpha, beta)"""
    verts = [v for v in range(6) if v != c]
    lab = {v: i for i, v in enumerate(verts)}
    edges, al, be = [], [], []
    for k, (i, j) in enumerate(k6.EDGES):
        if c in (i, j):
            continue
        edges.append((lab[i], lab[j]))
        al.append(q - b[k])
        be.append(b[k])
    return edges, al, be


def norm_J(be, q):
    n = F(factorial(q - 1))
    for x in be:
        n /= factorial(x - 1)
    return n


_per5 = None


def via_per5(b, q):
    global _per5
    if _per5 is None:
        import per5
        _per5 = per5
    for c in k6.k5_decompletions(b, q):
        edges, al, be = decomp(b, q, c)
        try:
            v = _per5.value(edges, al, q)
        except Exception:
            v = None
        if v is not None:
            n = norm_J(be, q)
            return {k: x * n for k, x in v.items() if x}, 'P:c=%d' % c
    return None


_seq = None


def _load_seq():
    global _seq
    if _seq is None:
        import seq
        import mzv
        from validate import to_labels
        _seq = (seq, mzv, to_labels)
    return _seq


def seq_plans(b, q):
    """sequential graphical-function plans for K6 - c (zero edges dropped, hub edges any integer weight)"""
    seq, mzv, to_labels = _load_seq()
    out = []
    for c in range(6):
        g2 = {}
        for k, (i, j) in enumerate(k6.EDGES):
            if c in (i, j) or b[k] == 0:
                continue
            g2[frozenset((i, j))] = b[k]
        for h0, h1, path in seq.decompositions(g2):
            pl = seq.plan(g2, q, h0, h1, path)
            if pl is not None:
                pl['route'] = 'c=%d' % c
                out.append(pl)
    out.sort(key=lambda pl: (pl['appends'], str(pl['route']), str(pl['path'])))
    return out


def via_seq(b, q, choice=0):
    """J = A_{K6-c} directly (graphical/seq.run returns the position-space period A_G = J)"""
    seq, mzv, to_labels = _load_seq()
    pls = seq_plans(b, q)
    if len(pls) <= choice:
        return None
    pl = pls[choice]
    A = seq.run(pl)
    v = to_labels(A)
    return {k: F(x) for k, x in v.items() if x}, 'S:%s h=%s path=%s app=%d' % (pl['route'], pl['hubs'], pl['path'],
                                                                              pl['appends'])


def value(b, q, methods=('E', 'P', 'S'), budget=60):
    for m in methods:
        if m == 'E':
            r = via_engine(b, q, budget=budget)
        elif m == 'P':
            r = via_per5(b, q)
        else:
            r = via_seq(b, q)
        if r is not None:
            return r
    return None
