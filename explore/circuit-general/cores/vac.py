"""Exact vacuum periods in integer dimension 2q without regulator.

Normalized period  p(G; a) = P(G; a) / prod_e Gamma(a_e),  P = int prod t^(a-1) / U^q
(projective, sum a = n q).  Rules, each used only if every term is a convergent period:

  series  (bivalent vertex):     p unchanged, powers add
  parallel (two edges u-v):      p = G_2q(s1,s2) * p(merged, s1+s2-q),
       G_2q(a,b) = Gamma(q-a)Gamma(q-b)Gamma(a+b-q)/(Gamma(a)Gamma(b)Gamma(2q-a-b))
  contraction (a_e = 0):         p(G; .., 0, ..) = p(G/e; ..)
  triangle rule (A,B cubic, L1=AB, L2=BC, L3=AC, L4 third at B, L5 third at A):
     (2q - 2a1 - a2 - a3) p(a) = a2 [p(a-e1+e2) - p(a+e2-e4)] + a3 [p(a-e1+e3) - p(a+e3-e5)]
  wheels W_n (n>=3, K4 = W3) by wheel.integrate.
A marked edge (the virtual cut) is never A-B, L4 or L5; this makes the sum of the powers
of the lowered edges a strictly decreasing potential.
Values: dicts {label: Fraction}, labels '1', 'z3', 'z5', ...
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
from math import factorial, prod

import gutil
import wheel


class Stuck(Exception):
    pass


def add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def scale(a, s):
    return {k: v * s for k, v in a.items() if v * s}


def canon(edges, powers, marks):
    """canonical key: sorted tuple of (u, v, power, mark) under vertex relabelling"""
    verts = sorted({x for e in edges for x in e})
    inc = {v: [] for v in verts}
    for (u, v), p, m in zip(edges, powers, marks):
        inc[u].append((p, m))
        inc[v].append((p, m))
    sig = {v: (len(inc[v]), tuple(sorted(inc[v]))) for v in verts}
    groups = {}
    for v in verts:
        groups.setdefault(sig[v], []).append(v)
    keys = sorted(groups)

    def rec(i, acc):
        if i == len(keys):
            yield acc
            return
        for perm in permutations(groups[keys[i]]):
            yield from rec(i + 1, acc + list(perm))
    best = None
    for order in rec(0, []):
        mp = {v: j for j, v in enumerate(order)}
        key = tuple(sorted((min(mp[u], mp[v]), max(mp[u], mp[v]), p, m)
                           for (u, v), p, m in zip(edges, powers, marks)))
        if best is None or key < best:
            best = key
    return best


def G2q(a, b, q):
    return F(factorial(q - a - 1) * factorial(q - b - 1) * factorial(a + b - q - 1),
             factorial(a - 1) * factorial(b - 1) * factorial(2 * q - a - b - 1))


STATS = {'nodes': 0, 'triangle': 0, 'wheel': 0}


def period(edges, powers, q, mark=0):
    """normalized period p(G; a); raises Stuck if no valid rule applies."""
    marks = [int(i == mark) for i in range(len(edges))]
    return _p(canon(edges, powers, marks), q)


def full_period(edges, powers, q, mark=0):
    p = period(edges, powers, q, mark)
    norm = prod(factorial(x - 1) for x in powers)
    return {k: v * norm for k, v in p.items()}


def unpack(key):
    return [(u, v) for u, v, p, m in key], [p for u, v, p, m in key], [m for u, v, p, m in key]


def contract(edges, pw, mk, i):
    u, v = edges[i]
    ne, npw, nmk = [], [], []
    for j, ((a, b), p, m) in enumerate(zip(edges, pw, mk)):
        if j == i:
            continue
        ne.append((u if a == v else a, u if b == v else b))
        npw.append(p)
        nmk.append(m)
    return ne, npw, nmk


def valid(edges, pw, q):
    if any(p < 0 for p in pw):
        return False
    ne, npw = list(edges), list(pw)
    while 0 in npw:
        i = npw.index(0)
        u, v = ne[i]
        if u == v:
            return False
        ne, npw, _ = contract(ne, npw, [0] * len(npw), i)
    if any(a == b for a, b in ne):
        return False
    return gutil.convergent(tuple(ne), npw, q)


@lru_cache(maxsize=None)
def _p(key, q):
    STATS['nodes'] += 1
    edges, pw, mk = unpack(key)
    if any(p == 0 for p in pw):
        i = pw.index(0)
        if edges[i][0] == edges[i][1]:
            raise Stuck('zero self-loop')
        ne, npw, nmk = contract(edges, pw, mk, i)
        return _p(canon(ne, npw, nmk), q)
    if any(u == v for u, v in edges):
        raise Stuck('self-loop')
    if not gutil.convergent(tuple(edges), pw, q):
        raise Stuck('divergent %s' % (key,))
    verts = sorted({x for e in edges for x in e})
    deg = {v: 0 for v in verts}
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    # series
    for v in verts:
        if deg[v] == 2:
            i, j = [t for t, e in enumerate(edges) if v in e]
            a = edges[i][0] if edges[i][1] == v else edges[i][1]
            b = edges[j][0] if edges[j][1] == v else edges[j][1]
            if a == b:
                raise Stuck('series to self-loop')
            keep = [t for t in range(len(edges)) if t not in (i, j)]
            ne = [edges[t] for t in keep] + [(a, b)]
            npw = [pw[t] for t in keep] + [pw[i] + pw[j]]
            nmk = [mk[t] for t in keep] + [max(mk[i], mk[j])]
            return _p(canon(ne, npw, nmk), q)
    # parallel
    seen = {}
    for i, (u, v) in enumerate(edges):
        pr = (min(u, v), max(u, v))
        if pr in seen:
            j = seen[pr]
            s1, s2 = pw[j], pw[i]
            if not (s1 + s2 >= q + 1):
                raise Stuck('divergent bubble')
            keep = [t for t in range(len(edges)) if t not in (i, j)]
            ne = [edges[t] for t in keep] + [pr]
            npw = [pw[t] for t in keep] + [s1 + s2 - q]
            nmk = [mk[t] for t in keep] + [max(mk[i], mk[j])]
            if len(ne) == 1:
                # one-loop vacuum banana closed by the merge: two edges between two vertices
                raise Stuck('banana')
            return scale(_p(canon(ne, npw, nmk), q), G2q(s1, s2, q))
        seen[pr] = i
    if len(verts) == 2:
        raise Stuck('two-vertex graph')
    w = wheel_form(edges, pw, verts, deg)
    if w is not None:
        spoke, rim = w
        STATS['wheel'] += 1
        val = wheel.integrate(spoke, rim, q)
        norm = prod(factorial(p - 1) for p in pw)
        return {('1' if k == 0 else 'z%d' % k): v / norm for k, v in val.items() if v}
    # triangle rule
    marked = [t for t in range(len(edges)) if mk[t]]
    forbidden = set()
    for t in marked:
        forbidden |= set(edges[t])
    best = None
    for i, (A0, B0) in enumerate(edges):
        if mk[i]:
            continue
        for A, B in ((A0, B0), (B0, A0)):
            if A in forbidden or B in forbidden or deg[A] != 3 or deg[B] != 3:
                continue
            incA = [t for t, e in enumerate(edges) if A in e and t != i]
            incB = [t for t, e in enumerate(edges) if B in e and t != i]
            for t3 in incA:
                C = edges[t3][0] if edges[t3][1] == A else edges[t3][1]
                for t2 in incB:
                    C2 = edges[t2][0] if edges[t2][1] == B else edges[t2][1]
                    if C2 != C or C in (A, B):
                        continue
                    t5 = [t for t in incA if t != t3][0]
                    t4 = [t for t in incB if t != t2][0]
                    if t4 == t5 or mk[t4] or mk[t5]:
                        continue
                    a1, a2, a3 = pw[i], pw[t2], pw[t3]
                    c0 = 2 * q - 2 * a1 - a2 - a3
                    if c0 == 0:
                        continue
                    terms = []
                    for (x, dx), (y, dy), coef in (((i, -1), (t2, 1), a2), ((t2, 1), (t4, -1), -a2),
                                                   ((i, -1), (t3, 1), a3), ((t3, 1), (t5, -1), -a3)):
                        np_ = list(pw)
                        np_[x] += dx
                        np_[y] += dy
                        terms.append((np_, coef))
                    if all(valid(edges, np_, q) for np_, c in terms):
                        cost = (a1 + pw[t4] + pw[t5], abs(c0))
                        if best is None or cost < best[0]:
                            best = (cost, terms, c0)
    if best is None:
        raise Stuck('no valid triangle %s' % (key,))
    STATS['triangle'] += 1
    _, terms, c0 = best
    out = {}
    for np_, c in terms:
        out = add(out, scale(_p(canon(edges, np_, mk), q), F(c)))
    return scale(out, F(1, c0))


def wheel_form(edges, pw, verts, deg):
    """recognize a simple wheel W_n: hub of degree n, rim a cycle of cubic vertices"""
    n = len(verts) - 1
    if n < 3 or len(edges) != 2 * n:
        return None
    if len({(min(u, v), max(u, v)) for u, v in edges}) != len(edges):
        return None
    if n == 3:
        if not all(deg[v] == 3 for v in verts):
            return None
        hub = verts[0]
    else:
        hubs = [v for v in verts if deg[v] == n]
        if len(hubs) != 1 or any(deg[w] != 3 for w in verts if w != hubs[0]):
            return None
        hub = hubs[0]
    rimv = [v for v in verts if v != hub]
    eidx = {frozenset(e): t for t, e in enumerate(edges)}
    if any(frozenset((hub, r)) not in eidx for r in rimv):
        return None
    adj = {r: [] for r in rimv}
    for t, (u, v) in enumerate(edges):
        if hub in (u, v):
            continue
        adj[u].append(v)
        adj[v].append(u)
    if any(len(adj[r]) != 2 for r in rimv):
        return None
    order = [rimv[0], adj[rimv[0]][0]]
    while len(order) < n:
        nxt = [w for w in adj[order[-1]] if w != order[-2]][0]
        if nxt == order[0]:
            return None
        order.append(nxt)
    if order[0] not in adj[order[-1]]:
        return None
    spoke = [pw[eidx[frozenset((hub, r))]] for r in order]
    rim = [pw[eidx[frozenset((order[k], order[(k + 1) % n]))]] for k in range(n)]
    return spoke, rim
