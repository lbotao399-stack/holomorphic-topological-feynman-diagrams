"""Massless propagator integrals in d = 2q - 2 eps by the triangle rule.

Graph: lines (u, v, n, k) = propagator 1/(p^2)^(n + k eps) between vertices u, v;
external momentum Q (Q^2 = 1) enters vertex 0 and leaves vertex 1.  Measure
d^d p / pi^(d/2) per loop.  Reductions: deletion of zero lines, series merge at
bivalent vertices, one-loop bubbles G(a,b), bridges and cut vertices, and the
triangle rule

  (d - 2 a1 - a2 - a3) I = a2 [I(1-,2+) - I(2+,4-)] + a3 [I(1-,3+) - I(3+,5-)],

L1 = A-B, L2 = B-C, L3 = A-C, L4 = third line at B, L5 = third line at A.
"""
from functools import lru_cache
from itertools import permutations
import sys

from fractions import Fraction
from laurent import (Series, linear, gamma_value, vmul, vadd, vscale_series, value_const)
F1 = Fraction(1)
F0 = Fraction(0)

sys.setrecursionlimit(100000)


class NotReducible(Exception):
    pass


def is_int(l):
    return l[3] == 0 and Fraction(l[2]).denominator == 1


Q = None   # set by set_q


def set_q(q):
    global Q
    Q = q
    reduce_graph.cache_clear()


def G_value(a, b):
    """G(a,b) = Gamma(d/2-a)Gamma(d/2-b)Gamma(a+b-d/2)/(Gamma(a)Gamma(b)Gamma(d-a-b)),
    a=(n1,k1), b=(n2,k2), d/2 = Q - eps."""
    n1, k1 = a
    n2, k2 = b
    v = gamma_value(Q - n1, -1 - k1)
    v = vmul(v, gamma_value(Q - n2, -1 - k2))
    v = vmul(v, gamma_value(n1 + n2 - Q, k1 + k2 + 1))
    v = vmul(v, gamma_value(n1, k1, -1))
    v = vmul(v, gamma_value(n2, k2, -1))
    v = vmul(v, gamma_value(2 * Q - n1 - n2, -2 - k1 - k2, -1))
    return v


def canonical(lines):
    """lines: iterable of (u,v,n,k). Returns canonical tuple key and nothing else."""
    lines = [tuple(l) for l in lines]
    verts = sorted({x for l in lines for x in l[:2]} | {0, 1})
    internal = [v for v in verts if v not in (0, 1)]
    best = None
    for swap in (False, True):
        for perm in permutations(range(2, 2 + len(internal))):
            mp = {0: 1 if swap else 0, 1: 0 if swap else 1}
            for v, w in zip(internal, perm):
                mp[v] = w
            key = tuple(sorted((min(mp[u], mp[v]), max(mp[u], mp[v]), n, k) for u, v, n, k in lines))
            if best is None or key < best:
                best = key
    return best


def components(lines, verts):
    adj = {v: set() for v in verts}
    for u, v, n, k in lines:
        adj[u].add(v)
        adj[v].add(u)
    seen = {}
    comps = []
    for v in verts:
        if v in seen:
            continue
        stack = [v]
        seen[v] = len(comps)
        comp = {v}
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen[y] = len(comps)
                    comp.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


def relabel(sub, a, b):
    """relabel vertices so that a -> 0 (in), b -> 1 (out), others 2,3,..."""
    vs = sorted({x for l in sub for x in l[:2]} | {a, b})
    mp = {a: 0, b: 1}
    nxt = 2
    for w in vs:
        if w not in mp:
            mp[w] = nxt
            nxt += 1
    return [(mp[x[0]], mp[x[1]], x[2], x[3]) for x in sub]


def articulation_split(lines, verts):
    """Return None if no cut vertex; 'zero' if a tadpole block exists; otherwise
    a list of (sublines, in, out) blocks along the 0-1 path for one cut vertex."""
    for c in verts:
        others = [v for v in verts if v != c]
        sub = [l for l in lines if c not in l[:2]]
        comps = components(sub, others)
        if len(comps) <= 1:
            continue
        if c in (0, 1):
            return 'zero'          # a block hanging at an external vertex without the other one
        c0 = [x for x in comps if 0 in x]
        c1 = [x for x in comps if 1 in x]
        if len(comps) != 2 or not c0 or not c1 or c0[0] is c1[0]:
            return 'zero'
        blocks = []
        for comp, a, b in ((c0[0], 0, c), (c1[0], c, 1)):
            bl = [l for l in lines if l[0] in comp or l[1] in comp]
            blocks.append((bl, a, b))
        return blocks
    return None


STAR_KITES = {'count': 0, 'unknown': 0}


def star_kite(lines):
    """Kite (two-loop propagator) whose only non-integer line is the central one,
    K(a,b,c,d; n + k eps).  If K(a,b,c,d; n) is finite at eps=0 then
    K(...; n + k eps) = K(...; n) + eps * U  with U regular: U is carried in the
    channel ('U',) so that any eps^0 dependence on it is detected."""
    if len(lines) != 5:
        return None
    verts = sorted({x for l in lines for x in l[:2]})
    if len(verts) != 4:
        return None
    internal = [v for v in verts if v not in (0, 1)]
    if len(internal) != 2:
        return None
    X, Y = internal
    mid = [l for l in lines if {l[0], l[1]} == {X, Y}]
    if len(mid) != 1:
        return None
    mid = mid[0]
    outer = [l for l in lines if l is not mid]
    if any(not is_int(l) for l in outer) or is_int(mid):
        return None
    if mid[2] < 0:
        return None
    base = [l if l is not mid else (l[0], l[1], l[2], 0) for l in lines]
    v = ev(base)
    STAR_KITES['count'] += 1
    from laurent import finite_part
    try:
        finite_part(v)
    except ArithmeticError:
        STAR_KITES['unknown'] += 1
        raise NotReducible('t1star kite with divergent integer limit: %s' % (lines,))
    out = dict(v)
    unk = Series(1, [F1] + [F0] * 13, 15)
    out[((999999, 1),)] = unk
    return out


@lru_cache(maxsize=None)
def reduce_graph(key):
    lines = [list(l) for l in key]
    return _reduce(lines)


def ev(lines):
    lines = [tuple(l) for l in lines]
    if not lines:
        raise NotReducible("empty graph")
    return reduce_graph(canonical(lines))


def _reduce(lines):
    # 1. zero lines: momentum space -> contract the line (identify its endpoints)
    while True:
        z = [l for l in lines if l[2] == 0 and l[3] == 0]
        if not z:
            break
        u, v = z[0][0], z[0][1]
        if u == v:
            return {}
        if {u, v} == {0, 1}:
            return {}           # external vertices identified: scaleless
        keep, drop = (u, v) if u in (0, 1) or (v not in (0, 1) and u < v) else (v, u)
        rest = [l for l in lines if l is not z[0]]
        lines = [(keep if a == drop else a, keep if b == drop else b, n, k) for a, b, n, k in rest]
    verts = sorted({x for l in lines for x in l[:2]} | {0, 1})
    if not lines:
        return {}
    # self loops -> scaleless
    if any(l[0] == l[1] for l in lines):
        return {}
    deg = {v: 0 for v in verts}
    for u, v, n, k in lines:
        deg[u] += 1
        deg[v] += 1
    if deg[0] == 0 or deg[1] == 0:
        return {}
    for v in verts:
        if v not in (0, 1) and deg[v] <= 1:
            return {}
    # 2. series merge at bivalent internal vertex
    for v in verts:
        if v in (0, 1) or deg[v] != 2:
            continue
        inc = [l for l in lines if v in l[:2]]
        a, b = inc
        ea = a[0] if a[1] == v else a[1]
        eb = b[0] if b[1] == v else b[1]
        if ea == eb:
            return {}       # tadpole bubble hanging at ea
        rest = [l for l in lines if l is not a and l is not b]
        rest.append((ea, eb, a[2] + b[2], a[3] + b[3]))
        return ev(rest)
    # 3. parallel lines -> bubble
    seen = {}
    for i, l in enumerate(lines):
        pair = (min(l[0], l[1]), max(l[0], l[1]))
        if pair in seen:
            j = seen[pair]
            a, b = lines[j], lines[i]
            rest = [x for t, x in enumerate(lines) if t not in (i, j)]
            na, ka = a[2], a[3]
            nb, kb = b[2], b[3]
            rest.append((pair[0], pair[1], na + nb - Q, ka + kb + 1))
            g = G_value((na, ka), (nb, kb))
            return vmul(g, ev(rest))
        seen[pair] = i
    # 4. connectivity
    comps = components(lines, verts)
    if len(comps) > 1:
        return {}
    # base case
    if len(lines) == 1:
        u, v, n, k = lines[0]
        if {u, v} == {0, 1}:
            return value_const(1)
        return {}
    # 5. bridges
    for i, l in enumerate(lines):
        rest = lines[:i] + lines[i + 1:]
        cs = components(rest, verts)
        if len(cs) > 1:
            u, v = l[0], l[1]
            cu = [c for c in cs if u in c][0]
            cv = [c for c in cs if v in c][0]
            if 0 in cu and 1 in cv:
                left_in, left_out = 0, u
                right_in, right_out = v, 1
                cl, cr = cu, cv
            elif 1 in cu and 0 in cv:
                left_in, left_out = 0, v
                right_in, right_out = u, 1
                cl, cr = cv, cu
            else:
                return {}
            res = value_const(1)
            for comp, a, b in ((cl, left_in, left_out), (cr, right_in, right_out)):
                sub = [x for x in rest if x[0] in comp]
                if a == b:
                    if sub:
                        return {}       # tadpole
                    continue
                mp = {}
                nxt = 2
                for w in sorted(comp):
                    if w == a:
                        mp[w] = 0
                    elif w == b:
                        mp[w] = 1
                    else:
                        mp[w] = nxt
                        nxt += 1
                sub = [(mp[x[0]], mp[x[1]], x[2], x[3]) for x in sub]
                res = vmul(res, ev(sub))
            return res
    # 6. cut vertices
    sp = articulation_split(lines, verts)
    if sp == 'zero':
        return {}
    if sp is not None:
        res = value_const(1)
        for bl, a, b in sp:
            res = vmul(res, ev(relabel(bl, a, b)))
        return res
    # 7. triangle rule
    best = None
    lineset = lines
    for i, L1 in enumerate(lineset):
        if not is_int(L1) or L1[2] < 1:
            continue
        for A, B in ((L1[0], L1[1]), (L1[1], L1[0])):
            if A in (0, 1) or B in (0, 1) or deg[A] != 3 or deg[B] != 3:
                continue
            others_A = [l for t, l in enumerate(lineset) if t != i and A in l[:2]]
            others_B = [l for t, l in enumerate(lineset) if t != i and B in l[:2]]
            for L3 in others_A:
                C = L3[0] if L3[1] == A else L3[1]
                L2s = [l for l in others_B if C in l[:2]]
                if not L2s:
                    continue
                L2 = L2s[0]
                L5 = [l for l in others_A if l is not L3][0]
                L4 = [l for l in others_B if l is not L2][0]
                if not is_int(L4) or not is_int(L5) or L4[2] < 1 or L5[2] < 1:
                    continue
                if L4 is L5:
                    continue
                c0 = 2 * Q - 2 * L1[2] - L2[2] - L3[2]
                c1 = -2 - L2[3] - L3[3]
                if c0 == 0 and c1 == 0:
                    continue
                cost = (L1[2] + L4[2] + L5[2], c0 == 0)
                if best is None or cost < best[0]:
                    best = (cost, L1, L2, L3, L4, L5, c0, c1)
    if best is None:
        kite = star_kite(lines)
        if kite is not None:
            return kite
        raise NotReducible(str(lines))
    _, L1, L2, L3, L4, L5, c0, c1 = best

    def mod(changes):
        out = []
        for l in lines:
            d = 0
            for target, delta in changes:
                if l is target:
                    d += delta
            out.append((l[0], l[1], l[2] + d, l[3]))
        return out

    a2n, a2k = L2[2], L2[3]
    a3n, a3k = L3[2], L3[3]
    res = {}
    # alpha2 terms
    t = vadd(ev(mod([(L1, -1), (L2, +1)])), ev(mod([(L2, +1), (L4, -1)])), -1)
    res = vadd(res, vscale_series(t, linear(a2n, a2k)))
    t = vadd(ev(mod([(L1, -1), (L3, +1)])), ev(mod([(L3, +1), (L5, -1)])), -1)
    res = vadd(res, vscale_series(t, linear(a3n, a3k)))
    inv = linear(c0, c1).inverse()
    return vscale_series(res, inv)


# ---------------- cores -----------------
# vertex labels of the uncut cores; edge lists in fixed order
CORES = {
    'K4': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    'W4': [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (1, 4)],   # hub 0, rim 1234
    'prism': [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)],
    'K33': [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)],
}


def cut_period(edges, powers, q, cut=None, order=None):
    """projective period  int prod t^(a-1) / U^q  of the vacuum graph via the cut of edge `cut`.
    Returns the finite-part polynomial dict."""
    from math import factorial, prod
    import laurent
    from laurent import finite_part
    from fractions import Fraction
    if order is not None:
        laurent.ORDER = order
    set_q(q)
    if cut is None:
        cut = 0
    u0, v0 = edges[cut]
    verts = sorted({x for e in edges for x in e})
    mp = {u0: 0, v0: 1}
    nxt = 2
    for w in verts:
        if w not in mp:
            mp[w] = nxt
            nxt += 1
    lines = [(mp[a], mp[b], Fraction(p), 0) for i, ((a, b), p) in enumerate(zip(edges, powers)) if i != cut]
    val = ev(lines)
    qf = Fraction(q)
    if qf.denominator == 1 and all(Fraction(x).denominator == 1 for x in powers):
        pref = Fraction(prod(factorial(int(p) - 1) for p in powers), factorial(int(qf) - 1))
        fp = finite_part(val)
        return {k: v * pref for k, v in fp.items()}
    import halfq
    pref, sk = Fraction(1), 0
    for x in powers:
        r, k = halfq.gamma(x)
        pref *= r
        sk += k
    r, k = halfq.gamma(qf)
    pref /= r
    sk -= k
    fp = finite_part(val)
    out = {}
    for key0, v in fp.items():
        d = dict(key0)
        d['sqrtpi'] = d.get('sqrtpi', 0) + sk
        key = tuple(sorted((a, b) for a, b in d.items() if b))
        out[key] = out.get(key, 0) + v * pref
    return {k: v for k, v in out.items() if v}
