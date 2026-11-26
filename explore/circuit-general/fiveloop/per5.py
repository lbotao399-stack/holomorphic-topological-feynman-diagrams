"""Exact projective periods I_q(G; a) = int prod t^(a-1) / U^q of graphs whose
series-parallel reduction ends in a supported core, including six-loop K5.

Reduction (msys.MultiSystem.classify, imported unchanged from ../cores): zero powers
(contraction), series, parallel (G_2q), banana, wheels W_n (wheel.py, any n, any integer q).
Hook `lower` (this file) for the remaining 2-connected simple graphs:
  2-vertex cut {u,v}, sides A, B:  I(G) = Gamma(q) I(A+e_uv(s_B)) I(B+e_uv(s_A)) / (Gamma(s_A) Gamma(s_B)),
        s_X = sum_X a - q h(X)                                   (2-sum factorization)
  K5-e (triangular bipyramid) = planar dual of the prism:  I_q(K5-e; a) = I_q(prism; q - a o pi)
  K5: read-only q=3..8 tables via ../period-dispatch/k5_bridge.py (six-loop core)
  prism, K33: tables (fiveloop/tables, assembly/tables, cores/tables; built with msys if absent)
  five-loop cores Z5, K33+e, V7_*, V8_*: tables built here (build5.py) when present.
Values: dict label -> Fraction, labels '1', 'z3', 'z5', 'z7', products 'z3*z3', ...
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from itertools import permutations
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
CORES = os.path.abspath(os.path.join(HERE, '..', 'cores'))
if CORES not in sys.path:
    sys.path.insert(0, CORES)

import gutil   # noqa: E402
import msys5 as msys    # noqa: E402  (copy of cores/msys.py + planar-duality hook)

TABLES_OUT = os.path.join(HERE, 'tables')
TABLES_IN = [TABLES_OUT, os.path.abspath(os.path.join(HERE, '..', 'assembly', 'tables')),
             os.path.join(CORES, 'tables')]
PRISM = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))
K33 = ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5))
FOUR = {'prism': PRISM, 'K33': K33}

# five-loop 3-connected cores (skel5.py, skel5.json) other than W5 and K5-e
_sk = json.load(open(os.path.join(HERE, 'skel5.json')))
SKEL = {r['name']: tuple(tuple(e) for e in r['edges']) for r in _sk}
FIVE_TABLED = ['Z5', 'K33+e', 'V7_0', 'V7_2', 'V7_3', 'V8_1', 'V8_2', 'V8_3', 'V8_4']

STRICT = [False]          # raise on unknown core instead of returning None


class UnknownCore(Exception):
    pass


# ---------------- labels ----------------
def lab_mul(a, b):
    if a == '1':
        return b
    if b == '1':
        return a
    parts = sorted(a.split('*') + b.split('*'), key=lambda s: int(s[1:]))
    return '*'.join(parts)


def vmul(x, y):
    out = {}
    for a, u in x.items():
        for b, v in y.items():
            k = lab_mul(a, b)
            out[k] = out.get(k, 0) + u * v
    return {k: v for k, v in out.items() if v}


def vscale(x, s):
    return {k: v * s for k, v in x.items() if v * s}


def vadd(x, y, s=1):
    out = dict(x)
    for k, v in y.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


# ---------------- tables ----------------
_tables = {}
HIDE = set(x for x in os.environ.get('PER5_HIDE', '').split(',') if x)   # tables ignored (validation runs)


def _load(name, q):
    for d in TABLES_IN:
        path = os.path.join(d, '%s_q%d.json' % (name, q))
        if os.path.exists(path):
            raw = json.load(open(path))
            return {tuple(r): {k: F(v) for k, v in val.items()} for r, val in raw['values']}, raw
    return None, None


def table(name, q, build=True, verbose=False):
    key = (name, q)
    if name in HIDE:
        return None
    if key in _tables:
        return _tables[key]
    tab, raw = _load(name, q)
    if tab is None:
        if name in FOUR and build:
            t = time.time()
            S = msys.MultiSystem([FOUR[name]], q)
            S.build(verbose=verbose)
            solved, piv = S.solve(verbose=verbose)
            if len(solved) != len(S.reps):
                raise ArithmeticError('%s q=%d not closed' % (name, q))
            tab = {S.reps[j][1]: v for j, v in solved.items()}
            os.makedirs(TABLES_OUT, exist_ok=True)
            json.dump({'core': name, 'edges': FOUR[name], 'q': str(q), 'seconds': round(time.time() - t, 1),
                       'values': [[list(r), {k: str(x) for k, x in v.items()}] for r, v in sorted(tab.items())]},
                      open(os.path.join(TABLES_OUT, '%s_q%d.json' % (name, q)), 'w'))
        else:
            _tables[key] = None
            return None
    else:
        if raw.get('edges') is not None and name not in FOUR:
            assert tuple(tuple(e) for e in raw['edges']) == SKEL[name], name
    _tables[key] = tab
    return tab


def save_table(name, q, tab, seconds, extra=None):
    os.makedirs(TABLES_OUT, exist_ok=True)
    d = {'core': name, 'edges': SKEL.get(name, FOUR.get(name)), 'q': str(q), 'seconds': round(seconds, 1),
         'values': [[list(r), {k: str(x) for k, x in v.items()}] for r, v in sorted(tab.items())]}
    if extra:
        d.update(extra)
    json.dump(d, open(os.path.join(TABLES_OUT, '%s_q%d.json' % (name, q)), 'w'))
    _tables[(name, q)] = tab


# ---------------- recognizers ----------------
_recog = {}


def recognizer(q):
    """MultiSystem used only to map a graph isomorphic to a tabulated core onto its orbit rep"""
    if q not in _recog:
        fams = [PRISM, K33] + [SKEL[nm] for nm in FIVE_TABLED]
        S = msys.MultiSystem(fams, q)
        S.lower = lower
        _recog[q] = S
    return _recog[q]


FAM_NAMES = ['prism', 'K33'] + FIVE_TABLED


def family_value(edges, pw, q):
    """value of a graph isomorphic to prism/K33/tabulated five-loop core, else None.
    Raises UnknownCore in STRICT mode if the core is recognized but not tabulated."""
    S = recognizer(q)
    sk = S.skey_cache.get(tuple(sorted(edges)))
    if sk is None:
        vs = sorted({x for e in edges for x in e})
        dg = {v: sum(1 for e in edges if v in e) for v in vs}
        if (len(vs), len(edges), tuple(sorted(dg.values()))) in S.fam_sigs:
            sk = msys.struct_key(edges)
        else:
            sk = ('nofamily',)
        S.skey_cache[tuple(sorted(edges))] = sk
    f = S.fam_by_key.get(sk)
    if f is None:
        return None
    name = FAM_NAMES[f]
    tab = table(name, q)
    if tab is None:
        if STRICT[0]:
            raise UnknownCore('%s at q=%d not tabulated' % (name, q))
        return None
    mapped = S.map_to_family(list(edges), list(pw), f)
    if mapped is None:
        return None
    rep = gutil.orbit_rep(tuple(mapped), S.fam[f]['autos'])
    val = tab.get(rep)
    if val is None:
        return None          # divergent (not in table)
    return dict(val)


# ---------------- structural moves ----------------
def two_cut(edges):
    """a 2-vertex cut (u, v, side_edge_indices) of a 2-connected simple graph, or None"""
    verts = sorted({x for e in edges for x in e})
    if len(verts) < 5:
        return None
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            u, v = verts[i], verts[j]
            rest = [w for w in verts if w not in (u, v)]
            adj = {w: [] for w in rest}
            for a, b in edges:
                if a in adj and b in adj:
                    adj[a].append(b)
                    adj[b].append(a)
            seen = {rest[0]}
            st = [rest[0]]
            while st:
                x = st.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        st.append(y)
            if len(seen) < len(rest):
                side = [k for k, (a, b) in enumerate(edges) if a in seen or b in seen]
                return u, v, side
    return None


def k5e_dual(edges, pw, q):
    """K5-e periods from the prism: I_q(K5-e; a) = I_q(prism; q - a o pi)"""
    verts = sorted({x for e in edges for x in e})
    if len(verts) != 5 or len(edges) != 9:
        return None
    deg = {v: sum(1 for e in edges if v in e) for v in verts}
    T = [v for v in verts if deg[v] == 3]
    S_ = [v for v in verts if deg[v] == 4]
    if len(T) != 2 or len(S_) != 3:
        return None
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    if frozenset(T) in idx:
        return None
    TA, TB = T
    S0, S1, S2 = S_
    Ss = [S0, S1, S2]
    # prism vertices a_i = i, b_i = 3 + i; PRISM edge order (0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(1,4),(2,5)
    dual = {}
    for i in range(3):
        ai, aj = i, (i + 1) % 3
        dual[frozenset((ai, aj))] = frozenset((TA, Ss[i]))
        dual[frozenset((3 + ai, 3 + aj))] = frozenset((TB, Ss[i]))
        dual[frozenset((i, 3 + i))] = frozenset((Ss[(i - 1) % 3], Ss[i]))
    ppw = []
    for e in PRISM:
        k = idx[dual[frozenset(e)]]
        ppw.append(q - pw[k])
    if min(ppw) < 1:
        return None
    return value(list(PRISM), ppw, q)


def _k5_bridge_value(edges, pw, q):
    bridge_path = os.path.abspath(os.path.join(HERE, '..', 'period-dispatch'))
    if bridge_path not in sys.path:
        sys.path.insert(0, bridge_path)
    from k5_bridge import k5_period, K5NotConvergent, K5Unavailable
    try:
        return k5_period(edges, pw, q)
    except K5NotConvergent:
        return None
    except K5Unavailable as error:
        if STRICT[0]:
            raise UnknownCore(str(error)) from error
        return None


def _is_k5(edges):
    return (len(edges) == 10 and len({v for edge in edges for v in edge}) == 5
            and all(a != b for a, b in edges) and len({frozenset(edge) for edge in edges}) == 10)


def lower(edges, pw, q):
    """hook for msys.classify: value of a reduced (2-connected, simple, min degree 3, no wheel) graph"""
    edges = [tuple(e) for e in edges]
    if q in (3, 4, 5, 6, 7, 8) and _is_k5(edges):
        return _k5_bridge_value(edges, pw, q)
    cut = two_cut(edges)
    if cut is not None:
        u, v, side = cut
        A = [edges[k] for k in side]
        Bi = [k for k in range(len(edges)) if k not in side]
        B = [edges[k] for k in Bi]
        pa = [pw[k] for k in side]
        pb = [pw[k] for k in Bi]
        hA = gutil.cycle_rank(A)
        hB = gutil.cycle_rank(B)
        sA = sum(pa) - q * hA
        sB = sum(pb) - q * hB
        if sA <= 0 or sB <= 0 or sA + sB != q:
            return None
        va = value(A + [(u, v)], pa + [sB], q)
        if va is None:
            return None
        vb = value(B + [(u, v)], pb + [sA], q)
        if vb is None:
            return None
        return vscale(vmul(va, vb), F(factorial(q - 1), factorial(sA - 1) * factorial(sB - 1)))
    verts = {x for e in edges for x in e}
    if len(verts) == 5 and len(edges) == 9:
        return k5e_dual(edges, pw, q)
    v = family_value(edges, pw, q)
    if v is None and USE_Y[0]:
        v = y_fallback(edges, pw, q)
    if v is None and _is_k5(edges) and STRICT[0]:
        raise UnknownCore('K5 at q=%s has no available table or successful reduction' % q)
    return v


USE_Y = [True]


def y_fallback(edges, pw, q):
    """unrecognized core: (Y) triangle -> star at a unique triangle (sum = 2q), loop number - 1
    (rel5.py (Y), check_rel5.log); first evaluable choice."""
    import rel5
    for tri in rel5.triangles(edges):
        if pw[tri[3]] + pw[tri[4]] + pw[tri[5]] == 2 * q:
            ne, npw = rel5.tri_to_star(edges, pw, tri, q)
            if min(npw) < 1:
                continue
            v = value(ne, npw, q)
            if v is not None:
                return v
    return None


# ---------------- public ----------------
_red = {}


def reducer(q):
    if q not in _red:
        S = msys.MultiSystem([], q)
        S.lower = lower
        _red[q] = S
    return _red[q]


def value(edges, pw, q):
    """exact period, or None if divergent / core not available"""
    # Complete K5 already has a proof-backed 15-subgraph convergence test and
    # exact orbit index.  Avoid the generic 1,022-subgraph test for every
    # monomial.  Zero weights still follow the ordinary contraction rules.
    if q in (3, 4, 5, 6, 7, 8) and all(a > 0 for a in pw) and _is_k5(edges):
        return _k5_bridge_value(edges, pw, q)
    S = reducer(q)
    res = S.classify([tuple(e) for e in edges], list(pw), q)
    if res is None:
        return None
    const, unk = res
    assert not unk
    return {k: v for k, v in const.items() if v}


def period(edges, pw, q):
    STRICT[0] = True
    try:
        v = value(edges, pw, q)
    finally:
        STRICT[0] = False
    if v is None:
        raise ValueError('divergent or unavailable: %s %s q=%d' % (edges, pw, q))
    return v


def fmt(v):
    return ' + '.join('%s%s' % (x, '' if k == '1' else '*' + k) for k, x in sorted(v.items(), key=lambda kv: (len(kv[0]), kv[0]))) or '0'
