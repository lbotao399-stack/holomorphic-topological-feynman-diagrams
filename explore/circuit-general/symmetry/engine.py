"""Symmetry engine for weighted core periods.

    I_q(G; a) = int_{RP_{>0}} prod_e t_e^(a_e-1) U_G^-q Omega ,   b_e = q - a_e   (d = 2q)

Completed graph Gc (one extra vertex 'inf', weights b, every vertex balanced: sum_{e at v} b_e = 2q).
Conformal period (position space, measure d^dx/pi^(d/2), three vertices at 0, 1, inf):

    Pi(Gc) = Gamma(q) I_q(Gc - c; a) / prod_{e in Gc - c} Gamma(b_e)     (any vertex c)

Moves on completed graphs (Pi(G) = factor * Pi(G')):
  M1  star -> triangle at a cubic vertex (always unique)       factor prod Gamma(q-b_i)/Gamma(b_i)
  M2  triangle (sum b = q) -> star                               inverse of M1
  FT  planar decompletion H = Gc - c  ->  completion of H*       factor prod_{e in H} Gamma(q-b_e)/Gamma(b_e)
  TW  twist at a 4-vertex cut (equal G2-weights in each pair)    factor 1
  PR  product at a 3-vertex cut                                  Pi(G) = Pi(G1) Pi(G2)
Terminal: a decompletion Gc - c evaluable by cores.period (SP reduction to banana, W_n, prism, K33;
prism/K33 only for q <= QTAB).
Values: dict  monomial -> Fraction,  monomial = sorted tuple of zeta labels ('z3', 'z5', ...).
"""
import heapq
import os
import sys
from fractions import Fraction as F
from itertools import combinations
from math import factorial

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
CORES = os.path.normpath(os.path.join(HERE, '..', 'cores'))
if CORES not in sys.path:
    sys.path.insert(0, CORES)
import cores  # noqa: E402
import gutil  # noqa: E402

QTAB = 7


class NotTabulated(Exception):
    pass


_orig_table = cores.table


def _guarded_table(name, q, verbose=False):
    if F(q).denominator != 1 or F(q) > QTAB:
        raise NotTabulated('%s q=%s' % (name, q))
    return _orig_table(name, q, verbose=False)


cores.table = _guarded_table


# ---------------------------------------------------------------- values
def gam(x):
    x = int(x)
    if x < 1:
        raise ZeroDivisionError('Gamma(%d)' % x)
    return F(factorial(x - 1))


def from_cores(d):
    return {(() if k == '1' else (k,)): F(v) for k, v in d.items() if v}


def vmul(a, b):
    out = {}
    for m1, x in a.items():
        for m2, y in b.items():
            m = tuple(sorted(m1 + m2))
            out[m] = out.get(m, 0) + x * y
    return {m: v for m, v in out.items() if v}


def vscale(a, s):
    return {m: v * s for m, v in a.items() if v * s}


def vstr(a):
    if not a:
        return '0'
    parts = []
    for m, v in sorted(a.items(), key=lambda t: (len(t[0]), t[0])):
        parts.append('%s%s' % (v, ''.join('*' + x for x in m)))
    return ' + '.join(parts)


ONE = {(): F(1)}


# ---------------------------------------------------------------- weighted graphs
def addw(E, i, j, w):
    if i == j:
        raise ValueError('self-loop')
    k = (i, j) if i < j else (j, i)
    x = E.get(k, 0) + w
    if x:
        E[k] = x
    else:
        E.pop(k, None)


def complete(edges, alpha, q):
    """completed graph of (edges, alpha); alpha >= 1 assumed.  Returns (V, E, inf)"""
    verts = sorted({x for e in edges for x in e})
    lab = {v: i for i, v in enumerate(verts)}
    E = {}
    for (u, v), a in zip(edges, alpha):
        addw(E, lab[u], lab[v], q - a)
    V = len(verts)
    deg = [0] * V
    for (i, j), w in E.items():
        deg[i] += w
        deg[j] += w
    for v in range(V):
        if 2 * q - deg[v]:
            addw(E, v, V, 2 * q - deg[v])
    return V + 1, E, V


def balanced(V, E, q):
    deg = [0] * V
    for (i, j), w in E.items():
        deg[i] += w
        deg[j] += w
    return all(x == 2 * q for x in deg)


def nbrs(V, E):
    adj = [dict() for _ in range(V)]
    for (i, j), w in E.items():
        adj[i][j] = w
        adj[j][i] = w
    return adj


def decomplete(V, E, c, q):
    """H = Gc - c as (edges, alpha, beta) with every b in [1, q-1]; None otherwise"""
    edges, al, be = [], [], []
    verts = [v for v in range(V) if v != c]
    lab = {v: i for i, v in enumerate(verts)}
    for (i, j), w in E.items():
        if c in (i, j):
            continue
        if not 1 <= w <= q - 1:
            return None
        edges.append((lab[i], lab[j]))
        al.append(q - w)
        be.append(w)
    used = {x for e in edges for x in e}
    if len(used) != len(verts):
        return None
    if not connected(len(verts), edges):
        return None
    return edges, al, be


def connected(n, edges):
    if n == 0:
        return False
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen = {0}
    st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen) == n


def cut_convergent(V, E, q):
    """cut(W) > 2q for every vertex set W, 2 <= |W| <= V-2"""
    verts = range(V)
    for k in range(2, V // 2 + 1):
        for W in combinations(verts, k):
            Ws = set(W)
            cut = sum(w for (i, j), w in E.items() if (i in Ws) != (j in Ws))
            if cut <= 2 * q:
                return False
    return True


# ---------------------------------------------------------------- canonical form
def _refine(V, adj, col):
    while True:
        sig = [(col[v], tuple(sorted((col[u], w) for u, w in adj[v].items()))) for v in range(V)]
        keys = sorted(set(sig))
        rank = {s: i for i, s in enumerate(keys)}
        new = [rank[s] for s in sig]
        if len(keys) == len(set(col)):
            return new
        col = new


def canon(V, E):
    """canonical key (V, ((i, j, w), ...)) of a weighted simple graph"""
    adj = nbrs(V, E)
    init = [tuple(sorted(adj[v].values())) for v in range(V)]
    keys = sorted(set(init))
    col = [keys.index(c) for c in init]
    best = [None]

    def rec(col):
        col = _refine(V, adj, col)
        if len(set(col)) == V:
            cert = tuple(sorted((min(col[i], col[j]), max(col[i], col[j]), w) for (i, j), w in E.items()))
            if best[0] is None or cert < best[0]:
                best[0] = cert
            return
        cells = {}
        for v in range(V):
            cells.setdefault(col[v], []).append(v)
        target = min(c for c, vs in cells.items() if len(vs) > 1)
        for v in cells[target]:
            nc = [2 * c + 1 for c in col]
            nc[v] = 2 * col[v]
            rec(nc)
    rec(col)
    return (V, best[0])


def unkey(key):
    V, cert = key
    return V, {(i, j): w for i, j, w in cert}


# ---------------------------------------------------------------- terminal evaluation
def sp_skeleton(edges):
    """unweighted series-parallel reduction of a multigraph; returns simple edge list"""
    E = [tuple(sorted(e)) for e in edges]
    while True:
        E = sorted(set(E))
        deg = {}
        for u, v in E:
            deg[u] = deg.get(u, 0) + 1
            deg[v] = deg.get(v, 0) + 1
        d2 = [v for v in deg if deg[v] == 2]
        if not d2 or len(deg) <= 2:
            return E
        v = d2[0]
        inc = [e for e in E if v in e]
        a = inc[0][0] if inc[0][1] == v else inc[0][1]
        b = inc[1][0] if inc[1][1] == v else inc[1][1]
        E = [e for e in E if v not in e]
        if a != b:
            E.append(tuple(sorted((a, b))))


def skeleton_type(edges):
    E = sp_skeleton(edges)
    verts = sorted({x for e in E for x in e})
    n = len(verts)
    if n <= 2:
        return 'SP'
    deg = {v: 0 for v in verts}
    for u, v in E:
        deg[u] += 1
        deg[v] += 1
    if len(E) == 2 * (n - 1) and n >= 4:
        hubs = [v for v in verts if deg[v] == n - 1]
        rest = [v for v in verts if deg[v] == 3]
        if hubs and len(rest) + 1 == n or (n == 4 and len(E) == 6):
            return 'W%d' % (n - 1)
    if n == 6 and len(E) == 9 and all(deg[v] == 3 for v in verts):
        G = nx.Graph(E)
        return 'prism' if sum(nx.triangles(G).values()) else 'K33'
    return None


STATS = {'eval': 0, 'eval_ok': 0, 'nonconv_terminal': 0}


def eval_I(edges, alpha, q):
    """I_q(edges; alpha) if tabulated, else None"""
    t = skeleton_type(edges)
    if t is None:
        return None
    if t in ('prism', 'K33') and q > QTAB:
        return None
    STATS['eval'] += 1
    try:
        val = cores.period(list(edges), list(alpha), q)
    except (ValueError, NotTabulated, ArithmeticError, ZeroDivisionError, AssertionError, KeyError):
        return None
    STATS['eval_ok'] += 1
    return from_cores(val), t


def pi_from_H(H, q):
    edges, al, be = H
    r = eval_I(edges, al, q)
    if r is None:
        return None
    val, t = r
    norm = gam(q)
    for b in be:
        norm /= gam(b)
    return vscale(val, norm), t


def terminal(V, E, q, skip=()):
    if V <= 3:
        return ONE, 'triangle'
    for c in range(V):
        if c in skip:
            continue
        H = decomplete(V, E, c, q)
        if H is None:
            continue
        r = pi_from_H(H, q)
        if r is not None:
            if not gutil.convergent(tuple(H[0]), H[1], q):
                STATS['nonconv_terminal'] += 1
                continue
            return r[0], 'c=%d:%s' % (c, r[1])
    return None


# ---------------------------------------------------------------- numerator (two-point) terminal
def poly_mul(a, b):
    out = {}
    for m1, x in a.items():
        for m2, y in b.items():
            m = tuple(u + v for u, v in zip(m1, m2))
            out[m] = out.get(m, 0) + x * y
    return out


def forest_poly(n, edges, a, b):
    """Psi^{ab}: sum over spanning 2-forests separating a, b of prod t_e (exponent tuples)"""
    ident = [((a if u == b else u), (a if v == b else v)) for u, v in edges]
    keep = [i for i, (u, v) in enumerate(ident) if u != v]
    sub = tuple(ident[i] for i in keep)
    out = {}
    for T in gutil.spanning_trees(sub):
        m = [0] * len(edges)
        for j in T:
            m[keep[j]] = 1
        out[tuple(m)] = out.get(tuple(m), 0) + 1
    return out


def pi_two_point(V, E, c, ab, q, kmax=3, mmax=20000):
    """Pi(G) with c at infinity and the edge ab between the points 0 and 1, b_ab = -k <= 0:
    Pi = Gamma(q+k)/prod_{H'} Gamma(b) sum_m c_m I_{q+k}(H'; q+k-b-m),  (Psi^{ab}_{H'})^k = sum c_m t^m,
    H' = G - c - ab;  every other edge of H' must have b in [1, q-1]."""
    a, b = ab
    k = -E.get((min(a, b), max(a, b)), 0)
    if k < 0 or k > kmax:
        return None
    verts = [v for v in range(V) if v != c]
    lab = {v: i for i, v in enumerate(verts)}
    edges, be = [], []
    for (i, j), w in E.items():
        if c in (i, j) or (i, j) == (min(a, b), max(a, b)):
            continue
        if not 1 <= w <= q - 1:
            return None
        edges.append((lab[i], lab[j]))
        be.append(w)
    n = len(verts)
    if not connected(n, edges):
        return None
    if skeleton_type(edges) is None:
        return None
    qq = q + k
    P = {tuple([0] * len(edges)): 1}
    if k:
        Fp = forest_poly(n, edges, lab[a], lab[b])
        for _ in range(k):
            P = poly_mul(P, Fp)
            if len(P) > mmax:
                return None
    tot = {}
    for m, cm in P.items():
        al = [qq - x - y for x, y in zip(be, m)]
        if min(al) < 1:
            return None
        if not gutil.convergent(tuple(edges), al, qq):
            return None
        r = eval_I(edges, al, qq)
        if r is None:
            return None
        for mon, v in r[0].items():
            tot[mon] = tot.get(mon, 0) + cm * v
    norm = gam(qq)
    for x in be:
        norm /= gam(x)
    return vscale({mm: v for mm, v in tot.items() if v}, norm), 'c=%d,ab=%d-%d,k=%d:%s' % (c, a, b, k, skeleton_type(edges))


def terminal_numer(V, E, q, kmax=3):
    """decompletion at c with all negative edges among {edges at c} + one edge ab"""
    for c in range(V):
        neg = [(i, j) for (i, j), w in E.items() if w < 0 and c not in (i, j)]
        if len(neg) > 1:
            continue
        cands = neg if neg else []
        for ab in cands:
            r = pi_two_point(V, E, c, ab, q, kmax)
            if r is not None:
                return r
    return None


# ---------------------------------------------------------------- moves
def relabel_drop(V, E, v):
    lab = {x: (x if x < v else x - 1) for x in range(V) if x != v}
    E2 = {}
    for (i, j), w in E.items():
        if v in (i, j):
            continue
        addw(E2, lab[i], lab[j], w)
    return V - 1, E2


def move_M1(V, E, v, q):
    adj = nbrs(V, E)
    nb = list(adj[v].items())
    if len(nb) != 3:
        return None
    if any(not 1 <= w <= q - 1 for _, w in nb):
        return None
    (a, ba), (b, bb), (c, bc) = nb
    fac = F(1)
    for w in (ba, bb, bc):
        fac *= gam(q - w) / gam(w)
    E2 = dict(E)
    for x in (a, b, c):
        E2.pop((min(v, x), max(v, x)))
    addw(E2, b, c, q - ba)
    addw(E2, a, c, q - bb)
    addw(E2, a, b, q - bc)
    V2, E3 = relabel_drop(V, E2, v)
    return V2, E3, fac


def move_M2(V, E, tri, gam3, q):
    """triangle (a,b,c) with parts g_ab, g_bc, g_ca (sum q, each >= 1, <= edge weight) -> star"""
    a, b, c = tri
    gab, gbc, gca = gam3
    E2 = dict(E)
    addw(E2, a, b, -gab)
    addw(E2, b, c, -gbc)
    addw(E2, c, a, -gca)
    s = V
    addw(E2, s, a, q - gbc)
    addw(E2, s, b, q - gca)
    addw(E2, s, c, q - gab)
    fac = F(1)
    for g in (gab, gbc, gca):
        fac *= gam(q - g) / gam(g)
    return V + 1, E2, fac


_PD = {}


def planar_dual(n, edges):
    """dual of a 2-connected planar simple graph: list of (fa, fb) per edge, number of faces;
    None if not planar or not 2-connected (cached on the labelled edge list)"""
    key = (n, tuple(edges))
    if key in _PD:
        return _PD[key]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    res = None
    if nx.is_biconnected(G):
        ok, emb = nx.check_planarity(G)
        if ok:
            faces = {}
            fid = 0
            visited = set()
            for u, v in emb.edges():
                if (u, v) in visited:
                    continue
                face = emb.traverse_face(u, v, mark_half_edges=visited)
                for i in range(len(face)):
                    faces[(face[i], face[(i + 1) % len(face)])] = fid
                fid += 1
            res = ([(faces[(u, v)], faces[(v, u)]) for u, v in edges], fid)
    if len(_PD) > 500000:
        _PD.clear()
    _PD[key] = res
    return res


def move_FT(V, E, c, q):
    H = decomplete(V, E, c, q)
    if H is None:
        return None
    edges, al, be = H
    order = sorted(range(len(edges)), key=lambda i: edges[i])
    edges = [edges[i] for i in order]
    be = [be[i] for i in order]
    pd = planar_dual(V - 1, edges)
    if pd is None:
        return None
    dual, nf = pd
    if any(a == b for a, b in dual):
        return None
    fac = F(1)
    for b in be:
        fac *= gam(q - b) / gam(b)
    # dual weights: a* = b, i.e. alpha_dual = beta
    V2, E2, _ = complete(dual, be, q)
    return V2, E2, fac


def masks(V, E):
    m = [0] * V
    for (i, j) in E:
        m[i] |= 1 << j
        m[j] |= 1 << i
    return m


def components_without(V, adj, S):
    """components of G - S; adj = list of neighbour bitmasks"""
    rest = ((1 << V) - 1) & ~sum(1 << s for s in S)
    comps = []
    while rest:
        low = rest & -rest
        comp = low
        frontier = low
        while frontier:
            b = frontier & -frontier
            frontier ^= b
            nb = adj[b.bit_length() - 1] & rest & ~comp
            comp |= nb
            frontier |= nb
        rest &= ~comp
        comps.append([v for v in range(V) if comp >> v & 1])
    return comps


def product_split(V, E, q):
    """first 3-vertex cut with both sides of >= 2 vertices (one-vertex sides = M1):
    returns ([(V1, E1), (V2, E2)], S) (completed pieces) or None"""
    adj = masks(V, E)
    for S in combinations(range(V), 3):
        comps = components_without(V, adj, S)
        if len(comps) < 2:
            continue
        if len(comps[0]) < 2 or sum(len(cp) for cp in comps[1:]) < 2:
            continue
        C1 = set(comps[0])
        C2 = set(x for cp in comps[1:] for x in cp)
        pieces = []
        for C in (C1, C2):
            verts = sorted(C) + list(S)
            lab = {x: i for i, x in enumerate(verts)}
            E2 = {}
            w = {s: 0 for s in S}
            for (i, j), x in E.items():
                if i in C or j in C:
                    addw(E2, lab[i], lab[j], x)
                    for s in S:
                        if s in (i, j):
                            w[s] += x
            a, b, c = S
            t_ab = q - F(w[a] + w[b] - w[c], 2)
            t_bc = q - F(w[b] + w[c] - w[a], 2)
            t_ca = q - F(w[c] + w[a] - w[b], 2)
            for (x, y), t in (((a, b), t_ab), ((b, c), t_bc), ((c, a), t_ca)):
                if t.denominator != 1:
                    raise ArithmeticError('half-integer triangle weight')
                if t:
                    addw(E2, lab[x], lab[y], int(t))
            pieces.append((len(verts), E2))
        return pieces, S
    return None


def twists(V, E, q):
    """all twisted graphs at 4-vertex cuts"""
    adj = masks(V, E)
    out = []
    for S in combinations(range(V), 4):
        comps = components_without(V, adj, S)
        if len(comps) < 2:
            continue
        for C in comps:
            Cs = set(C)
            w2 = {s: sum(x for (i, j), x in E.items() if s in (i, j) and ((i in Cs) or (j in Cs))) for s in S}
            a, b, c, d = S
            for sig in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
                (p1, p2), (p3, p4) = sig
                if w2[p1] != w2[p2] or w2[p3] != w2[p4]:
                    continue
                mp = {p1: p2, p2: p1, p3: p4, p4: p3}
                E2 = {}
                for (i, j), x in E.items():
                    if i in Cs or j in Cs:
                        i2 = mp.get(i, i)
                        j2 = mp.get(j, j)
                        addw(E2, i2, j2, x)
                    else:
                        addw(E2, i, j, x)
                out.append(((V, E2), (S, sig, tuple(sorted(C)))))
    return out


# ---------------------------------------------------------------- search
MEMO = {}
_ACTIVE = set()
_FLAG = [False]


def successors(V, E, q, use_m2=True, split_m2=False):
    adj = nbrs(V, E)
    out = []
    for v in range(V):
        if len(adj[v]) == 3:
            r = move_M1(V, E, v, q)
            if r is not None:
                out.append((r[0], r[1], r[2], 'M1(%d)' % v))
    for c in range(V):
        r = move_FT(V, E, c, q)
        if r is not None:
            out.append((r[0], r[1], r[2], 'FT(%d)' % c))
    for (G2, info) in twists(V, E, q):
        out.append((G2[0], G2[1], F(1), 'TW%s' % (info[0],)))
    if use_m2:
        for a, b, c in combinations(range(V), 3):
            if b in adj[a] and c in adj[a] and c in adj[b]:
                wab, wbc, wca = adj[a][b], adj[b][c], adj[c][a]
                if min(wab, wbc, wca) < 1:
                    continue
                if wab + wbc + wca == q:
                    r = move_M2(V, E, (a, b, c), (wab, wbc, wca), q)
                    out.append((r[0], r[1], r[2], 'M2(%d,%d,%d)' % (a, b, c)))
                elif split_m2 and wab + wbc + wca > q:
                    for gab in range(1, min(wab, q - 2) + 1):
                        for gbc in range(1, min(wbc, q - 1 - gab) + 1):
                            gca = q - gab - gbc
                            if 1 <= gca <= wca:
                                r = move_M2(V, E, (a, b, c), (gab, gbc, gca), q)
                                out.append((r[0], r[1], r[2], 'M2s(%d,%d,%d;%d,%d,%d)' % (a, b, c, gab, gbc, gca)))
    return out


NUMER = [0]


def solve(V, E, q, budget=200, vmax=None, use_m2=True, split_m2=False, depth=0, all_routes=False):
    """Pi(G) of a balanced completed graph: (value, chain) or None"""
    key0 = canon(V, E)
    mk = (key0, q, budget, vmax, use_m2, split_m2, NUMER[0])
    if mk in MEMO and not all_routes:
        return MEMO[mk]
    if mk in _ACTIVE:
        _FLAG[0] = True
        return None
    _ACTIVE.add(mk)
    saved = _FLAG[0]
    _FLAG[0] = False
    try:
        res = _solve(V, E, q, key0, mk, budget, vmax, use_m2, split_m2, depth, all_routes)
    finally:
        _ACTIVE.discard(mk)
    if not all_routes and (res is not None or not _FLAG[0]):
        MEMO[mk] = res
    _FLAG[0] = saved or _FLAG[0]
    return res


def _solve(V, E, q, key0, mk, budget, vmax, use_m2, split_m2, depth, all_routes):
    if vmax is None:
        vmax = V + (1 if use_m2 else 0)
    found = []
    heap = [(V, 0, 0, key0, F(1), ())]
    seen = {key0}
    cnt = 0
    tie = 0
    while heap and cnt < budget:
        _, d, _, key, fac, chain = heapq.heappop(heap)
        cnt += 1
        Vk, Ek = unkey(key)
        r = terminal(Vk, Ek, q)
        if r is None and NUMER[0]:
            r = terminal_numer(Vk, Ek, q, NUMER[0])
        if r is not None:
            res = (vscale(r[0], fac), chain + ('T:' + r[1],))
            if not all_routes:
                return res
            found.append(res)
            continue
        ps = product_split(Vk, Ek, q)
        if ps is not None and depth < 12:
            pieces, S = ps
            vals = []
            for (Vp, Ep) in pieces:
                if not balanced(Vp, Ep, q):
                    raise ArithmeticError('unbalanced piece')
                rr = solve(Vp, Ep, q, budget, None, use_m2, split_m2, depth + 1)
                if rr is None:
                    break
                vals.append(rr)
            if len(vals) == 2:
                res = (vscale(vmul(vals[0][0], vals[1][0]), fac),
                       chain + ('PR%s[%s | %s]' % (S, ' '.join(vals[0][1]), ' '.join(vals[1][1])),))
                if not all_routes:
                    return res
                found.append(res)
                continue
        for V2, E2, f2, desc in successors(Vk, Ek, q, use_m2, split_m2):
            if V2 > vmax:
                continue
            k2 = canon(V2, E2)
            if k2 in seen:
                continue
            seen.add(k2)
            tie += 1
            heapq.heappush(heap, (V2, d + 1, tie, k2, fac * f2, chain + (desc,)))
    if all_routes:
        return found
    return None


def period(edges, alpha, q, **kw):
    """I_q(G; alpha) through the engine: (value, chain) or None.  alpha_e = 0: contraction."""
    edges = [tuple(e) for e in edges]
    alpha = list(alpha)
    while 0 in alpha:
        i = alpha.index(0)
        u, v = edges[i]
        edges = [((u if a == v else a), (u if b == v else b)) for k, (a, b) in enumerate(edges) if k != i]
        alpha = [x for k, x in enumerate(alpha) if k != i]
    if any(a == b for a, b in edges):
        return None
    if not gutil.convergent(tuple(edges), alpha, q):
        return None
    V, E, inf = complete(edges, alpha, q)
    r = solve(V, E, q, **kw)
    if r is None:
        return None
    norm = F(1) / gam(q)
    for a in alpha:
        norm *= gam(q - a)
    if kw.get('all_routes'):
        return [(vscale(v, norm), ch) for v, ch in r]
    return vscale(r[0], norm), r[1]


def orbit_states(V, E, q, budget=60, vmax=None, use_m2=True):
    """all completed graphs reached from G by the moves (no product split, no terminal stop):
    list of (V', E', f) with Pi(G) = f * Pi(G')"""
    key0 = canon(V, E)
    if vmax is None:
        vmax = V + (1 if use_m2 else 0)
    out = []
    heap = [(V, 0, 0, key0, F(1))]
    seen = {key0}
    tie = 0
    while heap and len(out) < budget:
        _, d, _, key, fac = heapq.heappop(heap)
        Vk, Ek = unkey(key)
        out.append((Vk, Ek, fac))
        for V2, E2, f2, desc in successors(Vk, Ek, q, use_m2):
            if V2 > vmax:
                continue
            k2 = canon(V2, E2)
            if k2 in seen:
                continue
            seen.add(k2)
            tie += 1
            heapq.heappush(heap, (V2, d + 1, tie, k2, fac * f2))
    return out
