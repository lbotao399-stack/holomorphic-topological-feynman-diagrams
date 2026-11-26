"""Additional regulator-free relations for the period systems (all terms convergent):

(Y) triangle -> star at a unique triangle (powers c1+c2+c3 = 2q):
      I_q(G; c) = I_q(G - triangle + star; a),  a_{v x_k} = q - c(opposite edge of x_k)
    (inverse of star-triangle (S), startri.py).
(X) position-space vertex IBP (Fourier transform of every propagator, beta = q - alpha):
      vertex x, neighbour x1 adjacent to every other neighbour x_i of x, e1 = x x1, e_i = x x_i,
      f_i = x1 x_i;  v(a) = prod_e 1/Gamma(q - a_e)  (position integral X(q - a) = const v(a) I(a)):
      (2 a_1 + sum_i a_i - (k-1) q) v(a) I(a)
         = sum_i (q - a_i) [ v(a + e_1 - e_i) I(a + e_1 - e_i) - v(a - e_i + f_i) I(a - e_i + f_i) ]
    used only if every term has 1 <= a_e <= q-1 (a term with a_e = 0 carries a boundary term:
    check_rel5.log, 0-term cases fail when dropped).
(P) planar duality  I_q(G*; a) = I_q(G; q - a o pi)  (self-dual cores: Z5).
"""
from fractions import Fraction as F
from math import factorial

import networkx as nx


def vfac(a, q):
    """v(a) = prod 1/Gamma(q-a) as Fraction; 0 if some a_e = 0 (marker: divergent X); None if some a_e >= q"""
    r = F(1)
    for x in a:
        if x == 0:
            return F(0)
        if x >= q or x < 0:
            return None
        r /= factorial(q - x - 1)
    return r


def triangles(edges):
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    verts = sorted({x for e in edges for x in e})
    out = []
    for i, a in enumerate(verts):
        for j in range(i + 1, len(verts)):
            b = verts[j]
            if frozenset((a, b)) not in idx:
                continue
            for c in verts[j + 1:]:
                if frozenset((a, c)) in idx and frozenset((b, c)) in idx:
                    out.append((a, b, c, idx[frozenset((b, c))], idx[frozenset((a, c))], idx[frozenset((a, b))]))
    return out


def tri_to_star(edges, pw, tri, q):
    a, b, c, ebc, eac, eab = tri
    v = max(x for e in edges for x in e) + 1
    ne = [e for k, e in enumerate(edges) if k not in (ebc, eac, eab)]
    npw = [p for k, p in enumerate(pw) if k not in (ebc, eac, eab)]
    ne += [(v, a), (v, b), (v, c)]
    npw += [q - pw[ebc], q - pw[eac], q - pw[eab]]
    return ne, npw


def x_rules(edges):
    """list of (e1, [(e_i, f_i)]) for the position-space vertex IBP"""
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    verts = sorted({x for e in edges for x in e})
    nbr = {v: [] for v in verts}
    for u, w in edges:
        nbr[u].append(w)
        nbr[w].append(u)
    out = []
    for x in verts:
        for x1 in nbr[x]:
            others = [y for y in nbr[x] if y != x1]
            if all(frozenset((x1, y)) in idx for y in others):
                out.append((idx[frozenset((x, x1))], [(idx[frozenset((x, y))], idx[frozenset((x1, y))]) for y in others]))
    return out


def x_items(edges, pw, rule, q):
    e1, pairs = rule
    k = len(pairs) + 1
    c = 2 * pw[e1] + sum(pw[ei] for ei, fi in pairs) - (k - 1) * q
    v0 = vfac(pw, q)
    if v0 is None or v0 == 0:
        return None
    items = [(list(edges), list(pw), F(c))]
    for ei, fi in pairs:
        for sgn, shift in ((-1, ((e1, 1), (ei, -1))), (1, ((ei, -1), (fi, 1)))):
            b = list(pw)
            for e, d in shift:
                b[e] += d
            vb = vfac(b, q)
            if vb is None or vb == 0:
                return None          # a term with a_e = 0 (beta_e = q) or a_e = q: IBP boundary term, not used
            items.append((list(edges), b, sgn * (q - pw[ei]) * vb / v0))
    return items


def dual_map(edges):
    """planar dual: (dual edge list in the same order, face count); None if not planar"""
    G = nx.Graph([tuple(e) for e in edges])
    ok, emb = nx.check_planarity(G)
    if not ok:
        return None
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
    return [(faces[(u, v)], faces[(v, u)]) for u, v in edges]


def self_dual_perms(edges):
    """edge permutations sigma with I(G; a) = I(G; q - a o sigma) (G self-dual)"""
    d = dual_map(edges)
    if d is None:
        return []
    G = nx.Graph([tuple(e) for e in edges])
    Gd = nx.Graph(d)
    if Gd.number_of_edges() != len(edges) or not nx.is_isomorphic(G, Gd):
        return []
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    out = set()
    for phi in nx.algorithms.isomorphism.GraphMatcher(Gd, G).isomorphisms_iter():
        # dual edge k (faces) -> edge of G
        perm = tuple(idx[frozenset((phi[a], phi[b]))] for a, b in d)
        out.add(perm)
    return sorted(out)


def add_extra(S, use_Y=True, use_X=True, use_P=True, use_C=True):
    """append (Y), (X), (P), (C) rows for every family orbit of the MultiSystem S"""
    q = S.q
    for key in ('Y', 'X', 'P', 'C'):
        S.counts.setdefault(key, 0)
    for f, d in enumerate(S.fam):
        edges = list(d['edges'])
        tris = triangles(edges) if use_Y else []
        xr = x_rules(edges) if use_X else []
        sdp = self_dual_perms(edges) if use_P else []
        reps = [r for (ff, r) in list(S.reps) if ff == f]
        for r in reps:
            r = list(r)
            for tri in tris:
                if r[tri[3]] + r[tri[4]] + r[tri[5]] == 2 * q:
                    ne, npw = tri_to_star(edges, r, tri, q)
                    S.add_relation([(edges, r, F(1)), (ne, npw, F(-1))], 'Y')
            for rule in xr:
                items = x_items(edges, r, rule, q)
                if items is not None:
                    S.add_relation(items, 'X')
            for perm in sdp:
                b = [q - r[perm[k]] for k in range(len(r))]
                if min(b) >= 1:
                    S.add_relation([(edges, r, F(1)), (edges, b, F(-1))], 'P')
            if use_C:
                for items in completion_items(edges, r, q):
                    S.add_relation(items, 'C')


# ---------------- (C) conformal completion in position space ----------------
def completion_items(edges, pw, q):
    """(C): Xpos(G; q - a) = Xpos(G^ - u; ...) for every vertex u of G, G^ = G + vertex oo with
    position weights w_v = 2q - sum_{e at v} (q - a_e) on (oo, v)  (all vertices conformal).
    In periods: I(G; a) / prod Gamma(q - a) = I(G^-u; a^) / prod Gamma(q - a^),
    a^ = a on edges of G not at u, q - w_v on (oo, v).  Used if every w_v in [0, q-1]
    (w_v = 0: no edge).  Returns a list of item lists."""
    verts = sorted({x for e in edges for x in e})
    inf = max(verts) + 1
    w = {}
    for v in verts:
        s = sum(q - pw[k] for k, e in enumerate(edges) if v in e)
        w[v] = 2 * q - s
        if not 0 <= w[v] <= q - 1:
            return []
    v0 = vfac(pw, q)
    if not v0:
        return []
    out = []
    for u in verts:
        ne, npw = [], []
        for k, e in enumerate(edges):
            if u not in e:
                ne.append(e)
                npw.append(pw[k])
        if w[u] == 0:
            # u had no edge to oo; oo takes u's place
            pass
        for v in verts:
            if v != u and w[v] > 0:
                ne.append((v, inf))
                npw.append(q - w[v])
        vb = vfac(npw, q)
        if not vb:
            continue
        out.append([(list(edges), list(pw), F(1)), (ne, npw, -vb / v0)])
    return out


# ---------------- (L) parametric annihilators (syz5.py) ----------------
_SYZ = {}


def minimal_syzygies(E, kmax):
    """minimal generators of the syzygy module up to degree kmax (cached in syz_cache/)"""
    import os
    import pickle
    import syz5
    key = (tuple(E), kmax)
    if key in _SYZ:
        return _SYZ[key]
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'syz_cache')
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, '%s_k%d.pkl' % ('_'.join('%d-%d' % e for e in E), kmax))
    if os.path.exists(path):
        g = pickle.load(open(path, 'rb'))
    else:
        g = syz5.minimal_generators(E, kmax)
        pickle.dump(g, open(path, 'wb'))
    _SYZ[key] = g
    return g


def add_syzygy_rows(S, degrees=(1, 2), max_bases=None):
    import gutil
    import syz5
    q = S.q
    S.counts.setdefault('L', 0)
    for f, d in enumerate(S.fam):
        E = list(d['edges'])
        N = len(E)
        n = d['n']
        gens = minimal_syzygies(E, max(degrees))
        for k in degrees:
            sy = gens.get(k, [])
            for beta in gutil.compositions(n * q - k + 1, N, 1, q):
                for s in sy:
                    items = syz5.relation(E, s, list(beta), q)
                    if items is None or any(min(g) < 1 for e, g, c in items):
                        continue
                    S.add_relation(items, 'L')
