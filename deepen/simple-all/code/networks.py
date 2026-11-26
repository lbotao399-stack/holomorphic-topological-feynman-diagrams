"""Vacuum networks of f-vertices: algebra-independent identities by oriented graph isomorphism.

A network: list of vertices, each an ordered triple of edge labels (value f(e1,e2,e3) = kappa([e1,e2],e3), totally
antisymmetric), every label on exactly two slots (contracted with kappa^{-1}), times a prefactor sign.
If phi is an isomorphism of the underlying (simple) graphs, value(G1) = prefactor ratio * prod_v sgn(pi_v) * value(G2),
pi_v the permutation of the slots of v induced by phi (an identity of tensor networks, valid for every Lie algebra
with invariant form).  All isomorphisms are enumerated; equal signs for all of them are required.
Networks (slots A1..A4 = legs 0..3; index forms of simple-g/code/core.py, equal to t4.Engine.weights, logs/check_W.log):
  W223 = - f(x,y,z) f(A1,x,u) f(A2,y,w) f(A3,z,s) f(A4,s,r) f(u,w,r)
  W133 = - f(x,y,z) f(A1,y,u) f(A2,u,w) f(A3,z,s) f(A4,s,r) f(x,w,r)
  kk_pi probe: joins A_pi0-A_pi1, A_pi2-A_pi3;  ff_pi probe: f(A_pi0,A_pi1,e) f(e,A_pi2,A_pi3)
  ladders: Om(a,b,c,e) = f(x,a,c) f(x,b,e);  tr(Om^k): closure a_k = a_0, b_k = b_0;  tr(P Om^k): a_k = b_0, b_k = a_0.
Each network is also evaluated directly (dense einsum, exact integers via object arrays for small algebras) and
compared with the engine values (logs/b_data.json, form used = Kac form).
"""
import sys, itertools, json
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso
from fractions import Fraction as Fr
import t4

P3 = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]


def W_net(g):
    if g == '223':
        V = [('x', 'y', 'z'), ('A0', 'x', 'u'), ('A1', 'y', 'w'), ('A2', 'z', 's'), ('A3', 's', 'r'), ('u', 'w', 'r')]
    else:
        V = [('x', 'y', 'z'), ('A0', 'y', 'u'), ('A1', 'u', 'w'), ('A2', 'z', 's'), ('A3', 's', 'r'), ('x', 'w', 'r')]
    return V, -1


def probe(V, sign, kind, pi):
    V = [list(v) for v in V]
    ren = {}
    if kind == 'kk':
        ren['A%d' % pi[1]] = 'A%d' % pi[0]
        ren['A%d' % pi[3]] = 'A%d' % pi[2]
        V = [tuple(ren.get(l, l) for l in v) for v in V]
    else:
        V = [tuple(v) for v in V] + [('A%d' % pi[0], 'A%d' % pi[1], 'eN'), ('eN', 'A%d' % pi[2], 'A%d' % pi[3])]
    return V, sign


def ladder(k, twisted):
    V = []
    for j in range(k):
        V.append(('x%d' % j, 'a%d' % j, 'a%d' % (j + 1)))
        V.append(('x%d' % j, 'b%d' % j, 'b%d' % (j + 1)))
    ren = {'a%d' % k: 'b0', 'b%d' % k: 'a0'} if twisted else {'a%d' % k: 'a0', 'b%d' % k: 'b0'}
    V = [tuple(ren.get(l, l) for l in v) for v in V]
    return V, 1


def graph(V):
    G = nx.MultiGraph()
    G.add_nodes_from(range(len(V)))
    where = {}
    for i, v in enumerate(V):
        for l in v:
            where.setdefault(l, []).append(i)
    for l, ends in where.items():
        assert len(ends) == 2, (l, ends)
        G.add_edge(ends[0], ends[1], label=l)
    return G, where


def iso_signs(N1, N2):
    V1, s1 = N1
    V2, s2 = N2
    G1, w1 = graph(V1)
    G2, w2 = graph(V2)
    S1, S2 = nx.Graph(G1), nx.Graph(G2)
    if S1.number_of_edges() != G1.number_of_edges() or S2.number_of_edges() != G2.number_of_edges():
        return None, "multigraph"
    GM = iso.GraphMatcher(S1, S2)
    signs = set()
    n = 0
    for phi in GM.isomorphisms_iter():
        n += 1
        # edge map: label l (ends a,b) -> label of edge (phi a, phi b) in G2
        lab2 = {frozenset((u, v)): d['label'] for u, v, d in G2.edges(data=True)}
        emap = {l: lab2[frozenset((phi[e[0]], phi[e[1]]))] for l, e in w1.items()}
        sg = s1 * s2
        for i, v in enumerate(V1):
            img = [emap[l] for l in v]
            tgt = list(V2[phi[i]])
            sg *= t4.perm_sign([tgt.index(x) for x in img])
        signs.add(sg)
    return signs, n


def evaluate(V, sign, L):
    """dense evaluation in float64 (independent cross-check; integer result, |value| < 2^50 required and checked):
    the second slot carrying each label is raised with G = s kappa^{-1}; einsum over the 3-tensors."""
    import string
    d = L.dim
    a, b, c, v = L.f_low
    f = np.zeros((d, d, d), dtype=np.float64)
    for x, y, z, w in zip(a, b, c, v):
        f[x, y, z] = int(w)
    G = np.array(L.G.toarray(), dtype=np.float64)
    _, where = graph(V)
    letter = dict(zip(where, string.ascii_letters))
    seen = {}
    ops, subs = [], []
    for vv in V:
        T = f
        for pos, l in enumerate(vv):
            k = seen.get(l, 0); seen[l] = k + 1
            if k == 1:
                T = np.moveaxis(np.tensordot(T, G, axes=([pos], [0])), -1, pos)
        ops.append(T); subs.append(''.join(letter[l] for l in vv))
    # pairwise contraction: repeatedly contract the pair whose result has the fewest free indices
    items = [(T, list(sb)) for T, sb in zip(ops, subs)]
    while len(items) > 1:
        best = None
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                common = set(items[i][1]) & set(items[j][1])
                if not common:
                    continue
                free = (set(items[i][1]) | set(items[j][1])) - common
                if best is None or len(free) < best[0]:
                    best = (len(free), i, j, common)
        _, i, j, common = best
        (A, sa), (B, sb) = items[i], items[j]
        common = sorted(common)
        C = np.tensordot(A, B, axes=([sa.index(x) for x in common], [sb.index(x) for x in common]))
        sc = [x for x in sa if x not in common] + [x for x in sb if x not in common]
        items = [it for k, it in enumerate(items) if k not in (i, j)] + [(C, sc)]
    val = float(items[0][0])
    assert abs(val) < 2.0 ** 50 and abs(val - round(val)) < 1e-6, val
    return Fr(int(round(val)), L.s ** len(where)) * sign


if __name__ == "__main__":
    import lie
    tests = []
    lad = {('tr', k): ladder(k, False) for k in (2, 3, 4)}
    lad.update({('trP', k): ladder(k, True) for k in (2, 3, 4)})
    names = {('tr', 3): 'm3 = tr Om^3', ('tr', 4): 'm4 = tr Om^4', ('trP', 3): 'n3 = tr P Om^3', ('trP', 4): 'n4 = tr P Om^4',
             ('tr', 2): 'm2', ('trP', 2): 'n2'}
    nets = {}
    for g in ('133', '223'):
        for kind in ('kk', 'ff'):
            for i, p in enumerate(P3):
                nets['W%s.%s%d' % (g, kind, i)] = probe(*W_net(g), kind, p)
    out = {}
    for nm, N1 in nets.items():
        res = []
        for key, N2 in lad.items():
            if len(N1[0]) != len(N2[0]):
                continue
            signs, n = iso_signs(N1, N2)
            if signs is None:
                res.append('%s: %s' % (names[key], n)); continue
            if n:
                res.append('= %s * %s (%d isomorphisms, signs %s)' % ('+' if signs == {1} else ('-' if signs == {-1} else '?'), names[key], n, sorted(signs)))
        G, _ = graph(N1[0])
        simple = nx.Graph(G).number_of_edges() == G.number_of_edges()
        print("%-10s vertices %d, simple graph %s, triangle-free %s, bipartite %s : %s" % (
            nm, len(N1[0]), simple, sum(nx.triangles(nx.Graph(G)).values()) == 0, nx.is_bipartite(nx.Graph(G)),
            '; '.join(res) if res else 'no ladder isomorphism'), flush=True)
    # direct evaluation versus engine numbers
    recs = {r['key']: r for r in json.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/logs/b_data.json'))}
    for key in sys.argv[1:]:
        L = lie.build(key)
        r = Fr(recs[key]['r'])
        ok = True
        for nm, (V, s) in nets.items():
            val = evaluate(V, s, L)
            e = -3 if '.kk' in nm else -4
            ok &= (val * r ** e == Fr(recs[key]['W'][nm]))
        for (t, k), (V, s) in lad.items():
            val = evaluate(V, s, L)
            ok &= (val * r ** (-k) == Fr(recs[key]['moments'][('m' if t == 'tr' else 'n') + str(k)]))
        print("%s: direct dense float64 evaluation (integer-valued, rounded) of all 12 probe networks and 6 ladder networks = engine values: %s" % (L.name, ok), flush=True)
