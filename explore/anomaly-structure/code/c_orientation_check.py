"""Consistency of the orientation rule with the engine kernel signs, and the Wick sign of the
0-dimensional model.

 (1) For theta(1,3,3), theta(2,2,3) at (H,T) = (2,1) and (3,0): random relabellings of the vertices
     (any vertex may become the root) and random edge reversals; the engine kernel K' (symbolic
     momenta at the external vertices) and the colour weight W' of the relabelled graph give the same
     alternated cochain alt(K' W') as the original (sl2 and gl2 components, exact).
 (2) Wick sign of the model with odd generators theta^a_v, odd vertex measures mu_v, vertices
     mu_v f theta theta theta (even), contraction theta_{e0} theta_{e1} -> kappa p_{e0 e1}:
     sign of the permutation  mu_0 B_0 mu_1 B_1 ... -> (pairs)(mu_0 ... mu_{m-1})(externals)
     equals (-1)^{3 m (m-1)/2} times the orientation-rule sign, for random labellings.
"""
import sys, random
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import theta, Graph
import lie, cochain as cc

rng = random.Random(3)


def relabel_graph(G, ext, perm, flips):
    """new label of vertex v is perm[v]; edges flipped where flips[e]."""
    edges = []
    for e, (a, b) in enumerate(G.edges):
        a2, b2 = perm[a], perm[b]
        edges.append((b2, a2) if flips[e] else (a2, b2))
    return Graph(G.nv, edges, "relabelled"), [perm[v] for v in ext]


def cochain_of(G, ext, H, T, g, mom):
    P = mom.p
    lam = {v: (list(P[ext.index(v)]) if v in ext else [0] * H) for v in range(G.nv - 1)}
    K = sp.expand(htk.kernel_exact(G, lam, T=T) * (60 / sp.pi if T else 6))
    W, sc = cc.colour_weight(g, G.nv, G.edges, ext)
    return cc.alternate(mom, K, W, 4)


def same(c1, c2):
    keys = set(c1) | set(c2)
    return all(np.array_equal(c1.get(k, 0) * 1, c2.get(k, 0) * 1) for k in keys)


from htk import polygon
cases = []
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, interior = theta(*pqr)
    cases.append((2, 1, "theta%d%d%d" % pqr, G, [v for P_ in interior for v in P_]))
cases.append((3, 0, "box", polygon(4), [0, 1, 2, 3]))
for H, T, name, G, ext in cases:
    mom = cc.Momenta(6, H)
    if True:
        for g in [lie.sl(2), lie.gl(2)]:
            ref = cochain_of(G, ext, H, T, g, mom)
            oks = []
            for trial in range(4):
                perm = list(range(G.nv))
                rng.shuffle(perm)
                flips = [rng.random() < 0.5 for _ in G.edges]
                G2, ext2 = relabel_graph(G, ext, perm, flips)
                oks.append(same(ref, cochain_of(G2, ext2, H, T, g, mom)))
            print("(H,T)=(%d,%d) %s %s: alt(K W) invariant under 4 random relabellings/reversals: %s; reference zero: %s"
                  % (H, T, name, g.name, oks, cc.is_zero(ref)), flush=True)


def wick_sign(nv, edges, ext):
    L_rule = [('x', v) for v in ext]
    for e in range(len(edges)):
        L_rule += [('e', e, 0), ('e', e, 1)]
    local = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v],
                     key=lambda h: (h[1], h[2]))
        local.append(hs)
    Lp = [h for hs in local for h in hs]
    pos = {h: i for i, h in enumerate(L_rule)}
    s_rule = cc.perm_sign([pos[h] for h in Lp])
    toy = []
    for v in range(nv):
        toy += [('m', v)] + local[v]
    target = []
    for e in range(len(edges)):
        target += [('e', e, 0), ('e', e, 1)]
    target += [('m', v) for v in range(nv)] + [('x', v) for v in ext]
    tpos = {h: i for i, h in enumerate(target)}
    s_toy = cc.perm_sign([tpos[h] for h in toy])
    return s_toy, s_rule


res = []
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, interior = theta(*pqr)
    ext = [v for P_ in interior for v in P_]
    m = G.nv
    for trial in range(20):
        perm = list(range(G.nv)); rng.shuffle(perm)
        flips = [rng.random() < 0.5 for _ in G.edges]
        G2, ext2 = relabel_graph(G, ext, perm, flips)
        ext2 = list(ext2); rng.shuffle(ext2)
        s_toy, s_rule = wick_sign(G2.nv, G2.edges, ext2)
        res.append(s_toy == (-1) ** (3 * m * (m - 1) // 2) * s_rule)
print("(2) Wick sign = (-1)^{3m(m-1)/2} x orientation-rule sign (m = 6), 40 random labellings:", all(res))
