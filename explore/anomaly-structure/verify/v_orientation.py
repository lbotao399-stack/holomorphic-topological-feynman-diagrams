"""Partial check of (c3) (status derived-unchecked): alt(K_Gamma W_Gamma) is unchanged under random
vertex relabellings (root moved), edge reversals and edge-order permutations, (H,T) = (2,1) thetas,
own kernel (pointwise constant c with g = c prod_P w_P^{H+1/2-n_P}, K = c * Dirichlet factor, the
factor being labelling independent) and own implementation of the orientation rule.
The Wick-sign derivation (relative sign between different graphs) is not tested here."""
import sys, random, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import vk, sp_coords as spc, lie_eval as le, ce_eval as ce

rng = random.Random(99)
H = 2


def theta_const(nv, edges, paths, lam):
    """c with q |det J| = c prod_P w_P^{H - n_P} (checked at two random points)."""
    C = vk.positive_cycles(nv, edges, vk.dfs_tree(nv, edges))
    ns = [len(P) for P in paths]
    order = [e for P in paths for e in P]
    vals = []
    for _ in range(2):
        w = [Fr(rng.randint(1, 15)) for _ in range(3)]
        sw = sum(w)
        w = [x / sw for x in w]
        xs = []
        for n_ in ns:
            xx = [Fr(rng.randint(1, 15)) for _ in range(n_)]
            sx = sum(xx)
            xs.append([x / sx for x in xx])
        tl, J = spc.jac(ns, w, xs)
        t = [None] * len(edges)
        for pos, e in enumerate(order):
            t[e] = tl[pos]
        q, U = vk.density(nv, edges, t, lam, H, 1, C)
        assert U == 1 / (w[0] * w[1] * w[2])
        base = Fr(1)
        for P in range(3):
            base *= w[P] ** (H - ns[P])
        vals.append(q * abs(vk.det(J)) / base)
    assert vals[0] == vals[1], "not constant on the fibres"
    return vals[0]


def relabel(nv, edges, paths, ext):
    pi = list(range(nv))
    rng.shuffle(pi)
    eperm = list(range(len(edges)))
    rng.shuffle(eperm)  # new edge k = old edge eperm[k]
    rev = [rng.random() < 0.5 for _ in edges]
    new_edges = []
    for k in range(len(edges)):
        a, b = edges[eperm[k]]
        a, b = pi[a], pi[b]
        if rev[k]:
            a, b = b, a
        new_edges.append((a, b))
    inv = {old: k for k, old in enumerate(eperm)}
    new_paths = [[inv[e] for e in P] for P in paths]
    new_ext = [pi[v] for v in ext]
    return new_edges, new_paths, new_ext


def cochain_value(alg, nv, edges, paths, ext, W, args):
    tot = 0
    for s, sg in ce.PERMS[4]:
        ps = [args[s[i]][1] for i in range(4)]
        lam = {v: [0, 0] for v in range(nv - 1)}
        for i, v in enumerate(ext):
            if v != nv - 1:
                lam[v] = list(ps[i])
        c = theta_const(nv, edges, paths, lam)
        cs = [args[s[i]][0] for i in range(4)]
        r = W
        for cc in reversed(cs):
            r = np.tensordot(r, np.array(cc, dtype=object), axes=([r.ndim - 1], [0]))
        tot += sg * c * int(r)
    return tot


graphs = {
    "theta(1,3,3)": (6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [[0], [1, 2, 3], [4, 5, 6]], [1, 2, 3, 4]),
    "theta(2,2,3)": (6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [[0, 1], [2, 3], [4, 5, 6]], [1, 2, 3, 4]),
}
for alg in (le.sl(2), le.gl(2)):
    for name, (nv, edges, paths, ext) in graphs.items():
        W0, sc, sg0 = le.orientation_rule(alg, nv, edges, ext)
        samples = [ce.random_args(alg, rng, 4, H) for _ in range(2)]
        ref = [cochain_value(alg, nv, edges, paths, ext, W0, a) for a in samples]
        res = []
        for _ in range(4):
            e2, p2, x2 = relabel(nv, edges, paths, ext)
            W2, sc2, sg2 = le.orientation_rule(alg, nv, e2, x2)
            val = [cochain_value(alg, nv, e2, p2, x2, W2, a) for a in samples]
            res.append(val == ref)
        print("%s %s: alt(K W) on 2 samples equal for 4 random relabellings/reversals/edge orders: %s; reference nonzero: %s" % (
            alg.name, name, res, any(r != 0 for r in ref)), flush=True)
