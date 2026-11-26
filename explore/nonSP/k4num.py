"""Rough numerics for the K4 subdivisions on C^3 (H=3, T=0).
K = int_{RP^9} D iota_E vol, D from htk.oracle.scaled_density_batch (brute force, no reduction).
Coordinates: chain ratios x (exact Gauss cubature), K4 resistances R(xi), xi = (gA,gB,rA,rB,rC),
Y-Delta chart at the K4 vertex Dc:  R_{Dc v} = 1/g_v,  R_{w x} = (1-r_v)/(r_v g_w g_x)  (v opposite to wx)."""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk.graph import Graph
from htk import oracle
from graphs import K4SUB, skeleton_chains, deg

def setup(edges, Dc_index=0):
    nv = 8
    G = Graph(nv, edges)
    d = deg(edges, nv)
    d2 = [v for v in range(nv) if d[v] == 2]
    lam_all = {d2[0]: [1, 0, 0], d2[1]: [0, 1, 0], d2[2]: [0, 0, 1], d2[3]: [-1, -1, -1]}
    lam = np.array([lam_all.get(v, [0, 0, 0]) for v in range(nv - 1)], dtype=float)
    cubic, chains = skeleton_chains(edges, nv)
    Dc = cubic[Dc_index]
    others = [v for v in cubic if v != Dc]
    return G, lam, lam_all, cubic, chains, Dc, others

def gl01(m):
    x, w = np.polynomial.legendre.leggauss(m)
    return (x + 1) / 2, w / 2

def x_rule(chains, m):
    """nodes (Q, N) of t-multipliers x_e per edge, weights (Q,)"""
    per = []
    for path, es in chains:
        a = len(es)
        if a == 1:
            per.append((np.ones((1, 1)), np.ones(1)))
        elif a == 2:
            s, w = gl01(m)
            per.append((np.stack([1 - s, s], 1), w))
        elif a == 3:
            s, w = gl01(m)
            U, V = np.meshgrid(s, s, indexing='ij'); WU, WV = np.meshgrid(w, w, indexing='ij')
            U, V, W = U.ravel(), V.ravel(), (WU * WV).ravel() * (1 - U.ravel())
            per.append((np.stack([U, (1 - U) * V, (1 - U) * (1 - V)], 1), W))
        else:
            raise NotImplementedError
    nodes, wts = [], []
    for combo in itertools.product(*[range(len(p[1])) for p in per]):
        xv = np.zeros(10); wt = 1.0
        for (path, es), p, k in zip(chains, per, combo):
            for j, e in enumerate(es):
                xv[e] = p[0][k, j]
            wt *= p[1][k]
        nodes.append(xv); wts.append(wt)
    return np.array(nodes), np.array(wts)

def chainR_index(chains):
    """map K4 edge (frozenset of endpoints) -> chain index"""
    return {frozenset((path[0], path[-1])): i for i, (path, es) in enumerate(chains)}

def R_of_xi(xi, chains, Dc, others):
    gA, gB, rA, rB, rC = [xi[:, k] for k in range(5)]
    gC = 1 - gA - gB
    g = {others[0]: gA, others[1]: gB, others[2]: gC}
    r = {others[0]: rA, others[1]: rB, others[2]: rC}
    idx = chainR_index(chains)
    R = np.zeros((xi.shape[0], 6))
    for v in others:
        R[:, idx[frozenset((Dc, v))]] = 1 / g[v]
        w, x = [u for u in others if u != v]
        R[:, idx[frozenset((w, x))]] = (1 - r[v]) / (r[v] * g[w] * g[x])
    return R, g, r

def jac_logR(xi, chains, Dc, others):
    """det[1; d log R / d xi] (6x6)"""
    gA, gB, rA, rB, rC = [xi[:, k] for k in range(5)]
    B = xi.shape[0]
    idx = chainR_index(chains)
    # d log g_v / d(gA, gB):  g_C = 1 - gA - gB
    dlg = {others[0]: np.stack([1 / gA, 0 * gA], 1), others[1]: np.stack([0 * gA, 1 / gB], 1),
           others[2]: np.stack([-1 / (1 - gA - gB), -1 / (1 - gA - gB)], 1)}
    rr = {others[0]: (rA, 2), others[1]: (rB, 3), others[2]: (rC, 4)}
    J = np.zeros((B, 6, 6)); J[:, 0, :] = 1
    for v in others:
        c = idx[frozenset((Dc, v))]
        J[:, 1:3, c] = -dlg[v]
        w, x = [u for u in others if u != v]
        c = idx[frozenset((w, x))]
        J[:, 1:3, c] = -dlg[w] - dlg[x]
        rv, k = rr[v]
        J[:, 1 + k, c] = -1 / (1 - rv) - 1 / rv
    return np.linalg.det(J)

def SR(Rb, G, lam, chains, xn, xw, st):
    """S_R(R) = int_x S(t(x,R)) / prod x_e dx  for a batch of R (B,6)"""
    B = Rb.shape[0]
    Rt = np.zeros((B, 10))
    for i, (path, es) in enumerate(chains):
        for e in es:
            Rt[:, e] = Rb[:, i]
    Q = xn.shape[0]
    t = (Rt[:, None, :] * xn[None, :, :]).reshape(B * Q, 10)
    S = oracle.scaled_density_batch(G, t, lam, 0, _setup=st).reshape(B, Q)
    return (S / np.prod(xn, axis=1)[None, :]) @ xw

def density_xi(xi, G, lam, chains, Dc, others, xn, xw, st):
    R, g, r = R_of_xi(xi, chains, Dc, others)
    return SR(R, G, lam, chains, xn, xw, st) * np.abs(jac_logR(xi, chains, Dc, others))

if __name__ == "__main__":
    from scipy.stats import qmc
    which = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    npts = int(sys.argv[2]) if len(sys.argv) > 2 else 2**14
    m = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    name, edges = K4SUB[which]
    G, lam, lam_all, cubic, chains, Dc, others = setup(edges)
    st = oracle._Setup(G)
    xn, xw = x_rule(chains, m)
    print(name, "Dc", Dc, "others", others, "x nodes", len(xw))
    rng = np.random.RandomState(1)
    # cubature exactness check at a few R points
    Rtest = np.exp(rng.normal(size=(3, 6)))
    for mm in (m, m + 1, m + 3):
        a, b = x_rule(chains, mm)
        print("  m=%d S_R:" % mm, SR(Rtest, G, lam, chains, a, b, st))
    vals = []
    for rep in range(4):
        sob = qmc.Sobol(5, scramble=True, seed=rep).random(npts)
        p, q = sob[:, 0], sob[:, 1]
        xi = np.stack([p, (1 - p) * q, sob[:, 2], sob[:, 3], sob[:, 4]], 1)
        f = density_xi(xi, G, lam, chains, Dc, others, xn, xw, st) * (1 - p)
        vals.append(f.mean())
    vals = np.array(vals)
    print("  c estimate", vals.mean(), "+-", vals.std() / np.sqrt(len(vals)), vals)
