"""Periods I_q(G; alpha) via graphical functions in d = 2q for graphs with a sequential decompletion.

Position space: A_G(0,e) = int prod_{v != h0,h1} d^d x_v / pi^{d/2} prod_e |x_e|^{-2 beta_e}  (h0 = 0, h1 = e).
Dictionary (beta = q - alpha):  I_q(G; alpha) = prod_e Gamma(beta_e) / Gamma(q) * A_G(0, e).
Completion: Gbar = G + vertex 'oo', edges v-oo of weight 2q - deg_beta(v); A_{Gbar - c} = A_G for every c.
Sequential weighted graph G': hubs h0, h1, G' - {h0,h1} = path v_1..v_k, path edges 1 <= beta <= lam
(chain of alpha = q - beta edges of weight lam, each one append); hub edges any integer weight.
Chain: A_{G'} = A_{G'_sub} prod_{path e} Gamma(lam)^{alpha_e} Gamma(alpha_e) / Gamma(beta_e).
"""
import time
from fractions import Fraction as F
from math import factorial

import networkx as nx

import gf
import mzv


def gamma_int(n):
    assert n >= 1, n
    return factorial(n - 1)


def merge(wedges):
    """merge parallel edges (sum weights), drop zero weights"""
    d = {}
    for u, v, b in wedges:
        k = frozenset((u, v))
        assert len(k) == 2
        d[k] = d.get(k, 0) + b
    return {k: b for k, b in d.items() if b}


def decompositions(wd):
    G = nx.Graph()
    for k in wd:
        G.add_edge(*tuple(k))
    out = []
    V = list(G.nodes)
    for h0 in V:
        for h1 in V:
            if h0 == h1:
                continue
            H = G.copy()
            H.remove_nodes_from([h0, h1])
            n = H.number_of_nodes()
            if n == 0 or H.number_of_edges() != n - 1 or not nx.is_connected(H):
                continue
            if any(dd > 2 for _, dd in H.degree()):
                continue
            ends = [v for v, dd in H.degree() if dd <= 1]
            path = [ends[0]]
            while len(path) < n:
                path.append(next(u for u in H.neighbors(path[-1]) if u not in path))
            out.append((h0, h1, path))
            if n > 1:
                out.append((h0, h1, path[::-1]))
    return out


def plan(wd, q, h0, h1, path):
    lam = q - 1

    def b(u, v):
        return wd.get(frozenset((u, v)), 0)
    hub = [(b(v, h0), b(v, h1)) for v in path]
    chains = []
    for i in range(len(path) - 1):
        be = b(path[i], path[i + 1])
        if not 1 <= be <= lam:
            return None
        chains.append(q - be)
    cur = [e for e in [('z', '0', hub[0][0]), ('z', '1', hub[0][1])] if e[2]]
    cnt = 0
    ops = []
    for i, c in enumerate(chains):
        for _ in range(c):
            cnt += 1
            bd = gf.bounds_path(cur, lam) if cur else {'inf': 0}
            if -bd['inf'] <= 1:
                return None
            name = 'u%d' % cnt
            cur = [(name if u == 'z' else u, name if v == 'z' else v, bb) for u, v, bb in cur]
            cur.append((name, 'z', lam))
            ops.append(('append', gf.bounds_path(cur, lam)))
        h = hub[i + 1]
        if h[0] or h[1]:
            ops.append(('edge', h))
            cur += [e for e in [('z', '0', h[0]), ('z', '1', h[1])] if e[2]]
    bd = gf.bounds_path(cur, lam)
    if -bd['inf'] <= lam + 1:
        return None
    for s in ('0', '1'):
        # integrability of the final integrand at 0, 1: f ~ |z-s|^{2 M_s}, need 2 M_s + 2 lam + 2 > 0
        if bd[s] + lam + 1 <= 0:
            return None
    chainfac = F(1)
    for i, c in enumerate(chains):
        chainfac *= F(gamma_int(lam) ** c * gamma_int(c), gamma_int(q - c))
    return {'first': hub[0], 'ops': ops, 'chainfac': chainfac, 'lam': lam, 'path': path, 'hubs': (h0, h1),
            'appends': cnt}


def run(pl, verbose=False):
    lam = pl['lam']
    t = time.time()
    f = gf.VD.rational(lam, gf.edge_factor(*pl['first']))
    for op in pl['ops']:
        if op[0] == 'append':
            f = gf.append_particular(f)
            f = f.symmetrize()
            f, free = gf.fix_kernel(f, op[1], verbose=verbose)
            if verbose:
                print('   append: bounds %s, weight %d, size %d (%.1fs)' % (op[1], f.maxweight(), f.size(),
                                                                          time.time() - t), flush=True)
        else:
            f = f.mul_edge(gf.edge_factor(*op[1]))
    A = gf.period(f)
    return mzv.scal(A, pl['chainfac'])


def completion(wd, q):
    deg = {}
    for k, b in wd.items():
        for v in k:
            deg[v] = deg.get(v, 0) + b
    out = dict(wd)
    for v, dg in deg.items():
        if 2 * q - dg:
            out[frozenset((v, 'oo'))] = 2 * q - dg
    return out


def plans_for(edges, alphas, q, use_completion=True):
    wd = merge([(u, v, q - a) for (u, v), a in zip(edges, alphas)])
    cands = []
    graphs = [('direct', wd)]
    if use_completion:
        cb = completion(wd, q)
        verts = set()
        for k in cb:
            verts |= set(k)
        for c in sorted(verts, key=str):
            if c == 'oo':
                continue
            g2 = {k: b for k, b in cb.items() if c not in k}
            graphs.append(('decomplete %s' % c, g2))
    for tag, g in graphs:
        for h0, h1, path in decompositions(g):
            pl = plan(g, q, h0, h1, path)
            if pl is not None:
                pl['route'] = tag
                cands.append(pl)
    cands.sort(key=lambda pl: (pl['appends'], pl['route'] != 'direct', str(pl['path'])))
    return cands


def period(edges, alphas, q, verbose=False, choice=0, use_completion=True):
    cands = plans_for(edges, alphas, q, use_completion)
    if not cands:
        raise ValueError('no convergent sequential plan')
    pl = cands[choice]
    if verbose:
        print('  plan: %s hubs %s path %s, %d appends' % (pl['route'], pl['hubs'], pl['path'], pl['appends']))
    A = run(pl, verbose)
    pref = F(1, gamma_int(q))
    for a in alphas:
        pref *= gamma_int(q - a)
    return mzv.scal(A, pref), pl
