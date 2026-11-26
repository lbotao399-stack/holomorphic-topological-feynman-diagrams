"""Dimension of the response image of t_S -> (Lr(t_S)) and t_S -> (Lr, h) for k-terminal subnetworks:
generic rank of the Jacobian (exact, random rational point), fiber dimension f_S = N_S - rank,
degree of rho_S = pi_* Omega_S = n_S D - f_S (T = 0: D = H).  Rank-structured response of a star:
c_ab c_cd - c_ac c_bd = 0 (exact check).  Kirchhoff star-mesh identity K_Gamma = G K_{Gamma'}(x + g g^T / G)."""
import random
from fractions import Fraction as Fr
from itertools import combinations

import locality as L
from htk.graph import Graph
from htk.linalg import det


def rank(M):
    M = [[Fr(x) for x in r] for r in M]
    rk, rows, cols = 0, len(M), len(M[0]) if M else 0
    for c in range(cols):
        p = next((r for r in range(rk, rows) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(rows):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / M[rk][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[rk])]
        rk += 1
    return rk


def response_data(edges, terminals, rng):
    """S given by its own edges (vertices labelled freely), terminals list.  Embed S into a graph
    Gamma = S + a star on the terminals (outer vertex = root) so that Split applies."""
    V = sorted(set(v for e in edges for v in e))
    r = max(V) + 1
    G = Graph(r + 1, list(edges) + [(a, r) for a in terminals])
    sp = L.Split(G, list(range(len(edges))), terminal_order=list(terminals))
    t = L.rand_times(G.ne, rng)
    K = sp.kron(t)
    k, m = sp.k, sp.m
    rowsL, rowsLh = [], []
    for e in sp.S:
        dLr, dE, dZg = sp.dkron(t, K, e)
        cvec = [-dLr[a][b] for a, b in combinations(range(k), 2)]
        hvec = [dE[v][a] for v in range(m) for a in range(k - 1)]
        rowsL.append(cvec)
        rowsLh.append(cvec + hvec)
    # include the scale direction automatically (rows are d/dt_e)
    rL, rLh = rank(rowsL), rank(rowsLh)
    return dict(k=k, m=m, NS=len(sp.S), nS=sp.nS, rank_L=rL, rank_Lh=rLh, full=k * (k - 1) // 2, Lr=K['Lr'])


CASES = {
    'star k=3': ([(0, 1), (0, 2), (0, 3)], [1, 2, 3]),
    'star k=4': ([(0, 1), (0, 2), (0, 3), (0, 4)], [1, 2, 3, 4]),
    'star k=5': ([(0, a) for a in range(1, 6)], [1, 2, 3, 4, 5]),
    'star k=3, one leg subdivided': ([(0, 1), (0, 2), (0, 4), (4, 3)], [1, 2, 3]),
    'star k=3, all legs subdivided': ([(0, 4), (4, 1), (0, 5), (5, 2), (0, 6), (6, 3)], [1, 2, 3]),
    'paired stars K_{2,3}': ([(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)], [2, 3, 4]),
    'K4 minus vertex-star (triangle 1-2-3 + v0 adj 1,2,3), k=3 on {1,2,3}?': ([(1, 2), (1, 3), (2, 3), (0, 1), (0, 2), (0, 3)], [1, 2, 3]),
    'wheel-like: two adjacent centres, k=4': ([(0, 1), (0, 2), (0, 5), (5, 3), (5, 4)], [1, 2, 3, 4]),
    'two adjacent centres + rim edge, k=4': ([(0, 1), (0, 2), (0, 5), (5, 3), (5, 4), (0, 3)], [1, 2, 3, 4]),
    'K4 on terminals, k=4': ([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)], [1, 2, 3, 4]),
    'K_{2,4} (two deg-4 centres), k=4': ([(0, a) for a in range(2, 6)] + [(1, a) for a in range(2, 6)], [2, 3, 4, 5]),
    'K_{3,4}, k=4': ([(c, a) for c in range(3) for a in range(3, 7)], [3, 4, 5, 6]),
    'bridge K4-pq, k=2': ([(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], [0, 1]),
    'bridge, k=2, two subdivisions': ([(0, 2), (0, 3), (1, 2), (1, 4), (4, 3), (2, 5), (5, 3)], [0, 1]),
}

if __name__ == '__main__':
    rng = random.Random(17)
    print('%-62s %3s %3s %4s %3s %6s %7s %5s' % ('S', 'k', 'm', 'N_S', 'n_S', 'rk(Lr)', 'rk(Lr,h)', 'cone'))
    for name, (ed, term) in CASES.items():
        d = response_data(ed, term, rng)
        print('%-62s %3d %3d %4d %3d %6d %7d %5d' % (name, d['k'], d['m'], d['NS'], d['nS'], d['rank_L'], d['rank_Lh'], d['full']))
    # rank-structure of the star response: c_ab c_cd = c_ac c_bd
    for k in (4, 5):
        d = response_data([(0, a) for a in range(1, k + 1)], list(range(1, k + 1)), rng)
        Lr = d['Lr']
        c = {(a, b): -Lr[a][b] for a in range(k) for b in range(k) if a != b}
        ok = all(c[(a, b)] * c[(cc, dd)] == c[(a, cc)] * c[(b, dd)]
                 for a, b, cc, dd in [(0, 1, 2, 3), (0, 2, 1, 3), (1, 2, 0, 3)])
        print('star k=%d: c_ab c_cd = c_ac c_bd at a random point: %s' % (k, ok))
    # Kirchhoff star-mesh identity for K5 (v = 0 eliminated), conductances x (K4) and g (star), random rational
    import sympy as sp
    x = sp.symbols('x0:6', positive=True)
    gg = sp.symbols('g0:4', positive=True)
    K4E = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

    def kirchhoff(nv, edges, cond):
        Lm = sp.zeros(nv - 1, nv - 1)            # ground the last vertex
        for (a, b), c in zip(edges, cond):
            for u, w in ((a, b), (b, a)):
                if u < nv - 1:
                    Lm[u, u] += c
                    if w < nv - 1:
                        Lm[u, w] -= c
        return sp.expand(Lm.det())
    K5 = kirchhoff(5, [(0, a) for a in range(1, 5)] + K4E, list(gg) + list(x))
    G = sum(gg)
    xp = [x[i] + gg[a - 1] * gg[b - 1] / G for i, (a, b) in enumerate(K4E)]
    rel = [(a - 1, b - 1) for a, b in K4E]
    K4 = kirchhoff(4, rel, xp)
    print('Kirchhoff star-mesh K_{K5}(x, g) = G * K_{K4}(x + g g^T / G):', sp.simplify(K5 - G * K4) == 0)
