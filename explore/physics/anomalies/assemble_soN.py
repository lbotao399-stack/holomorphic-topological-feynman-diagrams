"""T_1, T_3 for so_N (N = 5..11) and sp_2n (2n = 4, 6, 8) by explicit contraction; fit in tr x^4, (tr x^2)^2;
then interpolate the coefficients as polynomials in N (O(N) family) and check sp_2n against N -> -2n."""
import json, random
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import colour as C
rng = random.Random(37)

def aut(nv, edges):
    g = nx.Graph(); g.add_nodes_from(range(nv)); g.add_edges_from(edges)
    return sum(1 for _ in GraphMatcher(g, g).isomorphisms_iter())

recs = json.load(open('kernels_n3.json'))
box = dict(idx='box', nv=4, edges=[[0, 1], [1, 2], [2, 3], [3, 0]], legs=[0, 1, 2], bare={'1': '-1/6'})
AUT = {r['idx']: aut(r['nv'], [tuple(e) for e in r['edges']]) for r in recs}

def fit(alg):
    rows, r1, r3, r3z = [], [], [], []
    for _ in range(4):
        x = [rng.randint(-2, 2) for _ in range(alg.dim)]
        M = sum((np.array(alg.B[e], dtype=object) * x[e] for e in range(alg.dim)), np.zeros_like(alg.B[0]))
        M2 = M.dot(M)
        rows.append([int(np.trace(M2.dot(M2))), int(np.trace(M2)) ** 2])
        w = Fr(C.einsum_W(alg, 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3], [x] * 4), alg.kinv_scale ** 4)
        r1.append(Fr(-1, 6) * w / 8)
        t, tz = Fr(0), Fr(0)
        for r in recs:
            c = {k: Fr(v) for k, v in r['bare'].items()}
            if not c:
                continue
            edges = [tuple(e) for e in r['edges']]
            w = Fr(C.einsum_W(alg, r['nv'], edges, r['legs'] + [r['nv'] - 1], [x] * 4), alg.kinv_scale ** len(edges))
            t += c.get('1', 0) * w / AUT[r['idx']]
            tz += c.get('z3', 0) * w / AUT[r['idx']]
        r3.append(t); r3z.append(tz)
    A = sp.Matrix(rows)
    out = []
    for b in (r1, r3, r3z):
        sol = A[:2, :].LUsolve(sp.Matrix(b[:2]))
        ok = all(sp.simplify(sum(A[i, j] * sol[j] for j in range(2)) - b[i]) == 0 for i in range(4))
        out.append(([sp.Rational(s) for s in sol], ok))
    return out

res = {}
for N in range(5, 12):
    o = fit(C.so(N))
    res[N] = o
    print('so%d: T_1 = %s (fit ok %s); T_3 = %s (ok %s); T_3 zeta3 part = %s (ok %s)' % (N, o[0][0], o[0][1], o[1][0], o[1][1], o[2][0], o[2][1]), flush=True)
t = sp.Symbol('t')
for k, nm in ((0, 'tr x^4'), (1, '(tr x^2)^2')):
    pts = [(N, res[N][1][0][k]) for N in range(5, 11)]
    poly = sp.interpolate(pts, t)
    print('T_3 coefficient of %s for so(t), interpolated on N=5..10: %s ; prediction N=11 %s, computed %s' % (nm, sp.factor(poly), poly.subs(t, 11), res[11][1][0][k]))
    pts1 = [(N, res[N][0][0][k]) for N in range(5, 11)]
    print('T_1 coefficient of %s: %s' % (nm, sp.factor(sp.interpolate(pts1, t))))
for n2 in (4, 6, 8):
    o = fit(C.sp_(n2))
    print('sp%d: T_1 = %s; T_3 = %s (ok %s); zeta3 part %s' % (n2, o[0][0], o[1][0], o[1][1], o[2][0]), flush=True)
