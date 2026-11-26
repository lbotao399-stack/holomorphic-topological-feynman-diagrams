"""T_1, T_3 (and the zeta(3) part of T_3) for so_N, sp_N by explicit exact contraction (einsum, integers),
x = random integer combination of the basis, fitted to tr x^4, (tr x^2)^2 (and Pf x for so8).
T_n(x) = sum_Gamma c_Gamma |Aut|^{-1} W_Gamma(x,x,x,x), kappa = trace form of the defining representation."""
import json, sys, random
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import colour as C
rng = random.Random(31)

def aut(nv, edges):
    g = nx.Graph(); g.add_nodes_from(range(nv)); g.add_edges_from(edges)
    return sum(1 for _ in GraphMatcher(g, g).isomorphisms_iter())

recs = json.load(open('kernels_n3.json'))
box = dict(idx='box', nv=4, edges=[[0, 1], [1, 2], [2, 3], [3, 0]], legs=[0, 1, 2], bare={'1': '-1/6'})

def Tvals(alg, rec_list, x):
    tot = {}
    for r in rec_list:
        nv, edges, legs = r['nv'], [tuple(e) for e in r['edges']], r['legs']
        c = {k: Fr(v) for k, v in r['bare'].items()}
        if not c:
            continue
        w = Fr(C.einsum_W(alg, nv, edges, legs + [nv - 1], [x] * 4), alg.kinv_scale ** len(edges))
        A = aut(nv, edges)
        for lab, cv in c.items():
            tot[lab] = tot.get(lab, 0) + cv * w / A
    return tot

def pf(m):
    M = sp.Matrix(m)
    n = M.shape[0]
    # Pfaffian by recursion
    def P(idx):
        if not idx:
            return sp.Integer(1)
        i = idx[0]; s = sp.Integer(0)
        for k in range(1, len(idx)):
            j = idx[k]
            rest = idx[1:k] + idx[k + 1:]
            s += (-1) ** (k - 1) * M[i, j] * P(rest)
        return s
    return P(list(range(n)))

for name, alg, extra in [('so5', C.so(5), False), ('so8', C.so(8), True), ('sp4', C.sp_(4), False), ('so6', C.so(6), False)]:
    rows, rhs1, rhs3, rhs3z = [], [], [], []
    nb = 3 if extra else 2
    for _ in range(nb + 2):
        x = [rng.randint(-2, 2) for _ in range(alg.dim)]
        M = sum((np.array(alg.B[e], dtype=object) * x[e] for e in range(alg.dim)), np.zeros_like(alg.B[0]))
        M2 = M.dot(M)
        basis = [int(np.trace(M2.dot(M2))), int(np.trace(M2)) ** 2]
        if extra:
            basis.append(int(pf(M.tolist())))
        rows.append(basis)
        t1 = Tvals(alg, [box], x)
        t3 = Tvals(alg, recs, x)
        rhs1.append(t1.get('1', 0)); rhs3.append(t3.get('1', 0)); rhs3z.append(t3.get('z3', 0))
    A = sp.Matrix(rows)
    def solve(b):
        sol = sp.Matrix(A[:nb, :]).LUsolve(sp.Matrix(b[:nb]))
        resid = [sp.nsimplify(sum(A[i, j] * sol[j] for j in range(nb)) - b[i]) for i in range(len(b))]
        return [str(s) for s in sol], all(r == 0 for r in resid)
    labels = ['tr x^4', '(tr x^2)^2'] + (['Pf x'] if extra else [])
    print('%s: basis %s' % (name, labels))
    s1, ok1 = solve(rhs1); print('   T_1         =', s1, ' fit consistent on extra samples:', ok1)
    s3, ok3 = solve(rhs3); print('   T_3 [1]     =', s3, ' fit consistent:', ok3)
    s3z, ok3z = solve(rhs3z); print('   T_3 [zeta3] =', s3z, ' fit consistent:', ok3z, flush=True)
