"""Class ratio: the unique c (if any) with omega_133 - c omega_223 in d(C^3), for sl3, gl3
(explicit components) and gl_N, N = 4..10 (formal trace basis).  Exact rational linear solve:
least-squares normal equations of [d(C3) | omega_223] x = omega_133 (consistent iff exact)."""
import sys
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import lie, cochain as cc, formal_gl as fg
from c_obstruction import kernel_poly, trace_tensors

mom = cc.Momenta(6, 2)
P = mom.p
monos3 = cc.all_monos(mom, 3, (2, 2))
monos4 = cc.all_monos(mom, 4, (2, 2))


def solve_ratio(A, v1, v2):
    """find c, x with A x + c v2 = v1 exactly; returns c or None."""
    M = np.concatenate([A, v2[:, None]], axis=1)
    keep = np.any(M != 0, axis=1) | (v1 != 0)
    M, b = M[keep], v1[keep]
    Mo = M.astype(object)
    G = sp.Matrix(Mo.T.dot(Mo).tolist())
    rhs = sp.Matrix(Mo.T.dot(b.astype(object)).tolist())
    sol, params = G.gauss_jordan_solve(rhs)
    sol = sol.subs({p_: 0 for p_ in params})
    x = np.array([int(sp.Rational(s).p) if sp.Rational(s).q == 1 else None for s in sol], dtype=object)
    # verify exactly
    solq = [sp.Rational(s) for s in sol]
    res = [sum(solq[j] * int(M[i, j]) for j in range(M.shape[1]) if M[i, j]) - int(b[i]) for i in range(M.shape[0])]
    if any(r != 0 for r in res):
        return None
    return solq[-1]


for g in [lie.sl(3), lie.gl(3)]:
    ws = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, K = kernel_poly(pqr)
        W, sc = cc.colour_weight(g, G.nv, G.edges, ext)
        ws[pqr] = cc.flatten(cc.alternate(mom, K, W, 4), monos4, g.dim, 4)
    cols = []
    for T in [T for s, T in trace_tensors(g, 3)] + [g.fl]:
        for mu in monos3:
            eta = cc.alternate(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), T, 3)
            cols.append(cc.flatten(cc.d_CE(mom, g, eta, 3), monos4, g.dim, 4) if eta else
                        np.zeros(len(monos4) * g.dim ** 4, dtype=np.int64))
    A = np.stack(cols, axis=1)
    c = solve_ratio(A, ws[(1, 3, 3)], ws[(2, 2, 3)])
    print("%s: omega_133 = c omega_223 + d(eta) with c = %s" % (g.name, c), flush=True)

Wf, Kp = {}, {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, K = kernel_poly(pqr)
    Kp[pqr] = K
    Wf[pqr] = fg.colour_weight_formal(G.nv, G.edges, ext)
om = {pqr: fg.alternate_formal(mom, Kp[pqr], Wf[pqr], 4) for pqr in Wf}
cols3 = []
for tau in fg.all_taus(3):
    for mu in monos3:
        eta = fg.alternate_formal(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), {tau: {0: 1}}, 3)
        cols3.append(fg.at_N(fg.d_formal(mom, eta, 3), 0) if eta else {})
for N in range(4, 11):
    keys = sorted({k for c in cols3 for k in c} | set(om[(1, 3, 3)]) | set(om[(2, 2, 3)]), key=str)
    A = np.array([[c.get(k, 0) for c in cols3] for k in keys], dtype=np.int64)
    v1 = np.array([fg.at_N(om[(1, 3, 3)], N).get(k, 0) for k in keys], dtype=np.int64)
    v2 = np.array([fg.at_N(om[(2, 2, 3)], N).get(k, 0) for k in keys], dtype=np.int64)
    c = solve_ratio(A, v1, v2)
    print("gl_%d (formal): omega_133 = c omega_223 + d(eta) with c = %s" % (N, c), flush=True)
