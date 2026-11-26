"""Cross-check of core.F_eval against the run-1 cochains (anomaly-structure/code, imported read-only)
on the run-1 sl2, sl3 bases (trace form): omega_133 = -3 F_133, omega_223 = 2 F_223 at slot 0 = b_theta."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import numpy as np, sympy as sp
import algebras as al, core
import lie, c_obstruction as co

print("momentum algebra (reduction identities, rank of the 6 quadratics):", core.check_momentum_algebra())
P = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(6)]
gens = [P[i][c] for i in range(3) for c in range(2)]
B = lambda a, b: P[a][0] * P[b][1] - P[a][1] * P[b][0]
u, v, w = B(0, 1), B(0, 2), B(1, 2)
qpolys = [sp.Poly(sp.expand(q), *gens).as_dict() for q in [u*u, v*v, w*w, u*v, u*w, v*w]]
for N in (2, 3):
    g = lie.sl(N)
    L = al.LA('run1-sl%d' % N, g.fu, g.kappa, 'run-1 trace form', cartan=list(range(N * (N - 1), g.dim)))
    assert L.s == g.kinv_scale, (L.s, g.kinv_scale)
    T = core.Tensors(L)
    th = L.theta
    res, om = co.run(g)
    for G, fac in (('133', -3), ('223', 2)):
        F = core.F_eval(T, G, th)
        mine = {}
        for mu in range(6):
            for m, c in qpolys[mu].items():
                mine[m] = mine.get(m, 0) + int(c) * F[mu]
        theirs = om[(1, 3, 3) if G == '133' else (2, 2, 3)]
        keys = set(mine) | set(theirs)
        ok = all(np.array_equal(fac * mine.get(m, np.zeros_like(F[0])),
                                theirs[m][th] if m in theirs else np.zeros_like(F[0])) for m in keys)
        print("sl%d: run-1 omega_%s(b_theta,.,.,.) = %d * F_%s(b_theta,.,.,.) for all %d monomials: %s"
              % (N, G, fac, G, len(keys), ok), flush=True)

# beta = d([q1q2]^2 f) against run-1 d_CE of alt([q0 q1]^2 f) (= 6 [q0q1]^2 f)
import cochain as cc
mom = cc.Momenta(6, 2)
Pm = mom.p
for N in (2, 3):
    g = lie.sl(N)
    L = al.LA('run1-sl%d' % N, g.fu, g.kappa, 'run-1 trace form', cartan=list(range(N * (N - 1), g.dim)))
    th = L.theta
    e = (Pm[0][0] * Pm[1][1] - Pm[0][1] * Pm[1][0]) ** 2
    eta = cc.alternate(mom, e, g.fl, 3)
    dq = cc.d_CE(mom, g, eta, 3)
    be = core.beta_eval(L, th)
    mine = {}
    for mu in range(6):
        for m, c in qpolys[mu].items():
            mine[m] = mine.get(m, 0) + int(c) * be[mu]
    keys = set(mine) | set(dq)
    ok = all(np.array_equal(6 * mine.get(m, np.zeros_like(be[0])), dq[m][th] if m in dq else np.zeros_like(be[0])) for m in keys)
    print("sl%d: run-1 d_CE(alt([q0q1]^2 f))(b_theta,...) = 6 beta(b_theta,...): %s" % (N, ok))

# kernels from the engine (htk), T = 1
import htk
from htk import theta as thg
lam = [[sp.Symbol('l%d_%d' % (i, c)) for c in range(2)] for i in range(5)]
for pqr in [(1, 3, 3), (2, 2, 3)]:
    Gr, interior = thg(*pqr)
    ext = [x for Pp in interior for x in Pp]
    lm = {x: (lam[ext.index(x)] if x in ext else [0, 0]) for x in range(Gr.nv - 1)}
    K = sp.factor(htk.kernel_exact(Gr, lm, T=1))
    Bl = lambda a, b: lam[a][0] * lam[b][1] - lam[a][1] * lam[b][0]
    print("theta%s: edges %s ext %s  K = %s ; K*60/pi / ([p1p2][p3p4]) = %s"
          % (pqr, Gr.edges, ext, K, sp.simplify(K * 60 / sp.pi / (Bl(0, 1) * Bl(2, 3)))))
