"""closedness.dF_exact against run-1 cochain.d_CE (read-only import) on a non-closed alternating cochain:
C = alt([p0 p1][p0 p2] W_133) for the run-1 sl2 and sl3 bases; comparison of d C(b_theta, ...)."""
import sys, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import numpy as np, sympy as sp
import algebras as al, core, closedness as cl
from algebras import perm_sign
import lie, cochain as cc
mom = cc.Momenta(6, 2)
Pm = mom.p
Bs = lambda a, b: Pm[a][0] * Pm[b][1] - Pm[a][1] * Pm[b][0]
gens5 = mom.gens(5)
mon_polys = [sp.Poly(sp.expand(Bs(*cl.PAIRS[i]) * Bs(*cl.PAIRS[j])).subs(
    {Pm[4][c]: -sum(Pm[t][c] for t in range(4)) for c in range(2)}), *gens5).as_dict() for (i, j) in cl.MON20]
for N in (2, 3):
    g = lie.sl(N)
    L = al.LA('run1-sl%d' % N, g.fu, g.kappa, 'run-1 trace form', cartan=list(range(N * (N - 1), g.dim)))
    T = core.Tensors(L); th = L.theta; d = L.dim
    W = T.W('133', [None] * 4)
    Tmu = np.zeros((6,) + (d,) * 4, dtype=np.int64)
    for sig in itertools.permutations(range(4)):
        c = core.quad(core.br(sig[0], sig[1]), core.br(sig[0], sig[2]))
        V = np.transpose(W, list(np.argsort(sig)))
        for mu in range(6):
            if c[mu]:
                Tmu[mu] += perm_sign(sig) * int(c[mu]) * V
    acc = cl.dF_exact(L, T, Tmu)
    mine = {}
    for nu in range(20):
        for m, c in mon_polys[nu].items():
            mine[m] = mine.get(m, 0) + int(c) * acc[nu]
    C = cc.alternate(mom, Bs(0, 1) * Bs(0, 2), W, 4)
    dC = cc.d_CE(mom, g, C, 4)
    keys = set(mine) | set(dC)
    z = np.zeros((d,) * 4, dtype=np.int64)
    ok = all(np.array_equal(mine.get(m, z), dC[m][th] if m in dC else z) for m in keys)
    print("sl%d: d C nonzero: %s ; dF_exact(C)(b_theta,...) = run-1 d_CE(C)(b_theta,...) on all %d monomials: %s"
          % (N, acc.any(), len(keys), ok), flush=True)

# non-closed input (formula applied to a non-alternating cochain): F_133 coefficient tensors with polynomial
# labels cyclically shifted; both codes evaluate sum_{i<j} (-1)^{i+j} C([x_i,x_j], ...).
gens4 = mom.gens(4)
B4 = lambda a, b: Pm[a][0] * Pm[b][1] - Pm[a][1] * Pm[b][0]
sub4 = {Pm[3][c]: -sum(Pm[t][c] for t in range(3)) for c in range(2)}
u, v, w = B4(0, 1), B4(0, 2), B4(1, 2)
q4 = [sp.Poly(sp.expand(x), *gens4).as_dict() for x in [u*u, v*v, w*w, u*v, u*w, v*w]]
for N in (2, 3):
    g = lie.sl(N)
    L = al.LA('run1-sl%d' % N, g.fu, g.kappa, 'run-1 trace form', cartan=list(range(N * (N - 1), g.dim)))
    T = core.Tensors(L); th = L.theta; d = L.dim
    Tmu, W = cl.T_mu_full(T, '133')
    bad = np.roll(Tmu, 1, axis=0)
    acc = cl.dF_exact(L, T, bad)
    mine = {}
    for nu in range(20):
        for m, c in mon_polys[nu].items():
            mine[m] = mine.get(m, 0) + int(c) * acc[nu]
    C = {}
    for mu in range(6):
        for m, c in q4[mu].items():
            C[m] = C.get(m, 0) + int(c) * bad[mu]
    dC = cc.d_CE(mom, g, C, 4)
    keys = set(mine) | set(dC)
    z = np.zeros((d,) * 4, dtype=np.int64)
    ok = all(np.array_equal(mine.get(m, z), dC[m][th] if m in dC else z) for m in keys)
    print("sl%d: shifted-label input: d C nonzero: %s ; dF_exact = run-1 d_CE at b_theta on all %d monomials: %s"
          % (N, acc.any(), len(keys), ok), flush=True)
