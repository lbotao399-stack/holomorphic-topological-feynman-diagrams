"""rho for the two-loop thetas with a current insertion y^i_e (z-derivative of K(lambda;z)) on a fixed edge,
since K_{theta(2,2,3)}(lambda;0) = 0 on sum lambda = 0 at H = 3 (the bare comparison in norm_check.py is vacuous).
Vertex relabellings (root moves) keep edge identities, so the inserted kernel has the same covariance as the bare one."""
import random, sys, io, contextlib
from fractions import Fraction as Fr
with contextlib.redirect_stdout(io.StringIO()):
    import norm_check as N
import htk
import colour as C
rng = random.Random(3)

def kern_ins(nv, edges, src, e, i):
    full = [src.get(v, [0, 0, 0]) for v in range(nv)]
    lam = {v: full[v] for v in range(nv - 1)}
    return {'1': Fr(str(htk.kernel_exact(htk.Graph(nv, edges), lam, powers={(e, i): 1})))}

for name, nv, edges, ext in [('theta133', 6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]),
                             ('theta223', 6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])]:
    for (e, i) in [(0, 0), (3, 1), (6, 2)]:
        ps = [[rng.randint(-5, 5) for _ in range(3)] for _ in range(3)]
        ps.append([-sum(p[a] for p in ps) for a in range(3)])
        mom = {v: ps[k] for k, v in enumerate(ext)}
        N.engine_kernel = lambda nv_, E_, src_, e=e, i=i: kern_ins(nv_, E_, src_, e, i)
        Fr_, _ = N.F_rule(nv, edges, ext, mom)
        rhos = [N.ratio(N.F_Q(nv, edges, r, mom, ext), Fr_) for r in range(nv)]
        nz = any(Fr_['1'].values())
        print('%s insertion y^%d_%d: F_rule nonzero %s; rho over 6 roots %s' % (name, i, e, nz, rhos), flush=True)
