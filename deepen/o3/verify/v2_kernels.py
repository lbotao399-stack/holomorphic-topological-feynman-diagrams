"""J_1-vertex kernels of theta(1,3,3), theta(2,2,3) at (H,T) = (2,1), own insertion construction.
K^(v)(lam;0) = int Omega_G [q_h1 q_h2], (h1,h2) = first two half-edges of the local order (choice A) at v.
Checks: base kernels; insertion with the other two pairs (h2,h3), (h3,h1) (conservation makes these equal as
integrands, so this only checks the construction); local order choice B gives sign sgn(perm) x kernel.
Saves verify/data/kernels.pkl."""
import pickle, os, time
import sympy as sp
from vgraph import *

t0 = time.time()
os.makedirs('data', exist_ok=True)
br = lambda a, b: a[0] * b[1] - a[1] * b[0]
sub4 = {PS[3][c]: -(PS[0][c] + PS[1][c] + PS[2][c]) for c in range(2)}
b = lambda i, j: sp.Symbol('b%d%d' % (i + 1, j + 1))
def to_b(expr):
    """write a (p4-eliminated) SL2-invariant polynomial in b12, b13, b23 (b_ij = [p_i p_j])."""
    e = sp.expand(expr.subs(sub4, simultaneous=True))
    B12, B13, B23 = sp.symbols('b12 b13 b23')
    # p1 = (1,0), p2 = (x,1)... generic substitution: solve by matching in a random basis
    import itertools
    cand = [B12 ** i * B13 ** j * B23 ** k for i in range(4) for j in range(4) for k in range(4) if i + j + k == sp.Poly(e, *[PS[i][c] for i in range(3) for c in range(2)]).total_degree() // 2] if e != 0 else []
    if e == 0:
        return 0
    cs = sp.symbols('c0:%d' % len(cand))
    bb = {B12: br(PS[0], PS[1]), B13: br(PS[0], PS[2]), B23: br(PS[1], PS[2])}
    ansatz = sum(c * m.subs(bb, simultaneous=True) for c, m in zip(cs, cand))
    eqs = sp.Poly(sp.expand(ansatz - e), *[PS[i][c] for i in range(3) for c in range(2)]).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    assert sol, "not SL2-invariant?"
    return sp.factor(sum(sol[0].get(c, 0) * m for c, m in zip(cs, cand)))

OUT = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, lam = theta_data(pqr)
    locA = local_orders(G, ext, 'A')
    locB = local_orders(G, ext, 'B')
    sA, sB = wick_sign(G, ext, locA), wick_sign(G, ext, locB)
    print("theta%s edges %s ext %s  sigma(A) = %+d  sigma(B) = %+d" % (pqr, G.edges, ext, sA, sB))
    K0 = kernel_over_pi(G, lam)
    print("  K/pi =", sp.factor(K0), "   K/pi / ([p1p2][p3p4]) =", sp.simplify(K0 / (br(PS[0], PS[1]) * br(PS[2], PS[3]))))
    Kv = {}
    for v in range(G.nv):
        vals = [kernel_over_pi(G, lam, insertion(G, locA, v, pr)) for pr in [(0, 1), (1, 2), (2, 0)]]
        same = all(sp.expand((x - vals[0]).subs(sub4, simultaneous=True)) == 0 for x in vals[1:])
        KB = kernel_over_pi(G, lam, insertion(G, locB, v, (0, 1)))
        # sign of the local-order change at v: sgn(perm A_v -> B_v)
        pos = {h: i for i, h in enumerate(locA[v])}
        sgn_v = seq_sign([pos[h] for h in locB[v]])
        covB = sp.expand((KB - sgn_v * vals[0]).subs(sub4, simultaneous=True)) == 0
        Kv[v] = vals[0]
        print("  v=%d L_A=%s: K^(v)/pi = %s;  pairs equal mod sum p: %s;  K_B = sgn(A->B) K_A: %s (sgn %+d)"
              % (v, locA[v], to_b(vals[0]), same, covB, sgn_v), flush=True)
    OUT[pqr] = dict(edges=G.edges, ext=ext, K0=K0, Kv=Kv, locA=locA, sigmaA=sA)
pickle.dump(OUT, open('data/kernels.pkl', 'wb'))
print("time %.1fs" % (time.time() - t0))
