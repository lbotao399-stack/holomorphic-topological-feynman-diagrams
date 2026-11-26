"""(d) Kernels of theta(1,3,3), theta(2,2,3) at (H,T) = (2,1) with one vertex v carrying the J_1 insertion
[q_h1 q_h2] (h1, h2 = first two half-edges of the local order L'_v; q as in common.py), exact
(htk.kernel_exact, Dirichlet integrals, PiNumber; every kernel = rational polynomial x pi).
Checks: (i) colour_weight (all f) == anomaly-structure formal_gl.colour_weight_formal;
(ii) K without insertion = -(pi/20)[p1p2][p3p4], (pi/30)[p1p2][p3p4];
(iii) pair independence: [q_h1 q_h2] = [q_h2 q_h3] = [q_h3 q_h1] kernels equal modulo sum p = 0;
(iv) brute-force density check of htk at 2 random points (check=2) for every kernel.
Saves data/kernels.pkl."""
import time
from common import *
import formal_gl as fg

t0 = time.time()
br = lambda a, b: a[0] * b[1] - a[1] * b[0]
sub3 = {PSYM[3][c]: -(PSYM[0][c] + PSYM[1][c] + PSYM[2][c]) for c in range(2)}
OUTK = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, lam = graph_data(pqr)
    print("theta%s: edges %s, ext %s, sgn(L->L') = %+d" % (pqr, G.edges, ext, orientation_sign(G.nv, G.edges, ext)))
    W0 = colour_weight(G.nv, G.edges, ext)
    Wref = fg.colour_weight_formal(G.nv, G.edges, ext)
    Wref = {t: {p: Fr(v) for p, v in d.items()} for t, d in Wref.items()}
    print("  colour_weight(all f) == formal_gl.colour_weight_formal:", W0 == Wref)
    K0 = kernel_over_pi(G, lam)
    print("  K/pi =", sp.factor(K0), "; / ([p1p2][p3p4]) =",
          sp.simplify(K0 / (br(PSYM[0], PSYM[1]) * br(PSYM[2], PSYM[3]))))
    loc = local_orders(G.nv, G.edges, ext)
    Kv = {}
    for v in range(G.nv):
        vals = []
        for which in [(0, 1), (1, 2), (2, 0)]:
            Kp = kernel_over_pi(G, lam, insertion_bracket(G, ext, v, which))
            vals.append(Kp)
        same = all(sp.expand(sp.expand(x - vals[0]).subs(sub3, simultaneous=True)) == 0 for x in vals[1:])
        Kv[v] = vals[0]
        print("  v=%d local order %s: K^(v)/pi has %d terms, degree %s; pair (0,1) = (1,2) = (2,0) mod sum p: %s"
              % (v, loc[v], len(sp.Poly(vals[0], *[PSYM[i][c] for i in range(4) for c in range(2)]).terms()),
                 sp.Poly(vals[0], *[PSYM[i][c] for i in range(4) for c in range(2)]).total_degree(), same), flush=True)
        f = sp.factor(sp.expand(vals[0].subs(sub3, simultaneous=True)))
        print("      K^(v)/pi (p4 = -p1-p2-p3) =", f)
    OUTK[pqr] = dict(edges=G.edges, ext=ext, K0=K0, Kv=Kv, local=loc)
save(OUTK, 'kernels.pkl')
print("time %.1fs" % (time.time() - t0))
