"""(d) Independent numerical check (htk.oracle, lattice rule in log Schwinger coordinates, no reduction code,
double precision) of two J_1-vertex kernels at numeric rational sources (vertex 0: 0; vertices 1..4: lambda; root:
minus the sum), against the exact spred values:
  theta(2,2,3), v = 3, L' = (x_3, (e4,1), (e5,0)):  K^(3) = lam3^1 K[y^2_4] - lam3^2 K[y^1_4]
  theta(1,3,3), v = 0, L' = ((e0,0), (e1,0), (e4,0)):  K^(0) = [(-y_0)(-y_1)] = K[y^1_0 y^2_1] - K[y^2_0 y^1_1]."""
import time
import numpy as np
from common import *
from htk import oracle

rng_l = [[Fr(3, 2), Fr(-1)], [Fr(1, 2), Fr(2)], [Fr(-2), Fr(1, 3)], [Fr(1), Fr(-3, 2)]]
cases = [((2, 2, 3), 3, [({(4, 1): 1}, ('lam', 3, 0)), ({(4, 0): 1}, ('lam', 3, 1))]),
         ((1, 3, 3), 0, [({(0, 0): 1, (1, 1): 1}, 1), ({(0, 1): 1, (1, 0): 1}, -1)])]
for pqr, v, terms in cases:
    t0 = time.time()
    G, ext, lamS = graph_data(pqr)
    lam = [[Fr(0), Fr(0)] if u not in ext else rng_l[ext.index(u)] for u in range(G.nv - 1)]
    exact = htk.sp_reduce(G, prefer="large").integrate(lam, 1, insertion=insertion_bracket(G, ext, v, lamnum={u: lam[u] for u in range(G.nv - 1)}), check=2)
    ex = float(exact.coeff(1)) * np.pi if hasattr(exact, 'coeff') else float(exact)
    tot, err = 0.0, 0.0
    for powers, w in terms:
        val, e = oracle.kernel(G, [[float(x) for x in r] for r in lam], T=1, powers=powers)
        if isinstance(w, tuple):
            lv = lam[w[1]][w[2]]
            wf = float(lv) * (1 if w[2] == 0 else -1)
        else:
            wf = float(w)
        tot += wf * val
        err += abs(wf) * e
    print("theta%s v=%d: exact %s = %.10f ; oracle %.10f +- %.1e ; rel. diff %.2e ; digits %.1f   (%.0fs)"
          % (pqr, v, exact, ex, tot, err, abs(tot - ex) / abs(ex), -np.log10(abs(tot - ex) / abs(ex)), time.time() - t0), flush=True)
