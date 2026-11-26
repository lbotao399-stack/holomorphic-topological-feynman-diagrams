"""(d) graph 2: QMC of the float brute-force density, with int |g| for scale (Dirichlet(1/2), own map)."""
import sys
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import hf_brute as hb, hf_numeric as hn
edges = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)]
G = htk.Graph(6, edges)
for prefer in ("small", "large"):
    red = spred.sp_reduce(G, prefer=prefer)
    plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
    C = hb.G(6, edges).cycles()
    for lam in ([[5.], [-5.], [-4.], [-5.], [5.]], [[1.], [-2.], [4.], [4.], [4.]]):
        mean, err, absint = hn.qmc_integral_general(6, edges, C, plan, np.array(lam), m=14, nrep=16, seed=4242)
        print("plan %s lambda %s: K = %+.5g +- %.2g ; int|g| = %.5g ; |K|/int|g| = %.1e" % (
            prefer, [int(x[0]) for x in lam], mean, err, absint, abs(mean) / absint), flush=True)
