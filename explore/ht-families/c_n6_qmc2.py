"""(c) n = 6, (H,T) = (1,1): QMC of the float brute-force density (own map, Dirichlet(1/2) sampling)
for the SP graphs with E >= 1 and for E = 0 graphs; prints K estimate, error, and int |g| (scale)."""
import sys, random, time
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import hf_brute as hb, hf_numeric as hn
import c_pi, c_pi2

gs = [Gx for Gx in eg.enumerate_tight(2, 6) if eg.is_sp(Gx)]
which = [int(x) for x in sys.argv[1].split(",")]
m = int(sys.argv[2]) if len(sys.argv) > 2 else 12
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 606)
for k in which:
    G, _ = htk.from_networkx(gs[k], None)
    red = spred.sp_reduce(G, prefer="large")
    plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
    L, p, _ = c_pi.loops_and_p(G.ne, plan)
    E = c_pi2.excess(G.ne, plan)
    Gm = hb.G(G.nv, G.edges); C = Gm.cycles()
    while True:
        lam = np.array([[float(rng.choice([-3, -2, -1, 1, 2, 3]))] for _ in range(G.nv - 1)])
        if lam.sum() != 0:
            break
    t0 = time.time()
    mean, err, absint = hn.qmc_integral_general(G.nv, G.edges, C, plan, lam, m=m, nrep=16, seed=900 + k)
    print("graph %2d E=%d p=%d degrees %s lambda %s: K/pi^3 = %+.6f +- %.6f ; K = %+.6g +- %.2g ; int|g| = %.4g ; |K|/int|g| = %.2e  (%.0fs)" % (
        k, E, p, sorted(G.degree(v) for v in range(G.nv)), [int(x) for x in lam[:, 0]], mean / np.pi ** 3, err / np.pi ** 3,
        mean, err, absint, abs(mean) / absint, time.time() - t0), flush=True)
