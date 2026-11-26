"""Driver: Monte Carlo (mc_chain) for the n = 4 non-SP product forests and the n = 4 tadpoles at
(H,T) = (1,1).  Sources: forest roots (2, -1, 3, 1) in root order, lambda_A = -2; tadpoles: lambda_q = 1.
Calibration: one SP forest with a star-triangle chain against its exact kernel.
Usage: python3 run_mc_chain.py Mper batches part nparts [names...]   Output: mc_chain_*.log"""
import math
import sys
import time

from forests import forests, forest_graph, tadpole_graph
from htk import catalog, sp_reduce, star
import mc_chain


def jobs():
    out = []
    src, lamA = [2, -1, 3, 1], -2
    for F in forests(4):
        G, roots = forest_graph(F)
        lam = [[0.0] for _ in range(4)] + [[float(lamA)]]
        for i, v in enumerate(roots):
            lam[v] = [float(src[i])]
        if G.name == "(A,(A,B)) (A,B) (A,B)":
            mv = star.SingleMove(G, star.find_single(G)[0][0])
            out.append((G, lam, [mv], True))
        if not catalog.is_series_parallel(G):
            out.append((G, lam, None, False))
    for (c, L) in ((5, 0), (4, 1), (3, 2)):
        G, q = tadpole_graph(c, L)
        lam = [[0.0] for _ in range(G.nv - 1)]
        lam[q] = [1.0]
        out.append((G, lam, None, False))
    return out


if __name__ == "__main__":
    Mper, batches = int(sys.argv[1]), int(sys.argv[2])
    part, nparts = int(sys.argv[3]), int(sys.argv[4])
    names = sys.argv[5:]
    for idx, (G, lam, chain, calib) in enumerate(jobs()):
        if idx % nparts != part or (names and G.name not in names):
            continue
        t0 = time.time()
        val, err, ch = mc_chain.estimate(G, lam, Mper, batches, seed=idx, chain=chain)
        ex = ""
        if calib:
            K = sp_reduce(G, prefer="large").kernel([[int(x[0])] for x in lam], T=1)
            ex = "  exact K/pi^2 = %.6f" % (float(K) / math.pi ** 2)
        print("%-28s chart %-26s K/pi^2 = %+.5f +- %.5f%s  (%d pts, %.0fs)" % (
            G.name, mc_chain.describe(ch), val / math.pi ** 2, err / math.pi ** 2, ex, Mper * batches, time.time() - t0), flush=True)
