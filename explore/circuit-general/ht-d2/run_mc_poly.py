"""Monte Carlo (mc_chain charts) of the source polynomial of the single-tree non-SP n = 4 product
forests at (H,T) = (1,1): K = sum_{k=0..4} c_k lambda_0^k lambda_A^(4-k) (lambda_B = -lambda_0 - lambda_A),
from 5 source vectors evaluated on common samples.  Output: mc_poly_*.log
Usage: python3 run_mc_poly.py Mper batches part nparts"""
import math
import sys
import time

import numpy as np

from forests import forests, forest_graph
from htk import catalog
import mc_chain

PTS = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2)]


def run(G, roots, Mper, batches, seed):
    chain = mc_chain.find_chain(G)
    last = chain[-1].red if chain else G
    spr = mc_chain.sp_reduce(last, prefer="large")
    rng = np.random.default_rng(seed)
    V = np.array([[l0 ** k * lA ** (4 - k) for k in range(5)] for (l0, lA) in PTS], dtype=float)
    Vinv = np.linalg.inv(V)
    coef_batches = []
    for _ in range(batches):
        bary, logq = mc_chain._sample_spred(spr, rng, Mper)
        T, J = mc_chain._spred_times(spr, last, bary, Mper)
        for mv in reversed(chain):
            f, lq = mc_chain._sample_fractions(mv, rng, Mper)
            logq = logq + lq
            T, Jl = (mc_chain._single_map if isinstance(mv, mc_chain.star.SingleMove) else mc_chain._paired_map)(mv, f, T)
            J = J * Jl
        vals = []
        for (l0, lA) in PTS:
            lam = [[0.0] for _ in range(G.nv - 1)]
            lam[roots[0]] = [float(l0)]
            lam[G.nv - 2] = [float(lA)]          # A = n
            D = mc_chain._safe_density(G, T, lam, None)
            w = D * J * np.exp(-logq)
            w[~np.isfinite(w)] = 0.0
            vals.append(np.mean(w))
        coef_batches.append(Vinv @ np.array(vals))
    C = np.array(coef_batches)
    return C.mean(axis=0), C.std(axis=0, ddof=1) / math.sqrt(len(C)), chain


if __name__ == "__main__":
    Mper, batches, part, nparts = (int(x) for x in sys.argv[1:5])
    jobs = []
    for F in forests(4):
        G, roots = forest_graph(F)
        if len(F) == 1 and (not catalog.is_series_parallel(G) or G.name == "(A,(A,(A,(A,B))))"):
            jobs.append((G, roots))   # (A,(A,(A,(A,B)))): SP calibration, exact -pi^2 l0^2 (l0+lA)^2/180
    for idx, (G, roots) in enumerate(jobs):
        if idx % nparts != part:
            continue
        t0 = time.time()
        m, e, ch = run(G, roots, Mper, batches, seed=100 + idx)
        s = "  ".join("c%d/pi^2 = %+.5f +- %.5f" % (k, m[k] / math.pi ** 2, e[k] / math.pi ** 2) for k in range(5))
        print("%-22s chart %-24s K = sum_k c_k l0^k lA^(4-k):  %s  (%.0fs)" % (G.name, mc_chain.describe(ch), s, time.time() - t0), flush=True)
