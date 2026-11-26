"""Float Monte Carlo (independent oracle, not exact): T = 1, H = 1 kernels of cubic product forests in
the chart  (single star-triangle move at a leaf, htk.star.SingleMove)  x  (series-parallel chart of the
reduced graph, htk.spred), with the brute-force float density htk.oracle.density_batch.

    K = int D_Gamma(t(f, T(xi))) |J_aff(f, T) det[T; dT/dxi]| df dxi,
sampled with f_i ~ Beta(1/2, 1/2) and Dirichlet(1/2, ..., 1/2) on every simplex of the reduced chart.
Usage: python3 nonsp_numeric.py N_per_batch batches [forest names...]   Output: nonsp_numeric.log
"""
import math
import sys
import time

import numpy as np

from forests import forests, forest_graph
from htk import star, catalog, sp_reduce, oracle
from htk.linalg import perm_sign


def chart(G):
    """First single-move site whose reduced graph is series-parallel."""
    for (v, u, w, x) in star.find_single(G):
        mv = star.SingleMove(G, v)
        if catalog.is_series_parallel(mv.red):
            return mv
    return None


def sampler(mv, spr, rng, M):
    """Samples: r (M x 2), bary dict name -> (M,) and the log proposal density."""
    logq = np.zeros(M)
    r = rng.beta(0.5, 0.5, size=(M, 2))
    logq += np.sum(-0.5 * np.log(r) - 0.5 * np.log1p(-r) - math.log(math.pi), axis=1)
    bary = {}
    for blk in spr.blocks:
        m = len(blk.bary)
        x = rng.dirichlet([0.5] * m, size=M)
        for j, b in enumerate(blk.bary):
            bary[b] = x[:, j]
        # Dirichlet(1/2) density on the simplex w.r.t. the independent coordinates
        logq += math.lgamma(0.5 * m) - m * math.lgamma(0.5) + np.sum(-0.5 * np.log(x), axis=1)
    return r, bary, logq


def integrand(G, mv, spr, lam, r, bary):
    M = r.shape[0]
    Nr = mv.red.ne
    # reduced times
    T = np.ones((M, Nr))
    for e in range(Nr):
        for b, k in spr.texp[e].items():
            T[:, e] *= bary[b] ** k
    # J_red = det[T ; dT/dxi] = prod T * det[1 ; dlogT/dxi]
    cols = []
    for c in spr.coords:
        blk = next(bb for bb in spr.blocks if c in bb.bary)
        b0 = blk.bary[0]
        col = np.zeros((M, Nr))
        for e in range(Nr):
            k1 = spr.texp[e].get(c, 0)
            k0 = spr.texp[e].get(b0, 0)
            col[:, e] = k1 / bary[c] - k0 / bary[b0]
        cols.append(col)
    Lmat = np.stack([np.ones((M, Nr))] + cols, axis=1)        # M x Nr x Nr (rows: 1, coords)
    Jred = np.prod(T, axis=1) * np.linalg.det(Lmat)
    # star map
    r1, r2 = r[:, 0], r[:, 1]
    Tuw, Tux, Twx = T[:, 0], T[:, 1], T[:, 2]
    Rux, Rwx = Tux / (1 - r1), Twx / (1 - r2)
    P = Tuw + Rux + Rwx
    t = np.zeros((M, G.ne))
    t[:, mv.e["vu"]] = Tuw * Rux / P
    t[:, mv.e["vw"]] = Tuw * Rwx / P
    t[:, mv.e["vx"]] = Rux * Rwx / P
    t[:, mv.e["ux"]] = Tux / r1
    t[:, mv.e["wx"]] = Twx / r2
    for j, e in enumerate(mv.others):
        t[:, e] = T[:, 3 + j]
    order_rows = [mv.e["ux"], mv.e["wx"], mv.e["vu"], mv.e["vw"], mv.e["vx"]] + mv.others
    sgn = perm_sign([order_rows.index(e) for e in range(G.ne)])
    Jaff = sgn * (-Tuw * Tux ** 2 * Twx ** 2 / (r1 ** 2 * r2 ** 2 * (1 - r1) ** 3 * (1 - r2) ** 3 * P ** 3))
    D = oracle.density_batch(G, t, lam, T=1)
    return D * np.abs(Jaff * Jred)


def estimate(G, lam, Mper, batches, seed=0):
    mv = chart(G)
    spr = sp_reduce(mv.red, prefer="large")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(batches):
        r, bary, logq = sampler(mv, spr, rng, Mper)
        g = integrand(G, mv, spr, lam, r, bary)
        w = g * np.exp(-logq)
        w = w[np.isfinite(w)]
        vals.append(np.mean(w))
    vals = np.array(vals)
    return float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(len(vals))), mv


if __name__ == "__main__":
    Mper = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    batches = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    names = sys.argv[3:]
    n = 4
    src = [2, -1, 3, 1]
    lamA = -2
    for F in forests(n):
        G, roots = forest_graph(F)
        if names and G.name not in names:
            continue
        mv = chart(G)
        if mv is None:
            print("%-28s no single-move chart with SP reduced graph" % G.name, flush=True)
            continue
        lam = [[0.0] for _ in range(n)] + [[float(lamA)]]
        for i, v in enumerate(roots):
            lam[v] = [float(src[i])]
        t0 = time.time()
        val, err, mv = estimate(G, lam, Mper, batches)
        ex = ""
        if catalog.is_series_parallel(G):
            K = sp_reduce(G, prefer="large").kernel([[int(x[0])] for x in lam], T=1)
            ex = " exact %.6f (%s)" % (float(K), K)
        print("%-28s centre %d  K/pi^2 = %.5f +- %.5f%s  (%.0fs)" % (G.name, mv.v, val / math.pi ** 2, err / math.pi ** 2,
                                                              (" exact/pi^2 %.5f" % (float(K) / math.pi ** 2)) if ex else "", time.time() - t0), flush=True)
