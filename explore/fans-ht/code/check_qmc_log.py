"""N3': n = 4, chart-free: K = int_{R^8} S(e^s) d^8 s (s_9 = 0), S = D prod t_e (htk.oracle.scaled_density_batch,
brute force).  s_j = c logit(x_j), c = 3, x scrambled Sobol; points with max|s_j| > smax are skipped (their
count and the largest |integrand| among kept points with max|s| > smax - 3 are printed)."""
import sys, time, math
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from htk import fan, oracle
from scipy.stats import qmc
c, smax, npts, nrep = 3.0, 26.0, 2 ** 17, 8
G = fan(4)
st = oracle._Setup(G)
exact = {(1, 0, 0, 0, 1): math.pi ** 2 / 180,
         (2, -1, 3, 1, -2): None}
l1, l2, l3, l4, l5 = 2, -1, 3, 1, -2
exact[(2, -1, 3, 1, -2)] = math.pi ** 2 * l1 * l5 * ((l1 + 2*l2 + 2*l3) * (2*l3 + 2*l4 + l5) + 8*l3**2) / 180
for lam, ex in exact.items():
    t0 = time.time()
    lamf = np.array([[float(x)] for x in lam])
    vals, skipped, edge = [], 0, 0.0
    for rep in range(nrep):
        x = qmc.Sobol(8, scramble=True, seed=500 + rep).random(npts)
        x = np.clip(x, 1e-300, 1 - 1e-16)
        s = c * np.log(x / (1 - x))
        jac = np.prod(c / (x * (1 - x)), axis=1)
        ok = np.abs(s).max(axis=1) <= smax
        skipped += int((~ok).sum())
        t = np.concatenate([np.exp(s[ok]), np.ones((ok.sum(), 1))], axis=1)
        S = oracle.scaled_density_batch(G, t, lamf, 1, None, None, None, st)
        f = S * jac[ok]
        near = np.abs(s[ok]).max(axis=1) > smax - 3
        if near.any():
            edge = max(edge, float(np.abs(f[near]).max()) / npts)
        vals.append(float(f.sum()) / npts)
    vals = np.array(vals)
    m, e = vals.mean(), vals.std(ddof=1) / math.sqrt(nrep)
    print("N3' n=4 lam=%s: log-space QMC %.6f +- %.1e (skipped %d of %d points; max weighted integrand near the cut %.1e), exact %.6f, rel.diff %.1e (%.0f s)" % (
        list(lam), m, e, skipped, nrep * npts, edge, ex, abs(m - ex) / abs(ex), time.time() - t0), flush=True)
