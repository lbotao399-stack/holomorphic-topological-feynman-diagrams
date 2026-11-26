"""Re-run of the own QMC (V2 of v_b1) for the two cases with |diff|/err > 2 (G3 and G4, first source set),
2^16 x 16 points, three independent seed blocks each."""
import sys, time
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
from v_b1 import GR, FORM, L, own_qmc
import sympy as sp
for k, lam in [(4, [2, -5, 5, 2, -7]), (3, [-6, -6, -7, -5, 2])]:
    ex = float(FORM[k].subs({s: v for s, v in zip(L, lam)}))
    for seed in (1000, 2000, 3000):
        t0 = time.time()
        mean, err, nex = own_qmc(6, GR[k], np.array([[float(v)] for v in lam]), 16, 16, seed)
        print("G%d lambda %s seed %d: own QMC 2^16 x 16 = %.8g +- %.2g ; claim %.10g ; diff/err %.2f ; rel diff %.1e ; %d exact (%.1fs)" % (
            k, lam, seed, mean, err, ex, (mean - ex) / err, abs(mean - ex) / abs(ex), nex, time.time() - t0))
