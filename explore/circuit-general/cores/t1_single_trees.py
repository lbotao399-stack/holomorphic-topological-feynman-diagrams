"""Exact (1,1) source polynomials of the n = 4 single product trees (ht-d2 notation):
K/pi^2 = sum_k c_k l0^k lA^(4-k), from kernels at lA = 1, l0 = 0..4 (all other sources 0)."""
import importlib.util
import os
import sys
import time
from fractions import Fraction as F

import cores

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('htd2_forests', os.path.join(HERE, '..', 'ht-d2', 'forests.py'))
fm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fm)

names = sys.argv[1:]
n = 4
for Fo in fm.forests(n):
    G, roots = fm.forest_graph(Fo)
    if len(roots) != 1 or (names and G.name not in names):
        continue
    t = time.time()
    vals = []
    for l0 in range(5):
        lam = [[0] for _ in range(n)] + [[1]]
        lam[roots[0]] = [l0]
        v = cores.kernel(G.nv, G.edges, lam, (), T=1)
        vals.append(v)
    # K(l0) = sum_k c_k l0^k  (lA = 1) per label; Lagrange on l0 = 0..4
    labels = sorted({k for v in vals for k in v}, key=str)
    res = {}
    for lab in labels:
        ys = [v.get(lab, F(0)) for v in vals]
        # Newton -> monomial coefficients
        import itertools
        coef = [F(0)] * 5
        for i in range(5):
            # Lagrange basis L_i(x) = prod_{j != i} (x - j)/(i - j)
            poly = [F(1)]
            den = F(1)
            for j in range(5):
                if j == i:
                    continue
                poly = [F(0)] + poly
                for k in range(len(poly) - 1):
                    poly[k] -= j * poly[k + 1]
                den *= (i - j)
            for k in range(5):
                coef[k] += ys[i] * poly[k] / den
        res[cores.enc(lab)] = [str(c) for c in coef]
    print('%-22s c_0..c_4 per label: %s  (%.0fs)' % (G.name, res, time.time() - t), flush=True)
