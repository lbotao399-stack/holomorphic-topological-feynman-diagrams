"""Independent numerical check (htk.oracle, double precision, no reductions) of the T = 1 theta
kernels: K(lam; 0) and one first-order insertion K[y^i_e] = int D(t) y^i_e."""
import math
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

import htk_path  # noqa: F401
import htk
from theta_kernels import ThetaKernel
from kernels_check import CASES

which = sys.argv[1:] or [c[0] for c in CASES]
for name, G, paths, H, zspec in CASES:
    if not any(w in name for w in which):
        continue
    TK = ThetaKernel(G, paths, H, name)
    TK.check_pointwise(npts=3, seed=1)
    rng = random.Random(21)
    lamn = [[rng.randint(-4, 4) for _ in range(H)] for _ in range(G.nv - 1)]
    sub = {TK.lam[v][i]: lamn[v][i] for v in range(G.nv - 1) for i in range(H)}
    D = TK.Dtheta().subs(sub)
    if D == 0:
        lamn[0] = [lamn[0][i] + 1 for i in range(H)]
        sub = {TK.lam[v][i]: lamn[v][i] for v in range(G.nv - 1) for i in range(H)}
        D = TK.Dtheta().subs(sub)
    exact0 = float(TK.s * D / 2 * TK.integrate_poly(sp.Integer(1))[1] * sp.pi)
    y = TK.currents()
    e_ins, i_ins = TK.paths[1][0], 0
    r1 = TK.integrate_poly(sp.expand(y[e_ins][i_ins].subs(sub)))
    exact1 = float(TK.s * D / 2 * r1.get(1, 0) * sp.pi)
    dim = G.ne - 1
    method = 'lattice' if dim <= 6 else 'qmc'
    for label, powers, exact in [("K(lam;0)", None, exact0), ("K[y^%d_e%d]" % (i_ins + 1, e_ins), {(e_ins, i_ins): 1}, exact1)]:
        t0 = time.time()
        val, err = htk.kernel_numeric(G, lamn, T=1, powers=powers, method=method)
        dig = -math.log10(abs(val - exact) / abs(exact)) if val != exact else 16.0
        print("%s, lam = %s: %s oracle(%s) = %.9g +- %.2g, exact = %.9g, relative difference %.1e (%.1f digits)  (%.0fs)"
              % (name, lamn, label, method, val, err, exact, abs(val - exact) / abs(exact), dig, time.time() - t0), flush=True)
