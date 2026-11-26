"""Numerical integral of the brute-force density pulled back to the ratio coordinates (floating
point): g(r) = eps_T det(B_T) det[d y^a_c/d r_k] sgn(omega_m(d_r)), d y/d r = (d y/d t)(d t/d r) with
mesh currents at t(r) (no polynomial formula, no Dirichlet integrals).  Domain: (0,1) for the r_1
of chain/bubble moves and the bubble ratios, ordered simplices sampled by sorting; Sobol points."""
import json
import math
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
from scipy.stats import qmc

from q3lib import GRAPHS, Reduction, Brute, ROOT, R, fdet

KJ = json.load(open('kernels_q3.json'))
rng = np.random.RandomState(7)          # same sources as numerics_q3.py
NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 18
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def exact_K0(name, lam):
    tot = Fr(0)
    for b1, b2, b3, c in KJ[name]['K0']:
        val = Fr(sp.Rational(c).p, sp.Rational(c).q)
        for b in (b1, b2, b3):
            val *= fdet([[lam[v][a] for v in b] for a in range(3)])
        tot += val
    return tot


for name in GRAPHS:
    outs = KJ[name]['outputs']
    lam = [[Fr(0)] * 3 for _ in range(ROOT)]
    for v in outs:
        lam[v] = [Fr(int(x)) for x in rng.randint(-4, 5, size=3)]
    ex = float(exact_K0(name, lam))
    t0 = time.time()
    red = Reduction(name)
    N = len(red.edges)
    B = Brute(red.edges, red.cotree())
    sgn = int(B.epsT * B.detBT)
    C = np.array(B.C, dtype=float)
    y0 = np.zeros((N, 3))
    lamf = np.array([[float(x) for x in v] for v in lam])
    BTinv = np.array([[float(x) for x in row] for row in B.BTinv])
    for j, e in enumerate(B.tree):
        y0[e] = BTinv[j] @ lamf
    tf = sp.lambdify(R, red.t, 'numpy')
    dtf = sp.lambdify(R, [[sp.diff(te, r) for r in R] for te in red.t], 'numpy')
    dom = red.domain()
    vals = []
    for rep in range(NREP):
        u = qmc.Sobol(9, scramble=True, seed=300 + rep).random(NPTS)
        cols, vol, k = {}, 1.0, 0
        for kind, rs in dom:
            if kind == 'cube':
                cols[rs[0]] = u[:, k]; k += 1
            else:
                m = len(rs)
                srt = -np.sort(-u[:, k:k + m], axis=1)
                for i, r in enumerate(rs):
                    cols[r] = srt[:, i]
                vol /= math.factorial(m); k += m
        args = [cols[r] for r in R]
        b = NPTS
        t = np.stack([np.broadcast_to(np.asarray(x, dtype=float), (b,)) for x in tf(*args)], axis=1)
        dt = np.array([[np.broadcast_to(np.asarray(x, dtype=float), (b,)) for x in row] for row in dtf(*args)])  # N x 9 x b
        dt = np.transpose(dt, (2, 0, 1))                                   # b x N x 9
        M = np.einsum('re,be,se->brs', C, t, C)
        Minv = np.linalg.inv(M)
        rhs = -np.einsum('re,be,ea->bra', C, t, y0)
        l = np.einsum('brs,bsa->bra', Minv, rhs)
        y = y0[None] + np.einsum('re,bra->bea', C, l)
        MC = np.einsum('brs,sj->brj', Minv, C)
        Dyt = -np.einsum('brj,bja->braj', MC, y)                            # b x n x 3 x N
        Dyr = np.einsum('braj,bjk->brak', Dyt, dt).reshape(b, 9, 9)
        om = np.linalg.det(np.concatenate([t[:, None, :], np.transpose(dt, (0, 2, 1))], axis=1))
        g = sgn * np.linalg.det(Dyr) * np.sign(om)
        vals.append(vol * g.mean())
    vals = np.array(vals)
    val, err = vals.mean(), vals.std(ddof=1) / math.sqrt(NREP)
    rel = abs(val - ex) / abs(ex)
    print('%s: exact %.10g, brute-force density in ratio coordinates (Sobol %d x %d) %.10g +- %.1e, '
          'rel. dev. %.1e (%.1f digits)  [%.0fs]' % (name, ex, NREP, NPTS, val, err, rel, -math.log10(rel), time.time() - t0),
          flush=True)
