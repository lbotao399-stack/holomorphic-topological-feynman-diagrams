"""Zero-displacement values: (a) exact, from the bracket expansion of kernels_q3.json;
(b) htk.kernel_exact (separate implementation of the series-parallel layer; I^S2 = K for n = 3);
(c) direct numerical integral over the Schwinger simplex of the brute-force density
    eps_T det(B_T) (-1)^{N-1} det[d y^a_c / d t_j]_{j < N-1} / t_{N-1}   (mesh currents, numpy),
    scrambled Sobol points with the smoothing u -> 3u^2 - 2u^3 and stick-breaking (own code).
Sources: random integers at the five output vertices, zero at the trivalent vertices."""
import json
import math
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
from scipy.stats import qmc

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from q3lib import GRAPHS, Reduction, Brute, ROOT, fdet

KJ = json.load(open('kernels_q3.json'))
rng = np.random.RandomState(7)
NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 18
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 8
import os
LOG = os.environ.get('LOGQMC', '1') == '1'
BSC = float(os.environ.get('BSC', '1.5'))
KOR = int(os.environ.get('KOR', '5'))
names = sys.argv[3:] or list(GRAPHS)


def exact_K0(name, lam):
    tot = Fr(0)
    for b1, b2, b3, c in KJ[name]['K0']:
        val = Fr(sp.Rational(c).p, sp.Rational(c).q)
        for b in (b1, b2, b3):
            val *= fdet([[lam[v][a] for v in b] for a in range(3)])
        tot += val
    return tot


def numeric(name, lam, npts, nrep):
    red = Reduction(name)
    E = red.edges
    N = len(E)
    B = Brute(E, red.cotree())
    sgn = int(B.epsT * B.detBT)
    C = np.array(B.C, dtype=float)                   # n x N
    tree = B.tree
    BTinv = np.array([[float(x) for x in row] for row in B.BTinv])
    lamf = np.array([[float(x) for x in v] for v in lam])   # 7 x 3
    y0 = np.zeros((N, 3))
    for j, e in enumerate(tree):
        y0[e] = BTinv[j] @ lamf
    cot = B.cot
    dim = N - 1

    def dens(t):
        # t: (b, N)
        M = np.einsum('re,be,se->brs', C, t, C)
        ok = np.abs(np.linalg.det(M)) > 1e-300
        M[~ok] = np.eye(M.shape[1])
        Minv = np.linalg.inv(M)
        rhs = -np.einsum('re,be,ea->bra', C, t, y0)
        l = np.einsum('brs,bsa->bra', Minv, rhs)
        y = y0[None] + np.einsum('re,bra->bea', C, l)
        # d l_r^a / d t_j = -(M^{-1} C_{:,j})_r y_j^a
        MC = np.einsum('brs,sj->brj', Minv, C)                 # b x n x N
        D = -np.einsum('brj,bja->braj', MC, y)                   # b x n x 3 x N
        J = D[:, :, :, :dim].reshape(t.shape[0], 3 * len(cot), dim)
        return sgn * (-1) ** (N - 1) * np.linalg.det(J) / t[:, N - 1]

    if LOG:
        # log coordinates t_j = exp(s_j) (j < N-1), t_{N-1} = 1; s ~ Laplace(0, b) via Sobol points
        vals = []
        for rep in range(nrep):
            u = qmc.Sobol(dim, scramble=True, seed=100 + rep).random(npts)
            tot = 0.0
            for a in range(0, npts, 65536):
                uu = np.clip(u[a:a + 65536], 1e-15, 1 - 1e-15)
                sv = np.where(uu < 0.5, BSC * np.log(2 * uu), -BSC * np.log(2 * (1 - uu)))
                logp = np.sum(-np.abs(sv) / BSC - np.log(2 * BSC), axis=1)
                sfull = np.concatenate([sv, np.zeros((uu.shape[0], 1))], axis=1)
                sfull = sfull - sfull.mean(axis=1, keepdims=True)
                t = np.exp(sfull)
                with np.errstate(all='ignore'):
                    d = dens(t) * np.prod(t, axis=1) * np.exp(-logp)   # S(t) = D(t) prod t_e, scale invariant
                d[~np.isfinite(d)] = 0.0
                tot += float(np.sum(d))
            vals.append(tot / npts)
        vals = np.array(vals)
        return vals.mean(), vals.std(ddof=1) / math.sqrt(nrep)
    vals = []
    for rep in range(nrep):
        u = qmc.Sobol(dim, scramble=True, seed=100 + rep).random(npts)
        if KOR == 3:
            v = u * u * (3.0 - 2.0 * u)
            jac = np.prod(6.0 * u * (1.0 - u), axis=1)
        else:
            v = u ** 3 * (10.0 - 15.0 * u + 6.0 * u * u)
            jac = np.prod(30.0 * u * u * (1.0 - u) ** 2, axis=1)
        tot = 0.0
        for a in range(0, npts, 65536):
            vv = v[a:a + 65536]
            jj = jac[a:a + 65536].copy()
            b = vv.shape[0]
            t = np.empty((b, N))
            rem = np.ones(b)
            for k in range(dim):
                t[:, k] = rem * vv[:, k]
                jj *= rem
                rem = rem * (1.0 - vv[:, k])
            t[:, dim] = rem
            with np.errstate(all='ignore'):
                d = dens(t)
            d[~np.isfinite(d)] = 0.0
            tot += float(np.sum(jj * d))
        vals.append(tot / npts)
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / math.sqrt(nrep)


for name in names:
    outs = KJ[name]['outputs']
    lam = [[Fr(0)] * 3 for _ in range(ROOT)]
    for v in outs:
        lam[v] = [Fr(int(x)) for x in rng.randint(-4, 5, size=3)]
    ex = exact_K0(name, lam)
    t0 = time.time()
    G = htk.Graph(8, GRAPHS[name]['edges'])
    kh = htk.kernel_exact(G, [[int(x) for x in v] for v in lam], T=0, check=2)
    th = time.time() - t0
    t0 = time.time()
    val, err = numeric(name, lam, NPTS, NREP)
    tn = time.time() - t0
    rel = abs(val - float(ex)) / abs(float(ex))
    digits = -math.log10(rel) if rel > 0 else 16
    print('%s: lam(outputs) = %s\n    exact (bracket expansion) = %s = %.10g\n    htk.kernel_exact          = %s  equal: %s  [%.1fs]\n'
          '    Sobol, %d x %d points    = %.6g +- %.2g   rel. dev. %.1e (%.1f digits)  [%.0fs]'
          % (name, [[int(x) for x in lam[v]] for v in outs], ex, float(ex), kh, Fr(kh) == ex, th,
             NREP, NPTS, val, err, rel, digits, tn), flush=True)
