"""Numerical integral of the brute-force T = 1 density in the chart t^P_j = x^P_j / w_P of RP^{N-1}_{>0}
(w in Delta^2, x^P in Delta^{n_P-1}): the integrand is D(t) |det[t ; dt/dxi]| with D(t) from the
vectorized Kirchhoff solve and forms (numeric_verify.density), not the claimed density.
Importance sampling: w ~ Dirichlet(1/2,1/2,1/2) (stick breaking with Beta inverse CDFs), x^P uniform;
scrambled Sobol points, R scrambles.  Compared with the exact K(lam;0) and int D y^0_{e} (e = spoke of B).
Usage: python3 numeric_wx.py case m R
"""
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from scipy.stats import qmc, beta
from scipy.special import gammaln
import sympy as sp

case = sys.argv[1]
m = int(sys.argv[2]) if len(sys.argv) > 2 else 20
R = int(sys.argv[3]) if len(sys.argv) > 3 else 8
sys.argv = sys.argv[:2]
import numeric_verify as nv      # builds TH, lam, exact values (with a tiny QMC run that is ignored)

TH, net, N = nv.TH, nv.net, nv.N
nP = TH.nP
paths = nv.paths
d = 2 + sum(n - 1 for n in nP)


def sample(u):
    P = u.shape[0]
    # w ~ Dirichlet(1/2, 1/2, 1/2): w_A = Beta(1/2, 1), w_B = (1 - w_A) Beta(1/2, 1/2)
    wA = beta.ppf(u[:, 0], 0.5, 1.0)
    r = beta.ppf(u[:, 1], 0.5, 0.5)
    w = np.stack([wA, (1 - wA) * r, (1 - wA) * (1 - r)], axis=1)
    logpdf = gammaln(1.5) - 3 * gammaln(0.5) - 0.5 * np.sum(np.log(w), axis=1)
    xs, k = [], 2
    for n in nP:
        if n == 1:
            xs.append(np.ones((P, 1)))
            continue
        # uniform on Delta^{n-1}: spacings of sorted uniforms; density (n-1)!
        s = np.sort(u[:, k:k + n - 1], axis=1)
        s = np.concatenate([np.zeros((P, 1)), s, np.ones((P, 1))], axis=1)
        xs.append(np.diff(s, axis=1))
        logpdf += gammaln(n)
        k += n - 1
    return w, xs, logpdf


def tmap(w, xs):
    P = w.shape[0]
    t = np.empty((P, N))
    for Pi in range(3):
        for j, e in enumerate(paths[Pi]):
            t[:, e] = xs[Pi][:, j] / w[:, Pi]
    return t


def jac(w, xs):
    P = w.shape[0]
    rows = [tmap(w, xs)]
    for Q in (1, 2):
        row = np.zeros((P, N))
        for Pi in range(3):
            dw = (1.0 if Pi == Q else 0.0) - (1.0 if Pi == 0 else 0.0)
            for j, e in enumerate(paths[Pi]):
                row[:, e] = -xs[Pi][:, j] / w[:, Pi] ** 2 * dw
        rows.append(row)
    for Pi in range(3):
        for l in range(1, nP[Pi]):
            row = np.zeros((P, N))
            for j, e in enumerate(paths[Pi]):
                row[:, e] = ((j == l) - (j == 0)) / w[:, Pi]
            rows.append(row)
    Mt = np.stack(rows, axis=1)
    return np.abs(np.linalg.det(Mt))


t0 = time.time()
v0, v1 = [], []
for rep in range(R):
    sob = qmc.Sobol(d, scramble=True, seed=500 + rep)
    a0 = a1 = 0.0
    nb = 0
    for ch in range(1 << max(0, m - 17)):
        u = np.clip(sob.random(1 << min(m, 17)), 1e-14, 1 - 1e-14)
        w, xs, logpdf = sample(u)
        t = tmap(w, xs)
        Dv, y = nv.density(t)
        g = Dv * jac(w, xs) / np.exp(logpdf)
        a0 += g.sum(); a1 += (g * y[:, nv.eB, 0]).sum(); nb += len(g)
    v0.append(a0 / nb); v1.append(a1 / nb)
v0, v1 = np.array(v0), np.array(v1)
print("case %s (chart (w, x), brute-force integrand), lam = %s:" % (case, nv.lam_int))
for nm, v, ex in [("K(lam;0)", v0, nv.K0), ("int D y^0_{e=%d}" % nv.eB, v1, nv.Ins)]:
    mean, err = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
    rel = abs(mean - ex) / abs(ex)
    print("   %-18s QMC %.10g +- %.2g (2^%d x %d)   exact %.10g   rel. diff %.1e (%.1f digits), |diff|/err = %.1f"
          % (nm, mean, err, m, R, ex, rel, -np.log10(rel), abs(mean - ex) / err))
print("   time %.0fs" % (time.time() - t0))
