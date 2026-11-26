"""C13, C14: S3 forms_check.mc_theta(1,2,2,H=1): value, sign, and the heavy tail of the sampled integrand.
- S3 mc_theta is run as is (S3 code imported read-only, no bytecode written) for the printed (mean, err, closed),
  and with nsamp = 1 at several seeds to read its single-sample summand.
- Own replay of the same random stream with the own brute-force density (vectorized float, Omega evaluated on the
  simplex tangent vectors e_j - e_{N-1}).
- Tail and Monte Carlo statistics with an own random stream; the density is the pointwise-verified closed form
  f = -D_theta / (2 U^{5/2}) (exact equality with the brute force established in v1_pointwise.py).
"""
import sys, time, math
from fractions import Fraction as Fr
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
from vcore import theta_std, bfs_cycles, incidence_full, full_sources, D_theta, ghalf, sg, density_bf, density_closed
import random

t0 = time.time()
th = theta_std(1, 2, 2)
ne, nv = th.ne, th.nv
C3 = np.array(th.C3(), float)
Bf = incidence_full(nv, th.edges)
KCL = np.array([[Bf[e][v] for e in range(ne)] for v in range(nv) if v != th.Np], float)
Cb = np.array(bfs_cycles(nv, th.edges, th.Np), float)


def brute_float(t, lam):
    """D(t) for H = 1, T = 1, cycle matrix C3, from Omega on the vectors e_j - e_4 (j = 0..3)."""
    P = t.shape[0]
    H = 1
    A = np.empty((P, ne, ne))
    A[:, :nv - 1, :] = KCL[None]
    A[:, nv - 1:, :] = Cb[None] * t[:, None, :]
    rhs = np.zeros((P, ne, 1))
    rhs[:, :nv - 1, 0] = lam[None, :]
    y = np.linalg.solve(A, rhs)[:, :, 0]
    Vt = np.zeros((ne, ne - 1))
    for j in range(ne - 1):
        Vt[j, j], Vt[ne - 1, j] = 1, -1
    F = np.stack([C3[0][None, :] * y, C3[1][None, :] * y], 1)       # (P, 2, ne)
    FV = F @ Vt                                                       # (P, 2, 4)
    M = np.einsum("ae,pe,be->pab", C3, t, C3)
    Ud = M[:, 0, 0] * M[:, 1, 1] - M[:, 0, 1] ** 2
    Minv = np.empty_like(M)
    Minv[:, 0, 0], Minv[:, 1, 1], Minv[:, 0, 1], Minv[:, 1, 0] = M[:, 1, 1] / Ud, M[:, 0, 0] / Ud, -M[:, 0, 1] / Ud, -M[:, 0, 1] / Ud
    dMV = np.einsum("ae,be,el->abl", C3, C3, Vt)
    om = np.zeros(P)
    m = ne - 1
    for k in range(m):
        for l in range(k + 1, m):
            a12 = np.zeros(P)
            for a in range(2):
                for b in range(2):
                    a12 += Minv[:, a, b] * (dMV[0, a, k] * dMV[b, 1, l] - dMV[0, a, l] * dMV[b, 1, k])
            a12 *= -0.5
            rest = [x for x in range(m) if x not in (k, l)]
            om += (-1) ** (k + l + 3) * a12 * np.linalg.det(FV[:, :, rest])
    om /= Ud ** 1.5
    dd = np.linalg.det(np.concatenate([t[:, None, :], np.tile(Vt.T[None], (P, 1, 1))], 1))
    return om / dd


# ---------------------------------------------------------------- S3 mc_theta itself
sys.path.insert(0, "/Users/libotao/Desktop/研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/code")
import forms_check as S3
ts = time.time()
mean, err, closed = S3.mc_theta(1, 2, 2, 1, 200000)
print("S3 mc_theta(1,2,2,1,200000) [seed 1]: mean %.10f  err %.10f  closed %.10f   (%.1f s)" % (mean, err, closed, time.time() - ts))

# ---------------------------------------------------------------- seed-1 sources, exact value
rg = np.random.default_rng(1)
lam = rg.normal(size=(nv - 1, 1))[:, 0]
print("seed-1 sources (N, B1, C1):", lam)
lamF = full_sources(nv, th.Np, {v: [Fr(float(lam[v]))] for v in range(nv - 1)})
c1 = Fr(1, 3)
Kq = sg(1) * c1 * D_theta(th, lamF, 1)   # K[C3] / pi
claimed = Fr(-22026136797979191689044246579867, 243388915243820045087367015432192)
print("own exact K[C3]/pi =", Kq, " = claimed fraction:", Kq == claimed, "  K =", float(Kq) * math.pi,
      "  -(pi/3) lam_B lam_C =", -math.pi / 3 * lam[1] * lam[2])
print("mc_theta 'closed' / K[C3] =", closed / (float(Kq) * math.pi))

# ---------------------------------------------------------------- replay of the stream with the own brute force
tt = np.array([rg.dirichlet(np.ones(ne)) for _ in range(200000)])
fv = brute_float(tt, lam)
vol = 1 / math.factorial(ne - 1)
m_own = fv.mean() * vol
e_own = math.sqrt(max((fv ** 2).mean() - fv.mean() ** 2, 0) / len(fv)) * vol
print("own replay: mean %.10f  err %.10f   rel. diff of mean vs S3: %.2e ; err rel diff: %.2e" %
      (m_own, e_own, abs(m_own / mean - 1), abs(e_own / err - 1)))
print("sampled f: min %.4f  max %.4f ; (exact - mean)/err = %.2f" % (fv.min(), fv.max(), (float(Kq) * math.pi - m_own) / e_own))
# exact check of the float brute force against the exact brute force at 3 sampled points
rr = random.Random(5)
okp = True
for i in rr.sample(range(200000), 3):
    tq = [Fr(float(x)) for x in tt[i]]
    bf, U, _ = density_bf(th, tq, lamF, 1, 1, rng=rr)
    ex = float(bf) / math.sqrt(float(U))
    okp &= abs(fv[i] / ex - 1) < 1e-10
print("float brute force = exact brute force (rel < 1e-10) at 3 sampled points:", okp)

# single-sample summands of S3 mc_theta vs own brute force
okS = True
for seed in (2, 3, 4, 5, 6):
    m1, _, _ = S3.mc_theta(1, 2, 2, 1, nsamp=1, seed=seed)
    g = np.random.default_rng(seed)
    l1 = g.normal(size=(nv - 1, 1))[:, 0]
    t1 = g.dirichlet(np.ones(ne))[None, :]
    own = brute_float(t1, l1)[0] * vol
    okS &= abs(m1 / own - 1) < 1e-12
print("S3 mc_theta single-sample summand = own brute-force D(t)/4! at 5 seeds (rel < 1e-12):", okS)

# ---------------------------------------------------------------- tails (own random stream)
c = abs(float(Kq) * 3)        # |D_theta| = |lam_B lam_C| ; f ~ -(|D_theta|/2) rho^{-5/2}
c = abs(lam[1] * lam[2]) / 2
g = np.random.default_rng(20261120)
n = 4_000_000
tt2 = g.dirichlet(np.ones(ne), size=n)
R = {"A": tt2[:, 0], "B": tt2[:, 1] + tt2[:, 2], "C": tt2[:, 3] + tt2[:, 4]}
U = R["A"] * R["B"] + R["B"] * R["C"] + R["C"] * R["A"]
f = -(lam[1] * lam[2]) / 2 / U ** 2.5            # closed density (sign of D_theta: lam_B lam_C > 0 here)
print("own stream n=%d: f < 0 everywhere: %s" % (n, bool(np.all(f < 0))))
xs = [1e2, 1e3, 1e4, 1e5]
emp = [float(np.mean(np.abs(f) > x)) for x in xs]
pred = [8 * (c / x) ** 1.2 for x in xs]
print("tail P(|f|>x): " + "  ".join("x=%g: emp %.3e pred 8(c/x)^(6/5) %.3e" % (x, a, b) for x, a, b in zip(xs, emp, pred)))
sl = np.polyfit(np.log(xs[1:]), np.log(emp[1:]), 1)[0]
print("fitted log-log slope (x = 1e3..1e5): %.3f  (claim -6/5)" % sl)
# corner asymptotics: U / rho -> 1 at the two corners
rhoAC = tt2[:, 0] + tt2[:, 3] + tt2[:, 4]
rhoAB = tt2[:, 0] + tt2[:, 1] + tt2[:, 2]
for nm_, rr_ in (("{tA1,tC1,tC2 -> 0}", rhoAC), ("{tA1,tB1,tB2 -> 0}", rhoAB)):
    print("corner %s: max |U/rho - 1| on rho < r for r = 0.04, 0.02, 0.01: %s" % (nm_, ["%.3e (%d pts)" % (np.max(np.abs(U[rr_ < r] / rr_[rr_ < r] - 1)), np.sum(rr_ < r)) for r in (0.04, 0.02, 0.01)]))
print("P(rho_AC < r)/(4 r^3) at r = 0.02, 0.05: %.3f %.3f" % (np.mean(rhoAC < 0.02) / (4 * 0.02 ** 3), np.mean(rhoAC < 0.05) / (4 * 0.05 ** 3)))
# second moment of the sample at growing n (own stream)
sm = {}
for k in (10 ** 4, 10 ** 5, 10 ** 6, 4 * 10 ** 6):
    sm[k] = float(np.mean(f[:k] ** 2))
print("sample second moments <f^2> at n = 1e4, 1e5, 1e6, 4e6:", ["%.4g" % sm[k] for k in sm])
# Monte Carlo runs: ratio mean/exact (estimator of K = vol * mean f)
exact = float(Kq) * math.pi
out = {}
for nrun, runs in ((2000, 2000), (20000, 600), (200000, 120)):
    ratios, zs = [], []
    for _ in range(runs):
        tr = g.dirichlet(np.ones(ne), size=nrun)
        RA, RB, RC = tr[:, 0], tr[:, 1] + tr[:, 2], tr[:, 3] + tr[:, 4]
        fr = -(lam[1] * lam[2]) / 2 / (RA * RB + RB * RC + RC * RA) ** 2.5
        mn = fr.mean() * vol
        er = fr.std() / math.sqrt(nrun) * vol
        ratios.append(mn / exact)
        zs.append((exact - mn) / er)
    ratios, zs = np.array(ratios), np.array(zs)
    out[nrun] = float(np.median(ratios))
    print("n=%7d runs=%4d  median(mean/exact) %.4f  fraction |mean|<|exact| %.3f  fraction |exact-mean|/err > 3: %.3f"
          % (nrun, runs, np.median(ratios), np.mean(ratios < 1), np.mean(np.abs(zs) > 3)), flush=True)
ns = np.array(list(out))
sl2 = np.polyfit(np.log(ns), np.log(1 - np.array([out[k] for k in ns])), 1)[0]
print("slope of log(1 - median ratio) vs log n: %.3f  (claim -0.184, stable-law exponent -1/6)" % sl2)
print("time %.1f s" % (time.time() - t0))
