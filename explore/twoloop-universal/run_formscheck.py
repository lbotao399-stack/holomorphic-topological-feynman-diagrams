"""Item (e): S3 code/forms_check.py, mc_theta(1, 2, 2, H=1, nsamp=200000, seed=1).

1. exact value of the quantity estimated by mc_theta (integrand U^{-H-1/2} <psi^tau_1 psi^tau_2> ^ prod psi,
   cycle matrix (A-B, B-C), orientation iota_E vol = K[C3] of twoloop.py) for the sources of seed 1;
2. pointwise identity of the mc_theta integrand (re-implemented, float) with the brute force and the
   closed density at random points;
3. deterministic integrals: own chart Gauss-Legendre, htk oracle lattice;
4. the Monte Carlo estimator itself: vectorized replay of the same random stream (same mean as
   mc_theta), tail of the sampled values, distribution of the estimator over independent runs.
Writes formscheck_results.json."""
import json
import math
import sys
import time
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import (theta_graph, Theta, full_sources, kernel_C3_closed, brute_density, closed_density_t,
                     D_theta)  # noqa
import numint

S3CODE = "/Users/libotao/Desktop/研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/code"

out = {}
pqr, H, T = (1, 2, 2), 1, 1
nv, edges = theta_graph(*pqr)
th = Theta(nv, edges)
st = numint.Setup(pqr)

# ---------------------------------------------------------------- 1. sources of mc_theta(seed=1), exact value
rg = np.random.default_rng(1)
lam = rg.normal(size=(nv - 1, H))
lamF = [[Fr(float(x)) for x in row] for row in lam]          # exact binary rationals of the doubles
lf = full_sources(nv, nv - 1, {v: lamF[v] for v in range(nv - 1)})
KC3 = kernel_C3_closed(th, lf, H, T)                          # coefficient of pi
KC3f = float(KC3) * math.pi
closed_printed = math.gamma(2.5 - 1) * math.gamma(2.5 - 2) ** 2 / (2 * math.gamma(2.5)) * float(D_theta(th, lf))
out["sources_seed1"] = lam.ravel().tolist()
out["K_C3_over_pi_exact"] = str(KC3)
out["K_C3"] = KC3f
out["minus_pi_over_3_lamB_lamC"] = -math.pi / 3 * lam[1, 0] * lam[2, 0]
out["closed_as_printed_by_mc_theta"] = closed_printed
print("sources (N, B1, C1):", lam.ravel())
print("K[C3] = %s * pi = %.10f ;  -(pi/3) lam_B lam_C = %.10f ; mc_theta 'closed' = %.10f" % (
    KC3, KC3f, out["minus_pi_over_3_lamB_lamC"], closed_printed))

# ---------------------------------------------------------------- 2. pointwise: mc_theta integrand
sys.dont_write_bytecode = True      # no __pycache__ in the source folder
sys.path.insert(0, S3CODE)
import forms_check as fc  # noqa: E402  (read-only import; module code under __main__ is not executed)


def mc_integrand(t):
    """the summand 'val' of forms_check.mc_theta at one point t (copied loop structure, float)"""
    p, q, r = pqr
    V, edg, _ = fc.theta(*pqr)
    E = len(edg)
    C = np.zeros((2, E)); C[0, :p] = 1; C[0, p:p + q] = -1; C[1, p:p + q] = 1; C[1, p + q:] = -1
    e0 = E - 1
    rem = [c for c in range(E) if c != e0]
    B = fc.incidence(V, edg)
    Lap = B.T @ (B / t[:, None])
    y = (B @ np.linalg.solve(Lap, lam)) / t[:, None]
    Mper = C @ (C * t).T; U = np.linalg.det(Mper)
    G = 0.5 * (np.diag(1 / t) - C.T @ np.linalg.inv(Mper) @ C)
    M = np.array([C[rr] * y[:, i] for i in range(H) for rr in range(2)])
    Dtot = 0.0
    for e in rem:
        for f in rem:
            if e == f:
                continue
            cols = [c for c in rem if c not in (e, f)]
            perm = [e, f] + cols
            pos = [rem.index(c) for c in perm]
            inv = sum(1 for a in range(len(pos)) for b in range(a + 1, len(pos)) if pos[a] > pos[b])
            Dtot += (-1) ** inv * C[0, e] * C[1, f] * G[e, f] * (np.linalg.det(M[:, cols]) if cols else 1.0)
    Dtot *= (-1) ** e0 / t[e0]
    return Dtot * U ** (-(1 + 2 * H) / 2)


prng = np.random.default_rng(5)
rows = []
for _ in range(5):
    tq = [Fr(int(prng.integers(1, 30)), int(prng.integers(1, 30))) for _ in edges]
    tf = np.array([float(x) for x in tq])
    b = brute_density(nv, edges, nv - 1, tq, lf, H, T, th.C3())
    c = closed_density_t(th, tq, lf, H, T)
    bf = float(b[0]) / math.sqrt(float(b[1]))
    mcv = mc_integrand(tf)
    vec = float(numint.density_batch(st, tf[None, :], lam, H, T)[0])
    rows.append({"brute==closed (exact)": b == c, "brute": bf, "mc_theta summand": mcv, "numint": vec,
                 "rel diff mc vs brute": abs(mcv - bf) / abs(bf)})
out["pointwise"] = rows
print("pointwise:", rows)

# ---------------------------------------------------------------- 3. deterministic integrals
det_res = []
for m in (12, 16, 24, 32):
    val, npts = numint.integrate(pqr, lam.tolist(), H, T, m)
    det_res.append({"m": m, "points": npts, "value": val, "rel_err_vs_exact": abs(val - KC3f) / abs(KC3f)})
    print("chart GL m=%d: %.15f (exact %.15f)" % (m, val, KC3f))
out["chart_GL"] = det_res
import htk  # noqa: E402
G = htk.Graph(nv, edges)
t0 = time.time()
ov, oe = htk.kernel_numeric(G, lam.tolist(), T)
s_htk = 1 if (pqr[1] + 1) * (H + T) % 2 == 0 else -1          # sgn det[B^T; C3]^D = (-1)^{(q+1)D}
out["oracle"] = {"K_htk": ov, "error_estimate": oe, "K_C3 = sgn*K_htk": s_htk * ov,
                 "rel_err_vs_exact": abs(s_htk * ov - KC3f) / abs(KC3f), "time_s": round(time.time() - t0, 1)}
print("oracle:", out["oracle"])

# ---------------------------------------------------------------- 4. the Monte Carlo estimator
def mc_values(rng, n, chunk=500000):
    vals = []
    for a in range(0, n, chunk):
        k = min(chunk, n - a)
        t = rng.dirichlet(np.ones(len(edges)), size=k)
        vals.append(numint.density_batch(st, t, lam, H, T))
    return np.concatenate(vals)


vol = 1.0 / math.factorial(len(edges) - 1)
# 4a. replay of mc_theta's stream: seed 1, the first normal(3,1) draw are the sources
rg = np.random.default_rng(1)
_ = rg.normal(size=(nv - 1, H))
v = mc_values(rg, 200000)
mean, err = v.mean() * vol, v.std() / math.sqrt(len(v)) * vol
out["replay_seed1"] = {"mean": mean, "err": err, "closed_printed": closed_printed, "exact_K_C3": KC3f,
                       "(exact-mean)/err": (KC3f - mean) / err}
print("replay of mc_theta seed 1, n=2e5: mean %.4f err %.4f ; exact %.4f ; (exact-mean)/err = %.1f" % (
    mean, err, KC3f, (KC3f - mean) / err))
# sign of the integrand
out["sign_of_integrand"] = {"min": float(v.min()), "max": float(v.max())}
# 4b. tail of the sampled values: P(|f| > x) vs  8 (c/x)^{6/5}, c = |D_theta|/2 (corner asymptotics)
rng = np.random.default_rng(12345)
vb = np.abs(mc_values(rng, 20000000))
c = abs(float(D_theta(th, lf))) / 2
tail = []
for x in (1e1, 1e2, 1e3, 1e4, 1e5):
    emp = float(np.mean(vb > x))
    tail.append({"x": x, "P_emp": emp, "8(c/x)^(6/5)": 8 * (c / x) ** 1.2})
out["tail"] = tail
print("tail:", tail)
# local slope
xs = np.array([1e2, 1e3, 1e4, 1e5])
ps = np.array([np.mean(vb > x) for x in xs])
slope = np.polyfit(np.log(xs), np.log(ps), 1)[0]
out["tail_slope_fit_1e2_1e5"] = float(slope)
print("fitted tail exponent:", slope, "(expected -6/5)")
# second moment grows with n: running sample variance
out["sample_second_moment"] = {str(n_): float(np.mean(vb[:n_] ** 2)) for n_ in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 2 * 10 ** 7)}
print("sample second moments:", out["sample_second_moment"])
# 4c. distribution of the estimator over independent runs
dist = {}
for n_, runs in ((2000, 2000), (20000, 1000), (200000, 200), (2000000, 40)):
    rng = np.random.default_rng(777 + n_)
    ratios, zs = [], []
    for _ in range(runs):
        vv = mc_values(rng, n_)
        m_, e_ = vv.mean() * vol, vv.std() / math.sqrt(n_) * vol
        ratios.append(m_ / KC3f)
        zs.append((KC3f - m_) / e_)
    ratios = np.array(ratios)
    zs = np.array(zs)
    dist[str(n_)] = {"runs": runs, "median(mean/exact)": float(np.median(ratios)),
                     "mean(mean/exact)": float(np.mean(ratios)),
                     "q10": float(np.quantile(ratios, 0.1)), "q90": float(np.quantile(ratios, 0.9)),
                     "fraction with |mean| < |exact|": float(np.mean(ratios < 1)),
                     "median of (exact-mean)/err": float(np.median(zs)),
                     "fraction |(exact-mean)/err| > 3": float(np.mean(np.abs(zs) > 3)),
                     "1 - median ratio": float(1 - np.median(ratios))}
    print(n_, dist[str(n_)], flush=True)
out["estimator_distribution"] = dist
ns = np.array([2000, 20000, 200000, 2000000], dtype=float)
defs = np.array([dist[str(int(n_))]["1 - median ratio"] for n_ in ns])
out["deficit_slope_vs_n"] = float(np.polyfit(np.log(ns), np.log(defs), 1)[0])
print("slope of log(1 - median ratio) vs log n:", out["deficit_slope_vs_n"], "(expected -1/6)")

with open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/formscheck_results.json", "w") as f:
    json.dump(out, f, indent=1)
