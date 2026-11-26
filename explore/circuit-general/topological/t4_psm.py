"""t4_psm.py -- stretch (b): Poisson sigma model on the upper half plane, heat-kernel gauge with images.

Weights I_Gamma = (2 pi)^{-2n} int_{H^n} wedge_e d phi_e  (edge order = Buring-Kiselev encoding order,
orientation dx_0 dy_0 dx_1 dy_1), computed from the circuit representation (psm.py, psm2.py):
  parallel direction: ordinary incidence, Dirichlet data x_L = 0, x_R = 1;
  perpendicular direction: signed incidence (image edge: +1, +1), ground y = 0, plus one T = 1 edge
  v -> ground per aerial vertex (sgn y_v / 2), integrated analytically;
  edge times t = lambda u: radial integral analytic (P1 = 0), Hepp sectors + tanh-sinh on the rest.
Benchmark: Buring-Kiselev arXiv:2209.14438v2 (HTML), App. A.2, star product
  h^1: "2 1 1  0 1  1";  h^2: "1 3 0 2  -1/6", "1 3 0 1  -1/3", "0 3 0 1  1/3", "0 1 0 1  1/2";
  coefficient c_Gamma = (# labelled graphs / n!) I_Gamma  (2/1 for n = 1; 8/2 for n = 2 without
  automorphism; 4/2 for the product graph).
"""
import math, time
import numpy as np
import psm, psm2

GRAPHS = [
    ("h^1 wedge        '0 1'", 1, [(0, 'L'), (0, 'R')], 2, 1.0),
    ("h^2 product  '0 1 0 1'", 2, [(0, 'L'), (0, 'R'), (1, 'L'), (1, 'R')], 4, 1 / 2),
    ("h^2 Bernoulli '0 3 0 1'", 2, [(0, 'L'), (0, 1), (1, 'L'), (1, 'R')], 8, 1 / 3),
    ("h^2 Bernoulli '1 3 0 1'", 2, [(0, 'R'), (0, 1), (1, 'L'), (1, 'R')], 8, -1 / 3),
    ("h^2 two-cycle '1 3 0 2'", 2, [(0, 'R'), (0, 1), (1, 'L'), (1, 0)], 8, -1 / 6),
]


def main():
    t0 = time.time()
    print(__doc__)
    rng = np.random.default_rng(11)
    # (1) P1 = 0 pointwise
    print("(1) P1 / P0 at 20 random points (parallel covariance terms cancel):")
    for name, n, ed, lab, bk in GRAPHS:
        if n == 2:
            u = rng.dirichlet([1] * 4, size=20)
            P0, P1, E = psm2.Jparts(2, ed, u, noexp=True)
            print("   %-26s max |P1/P0| = %.1e" % (name, np.max(np.abs(P1 / P0))))
    # (2) tau-integration check: J_tau(t) (psm2) = int dtau J(t, tau) (psm, aux edges explicit)
    print("\n(2) auxiliary edges: int dtau J(t,tau) (2-d quadrature in log tau) vs closed form with sgn moments:")
    x, w = psm2.tanh_sinh(0.08, 60)
    for name, n, ed, lab, bk in GRAPHS[4:5] + GRAPHS[2:3]:
        for rep in range(2):
            t = rng.random((1, 4)) * 2 + 0.2
            J0, J1, E = psm2.Jparts(2, ed, t)
            # tau_v = s/(1-s), s in (0,1)
            S1, S2 = np.meshgrid(x, x, indexing='ij')
            W1, W2 = np.meshgrid(w, w, indexing='ij')
            ta = np.stack([S1.ravel() / (1 - S1.ravel()), S2.ravel() / (1 - S2.ravel())], axis=1)
            jac = 1 / (1 - S1.ravel()) ** 2 / (1 - S2.ravel()) ** 2
            tt = np.repeat(t, ta.shape[0], axis=0)
            J = psm.J_batch(2, ed, tt, ta)
            val = np.sum(J * jac * W1.ravel() * W2.ravel())
            print("   %-26s t=%s: quadrature %.10f, closed form %.10f" % (name, np.round(t[0], 3), val, (J0 + J1)[0]))
    # (3) weights
    print("\n(3) weights (Hepp sectors, tanh-sinh step h):")
    res = {}
    for name, n, ed, lab, bk in GRAPHS:
        vals = []
        for h in (0.3, 0.2):
            vals.append(psm2.weight_sectors(n, ed, h=h))
        I = vals[-1]
        c = lab / math.factorial(n) * I
        res[name] = (I, c)
        print("   %-26s I(h=0.3) = %.9f, I(h=0.2) = %.9f ;  c = %d/%d! * I = %.7f ; Buring-Kiselev %s ; rel. dev %.1e"
              % (name, vals[0], vals[1], lab, n, c, ("%.7f" % bk), abs(c - bk) / abs(bk)))
    print("\nexact values (wedge analytic; others identified): I = 1/2, 1/4, 1/12, -1/12, -1/24")
    print("time %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
