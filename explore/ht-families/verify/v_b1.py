"""Verification of b1 (and of the n=4 part of a5 for G0): (H,T) = (1,1), n = 4, graphs G0..G4 as stated.
 V1  exact: htk RatioDensity at integer sources (plan 'large'), checked pointwise at 4 rational points against
     the own brute force (vcore: DFS tree, own Kirchhoff solve, own Pfaffian of forms) with the own exact
     Jacobian det[t; dt/dxi] of the engine's t-map; then integrated with the own Dirichlet code (vdir).
     Compared with the claimed formula at the same sources.  3 random source sets per graph.
 V2  numeric, no engine input: own SP decomposition and ratio map (vsp), own float density, scrambled Sobol
     QMC with Dirichlet(alpha) sampling, alpha_b = 1/2 if the U-exponent of b is odd, else 1.
"""
import sys, random, time, math
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import vcore as vc
import vsp
import vdir
import htk
from htk import spred

GR = {
    0: [(0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)],
    1: [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (2, 5), (3, 4)],
    2: [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)],
    3: [(0, 1), (0, 2), (0, 5), (1, 2), (1, 3), (1, 4), (1, 5), (3, 4), (4, 5)],
    4: [(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (2, 3), (3, 4), (4, 5)],
}
l0, l1, l2, l3, l4 = L = sp.symbols("l0:5")
S = l0 + l1 + l2 + l3 + l4
FORM = {
    0: sp.pi ** 2 * l0 * l1 * l2 * l3 / 5,
    1: -sp.pi ** 2 * l1 * l4 * (3 * l0 - 3 * l2 + 7 * l3 + 5 * l4) * S / 90,
    2: sp.Integer(0),
    3: -sp.pi ** 2 * l2 * l3 * (8 * l0 ** 2 + 20 * l0 * l1 + 18 * l0 * l2 + 20 * l0 * l3 + 20 * l0 * l4 + 12 * l1 ** 2 + 22 * l1 * l2
                                + 22 * l1 * l3 + 20 * l1 * l4 + 10 * l2 ** 2 + 21 * l2 * l3 + 20 * l2 * l4 + 10 * l3 ** 2 + 18 * l3 * l4 + 8 * l4 ** 2) / 180,
    4: -sp.pi ** 2 * l1 * S * (l0 * l1 + 2 * l0 * l3 - 9 * l1 ** 2 - 17 * l1 * l2 - 19 * l1 * l3 - l1 * l4 - 8 * l2 ** 2 - 18 * l2 * l3
                               - 10 * l3 ** 2 - 2 * l3 * l4) / 180,
}


def engine_exact(nv, edges, lam, rng, npts=4):
    G = htk.Graph(nv, edges)
    red = spred.sp_reduce(G, prefer="large")
    dens = red.density(lam, T=1, check=0)
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    ok = True
    for _ in range(npts):
        pt = red.complete_point({c: v for c, v in red.random_point(rng).items() if c in red.coords})
        t = red.tmap(pt)
        rows = []
        for c in red.coords:
            blk = next(b for b in red.blocks if c in b.bary)
            b0 = blk.bary[0]
            rows.append([t[e] * (Fr(red.texp[e].get(c, 0)) / pt[c] - Fr(red.texp[e].get(b0, 0)) / pt[b0]) for e in range(net.N)])
        J = vc.det_fr([t] + rows)
        coef, U, _ = vc.density_slice(net, C, t, lam)
        mine = coef * abs(J)                         # g = mine * U^{-1/2}
        rat, rad = dens.evaluate(tuple(pt[c] for c in red.coords))
        ok &= ((mine > 0) == (rat > 0)) and ((mine < 0) == (rat < 0)) and mine ** 2 / U == rat ** 2 * rad
    K = vdir.integrate_engine_density(red, dens)
    return K, ok, dens


def own_qmc(nv, edges, lamf, m, nrep, seed, shift=1.0, thresh=None):
    """own QMC; alpha_b = max(1/2, e_b + shift), e_b = probed leading exponent of b in g;
    exact re-evaluation (own exact density) where max t / min t > thresh (1e7 for |I| = 9, 1e6 otherwise)."""
    if thresh is None:
        thresh = 1e7 if len(edges) <= 9 else 1e6
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    rc = vsp.RatioCoords(vsp.sp_decompose(nv, edges), net.N)
    f = lambda t, J, full: vc.density_batch(nv, edges, C, t, lamf) * J
    lamF = [[Fr(float(x)) for x in row] for row in lamf]

    def exact_cb(xi):
        full = rc.complete(xi)
        tt = rc.t_of(full)
        Jx = vc.det_fr([tt] + rc.jac_rows(full, tt))
        coef, U, _ = vc.density_slice(net, C, tt, lamF)
        return float(coef * abs(Jx)) / math.sqrt(float(U))

    ex = vsp.probe_exponents(rc, f, random.Random(seed))
    alpha = {b: max(0.5, ex[b] + shift) for b in rc.names}
    mean, err, vals, nex = vsp.qmc_integrate(rc, f, alpha, m=m, nrep=nrep, seed=seed, exact_cb=exact_cb, thresh=thresh)
    return mean, err, nex


if __name__ == "__main__":
    which = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else range(5)
    m = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    rng = random.Random(2026)
    for k in which:
        edges = GR[k]
        print("=" * 90)
        print("G%d edges %s root 5" % (k, edges))
        for trial in range(3):
            lam = [[Fr(rng.choice([-7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7]))] for _ in range(5)]
            t0 = time.time()
            K, ok, dens = engine_exact(6, edges, lam, rng)
            claim = FORM[k].subs({s: int(l[0]) for s, l in zip(L, lam)})
            Kexpr = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())
            print("  lambda %s" % [int(l[0]) for l in lam])
            print("    V1 pointwise engine g == own brute force (4 pts): %s ; own Dirichlet K = %s ; claim = %s ; equal: %s  (%.1fs)" % (
                ok, Kexpr, claim, sp.simplify(Kexpr - claim) == 0, time.time() - t0))
            if trial < 2:
                t0 = time.time()
                lamf = np.array([[float(l[0])] for l in lam])
                mean, err, nex = own_qmc(6, edges, lamf, m, 16, 100 * k + trial)
                ex = float(claim)
                print("    V2 own QMC 2^%d x 16: %.8g +- %.2g ; claim %.10g ; diff/err = %.2f ; rel diff %.1e ; %d points exact  (%.1fs)" % (
                    m, mean, err, ex, (mean - ex) / err, abs(mean - ex) / abs(ex) if ex else float("nan"), nex, time.time() - t0))
