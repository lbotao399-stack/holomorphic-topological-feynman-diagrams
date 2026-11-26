"""Self-tests of the own brute-force layer vcore/vsp.
 S1  slice formula with e0 = N-1 and e0 = 0 agree (exact)
 S2  pullback of Omega to own ratio coordinates == D(t(xi)) det[t; dt/dxi] (exact)
 S3  float batch density == exact density (rel. error)
 S4  projective volume test: int |det[t;dt/dxi]| (sum t)^{-N} d xi = 1/(N-1)! (QMC, own coordinates)
 S5  T=0 control: 6d box of S2 (H=3), K = (1/6) det(l1,l2,l3) up to the sign (-1)^{n(n+1)/2} of CONVENTIONS 3;
     own QMC
 S6  cross-check of the density against htk.brute.density at 2 rational points (engine as reference only)
"""
import sys, random, math, time
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import vcore as vc
import vsp

G2 = (6, [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)])
G4 = (6, [(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (2, 3), (3, 4), (4, 5)])
rng = random.Random(11)

for name, (nv, edges) in [("G2", G2), ("G4", G4)]:
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    lam = [[Fr(rng.randint(-6, 6))] for _ in range(nv - 1)]
    t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(net.N)]
    c1, U1, _ = vc.density_slice(net, C, t, lam, e0=net.N - 1)
    c0, U0, _ = vc.density_slice(net, C, t, lam, e0=0)
    c3, U3, _ = vc.density_slice(net, C, t, lam, e0=3)
    print("S1 %s: D(e0=N-1) == D(e0=0) == D(e0=3): %s   (coef %s, U %s)" % (name, c1 == c0 == c3, c1, U1))
    # S2
    root = vsp.sp_decompose(nv, edges)
    rc = vsp.RatioCoords(root, net.N)
    ok = True
    for _ in range(3):
        xi = rc.random_point(rng)
        full = rc.complete(xi)
        tt = rc.t_of(full)
        rows = rc.jac_rows(full, tt)
        J = vc.det_fr([tt] + rows)
        Jt = [[rows[c][e] for c in range(net.N - 1)] for e in range(net.N)]
        dts = [vc.oneform(Jt[e]) for e in range(net.N)]
        f, U, y = vc.omega(net, C, tt, lam, dts, T=1)
        gsig = f.get((1 << (net.N - 1)) - 1, Fr(0))
        cD, UD, _ = vc.density_slice(net, C, tt, lam)
        ok &= (gsig == cD * J) and U == UD
    print("S2 %s: pullback to own ratio coordinates == D det[t; dt/dxi] at 3 points: %s" % (name, ok))
    # S3
    tf = np.array([[float(x) for x in t]])
    Df = vc.density_batch(nv, edges, C, tf, np.array([[float(l[0])] for l in lam]))[0]
    De = float(c1) / math.sqrt(float(U1))
    print("S3 %s: float %.15g exact %.15g rel %.1e" % (name, Df, De, abs(Df - De) / abs(De)))
    # S4
    alpha = {b: 1.0 for b in rc.names}
    t0 = time.time()
    mean, err, _, _ = vsp.qmc_integrate(rc, lambda t, J, full: J / t.sum(axis=1) ** net.N, alpha, m=13, nrep=8, seed=5)
    print("S4 %s: int |J| (sum t)^-N = %.6g +- %.1g ; 1/(N-1)! = %.6g  (%.1fs)" % (
        name, mean, err, 1 / math.factorial(net.N - 1), time.time() - t0))

# S6: engine reference
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import brute
for name, (nv, edges) in [("G2", G2), ("G4", G4)]:
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    G = htk.Graph(nv, edges)
    for _ in range(2):
        lam = [[Fr(rng.randint(-6, 6))] for _ in range(nv - 1)]
        t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(net.N)]
        c1, U1, _ = vc.density_slice(net, C, t, lam)
        D = brute.density(G, t, lam, T=1)
        same = (c1 > 0) == (D.coeff > 0) and c1 ** 2 / U1 == D.coeff ** 2 / D.radicand
        print("S6 %s: own density == htk.brute.density (sign and square): %s" % (name, same))

# S7: bijectivity of the own ratio map (exact): ratios_of_t(t(xi)) == xi, and t(ratios_of_t(t)) == t / R_root
for name, (nv, edges) in [("G2", G2), ("G4", G4)]:
    net = vc.Net(nv, edges)
    ok = True
    for trial in range(5):
        root = vsp.sp_decompose(nv, edges, rng=random.Random(trial))
        rc = vsp.RatioCoords(root, net.N)
        xi = rc.random_point(rng)
        full = rc.complete(xi)
        tt = rc.t_of(full)
        back, Rr = rc.ratios_of_t(tt)
        ok &= all(back[b] == full[b] for b in rc.names) and Rr == 1
        t = [Fr(rng.randint(1, 50), rng.randint(1, 50)) for _ in range(net.N)]
        r, Rr = rc.ratios_of_t(t)
        t2 = rc.t_of(r)
        ok &= all(a * Rr == b for a, b in zip(t2, t)) and all(sum(r[b] for b in nm) == 1 for _, nm in rc.blocks)
    print("S7 %s: own ratio map bijective (5 random decompositions, exact): %s" % (name, ok))

# S4': projective volume with Dirichlet(1/2) sampling, more points
for name, (nv, edges) in [("G2", G2), ("G4", G4)]:
    net = vc.Net(nv, edges)
    rc = vsp.RatioCoords(vsp.sp_decompose(nv, edges), net.N)
    alpha = {b: 0.5 for b in rc.names}
    mean, err, _, _ = vsp.qmc_integrate(rc, lambda t, J, full: J / t.sum(axis=1) ** net.N, alpha, m=16, nrep=16, seed=9)
    print("S4' %s: int |J| (sum t)^-N = %.7g +- %.1g ; 1/(N-1)! = %.7g ; rel %.1e" % (
        name, mean, err, 1 / math.factorial(net.N - 1), abs(mean * math.factorial(net.N - 1) - 1)))
