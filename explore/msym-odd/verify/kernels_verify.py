"""Independent check of the claim 'theta-kernels' (and the Q1-zero density) with vkern.py:
fresh Kirchhoff solver, BFS-tree cycle basis made positive, rho_2 and psi as exact forms on dt-space,
coefficient extraction checked for every omitted edge (basicness).
(1) pointwise: D(t(w,x)) |det[t; dt/dxi]| = s (D_Gamma/2) prod_P w_P^{H+1/2-n_P}  (exact, square and sign);
(2) D_Gamma against the explicit determinants of the claim;
(3) K(lam;0) / (pi D_Gamma) from the density by Dirichlet integrals;
(4) displacements -sum_e y_e.z_e = sum_v lam_v.X_v with the fresh currents (general rule and the
    explicit X maps of the claim for (A), (B)).
"""
import random
import time
from fractions import Fraction as Fr

from vkern import Net, Theta, theta_rootN, S2_A, S2_B, det, dirichlet

T0 = time.time()


def d2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def d3(a, b, c):
    return det([[a[i], b[i], c[i]] for i in range(3)])


cases = []
cases.append(("5d (A) S2 data", S2_A, 2, Fr(1, 20), +1,
              lambda L: -d2(L[0], L[1]) * d2(L[3], L[4])))
cases.append(("5d (B) S2 data", S2_B, 2, Fr(-1, 30), -1,
              lambda L: -d2(L[0], L[1]) * d2(L[3], L[4])))
for pqr in [(1, 3, 3), (2, 2, 3)]:
    nv, ed, pa = theta_rootN(*pqr)
    cases.append(("5d theta_rootN%s" % (pqr,), (nv, ed, pa), 2, None, None, None))
nv, ed, pa = theta_rootN(1, 4, 4)     # interior: B = 0,1,2 ; C = 3,4,5
cases.append(("7d theta(1,4,4)", (nv, ed, pa), 3, Fr(1, 252), +1, lambda L: d3(L[0], L[1], L[2]) * d3(L[3], L[4], L[5])))
nv, ed, pa = theta_rootN(2, 3, 4)     # A = 0; B = 1,2; C = 3,4,5
cases.append(("7d theta(2,3,4)", (nv, ed, pa), 3, Fr(1, 420), +1, lambda L: d3(L[0], L[1], L[2]) * d3(L[3], L[4], L[5])))
nv, ed, pa = theta_rootN(3, 3, 3)     # A = 0,1; B = 2,3; C = 4,5
cases.append(("7d theta(3,3,3)", (nv, ed, pa), 3, Fr(1, 840), +1,
              lambda L: d3(L[0], L[1], L[2]) * d3(L[3], L[4], L[5]) - d3(L[0], L[1], L[3]) * d3(L[2], L[4], L[5])))
nv, ed, pa = theta_rootN(1, 2, 2)     # B = 0; C = 1
cases.append(("3d theta(1,2,2)", (nv, ed, pa), 1, Fr(-1, 3), -1, lambda L: L[0][0] * L[1][0]))

rng = random.Random(20261120)
for name, (nv, ed, pa), H, cclaim, sclaim, Dexpl in cases:
    t0 = time.time()
    TH = Theta(nv, ed, pa, H)
    signs = TH.check(npts=8, seed=rng.randint(0, 10 ** 6))
    assert len(signs) == 1
    s = signs.pop()
    # (2) explicit determinant
    okD = None
    if Dexpl is not None:
        okD = True
        for _ in range(5):
            lam = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)]
            okD &= TH.Dtheta(lam) == Dexpl(lam)
    # (3) K(lam;0) / (pi D) = s/2 * Dirichlet_w(alpha) * prod_P 1/(n_P-1)!
    q, k = dirichlet(TH.alpha)
    assert k == 2, k                     # sqrt(pi)^2 = pi
    c = Fr(s, 2) * q
    for n in TH.nP:
        qq, kk = dirichlet([0] * n)
        assert kk == 0
        c *= qq
    # (4) displacement rule with z on the spokes (first edge of each path, oriented from N), random lam incl. N'
    okX = True
    for _ in range(4):
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)]
        w, x = TH.point(rng)
        t = TH.t_of(w, x)
        y = TH.net.currents(t, lam)
        zP = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(3)]
        lhs = Fr(0)
        for P in range(3):
            e = pa[P][0]
            ze = [TH.sg[P][0] * c_ for c_ in zP[P]]       # edge displacement in the edge's own orientation
            lhs -= sum(y[e][i] * ze[i] for i in range(H))
        X = TH.X_rule(w, x, zP)
        rhs = sum(lam[v][i] * X[v][i] for v in X for i in range(H))
        okX &= lhs == rhs
    print("%-24s |I|=%d |V|=%d n_P=%s alpha=%s: pointwise density (8 points) single sign s=%+d%s;"
          " D_Gamma explicit: %s; K(lam;0)/(pi D) = %s%s; X rule (4 points): %s  (%.1fs)"
          % (name, TH.net.N, nv, TH.nP, [str(a) for a in TH.alpha], s,
             "" if sclaim is None else (" (claim %+d: %s)" % (sclaim, s == sclaim)), okD, c,
             "" if cclaim is None else (" (claim %s: %s)" % (cclaim, c == cclaim)), okX, time.time() - t0))

# explicit X maps of the claim for (A) and (B), z_1 on edge 0, z_2 on edge 3 (A) / 2 (B), -z_3 on edge 6
print("explicit X maps of the claim (5d):")
for name, data, zedges in [("A", S2_A, (0, 3, 6)), ("B", S2_B, (0, 2, 6))]:
    nv, ed, pa = data
    TH = Theta(nv, ed, pa, 2)
    ok = True
    for _ in range(6):
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(nv - 1)]
        lam[2] = [Fr(0), Fr(0)]                          # N' carries no source (no output)
        w, x = TH.point(rng)
        wA, wB, wC = w
        xA, xB, xC = x                                   # x^P_j = x[P][j-1]
        t = TH.t_of(w, x)
        y = TH.net.currents(t, lam)
        z = [[Fr(rng.randint(-9, 9)) for _ in range(2)] for _ in range(3)]
        ze = {zedges[0]: z[0], zedges[1]: z[1], zedges[2]: [-c_ for c_ in z[2]]}
        lhs = -sum(y[e][i] * ze[e][i] for e in ze for i in range(2))
        if name == "A":
            Zb = [wB * z[0][i] + wA * z[1][i] + wC * z[2][i] for i in range(2)]
            s12 = xB[1] + xB[2]; s3 = xB[2]; c3 = xC[2]; c23 = xC[1] + xC[2]
            X = {0: [s12 * z[0][i] + (1 - s12) * Zb[i] for i in range(2)],
                 1: [s3 * z[0][i] + (1 - s3) * Zb[i] for i in range(2)],
                 3: [c3 * z[2][i] + (1 - c3) * Zb[i] for i in range(2)],
                 4: [c23 * z[2][i] + (1 - c23) * Zb[i] for i in range(2)]}
        else:
            Zb = [wA * z[0][i] + wB * z[1][i] + wC * z[2][i] for i in range(2)]
            a2 = xA[1]; b2 = xB[1]; c3 = xC[2]; c23 = xC[1] + xC[2]
            X = {0: [a2 * z[0][i] + (1 - a2) * Zb[i] for i in range(2)],
                 1: [b2 * z[1][i] + (1 - b2) * Zb[i] for i in range(2)],
                 3: [c3 * z[2][i] + (1 - c3) * Zb[i] for i in range(2)],
                 4: [c23 * z[2][i] + (1 - c23) * Zb[i] for i in range(2)]}
        rhs = sum(lam[v][i] * X[v][i] for v in X for i in range(2))
        ok &= lhs == rhs
    print("   (%s): -sum y.z = sum lam.X at 6 rational points: %s" % (name, ok))

# Q1-zero: n = 1 polygon, T = 1: rho_1 = (det M)^{-1/2} Pf(1x1) = 0
print("Q1-zero: n = 1 Pfaffian of a 1x1 antisymmetric matrix of 2-forms is 0 by definition; the 2-form")
for H in (1, 2, 3):
    m = H + 2
    net = Net(m, [(i, (i + 1) % m) for i in range(m)])
    C = net.C
    tt = [Fr(random.Random(H).randint(1, 20), 7) for _ in range(m)]
    Mv = sum(C[0][e] ** 2 * tt[e] for e in range(m))
    # -1/2 dM M^{-1} dM for n = 1: dM ^ dM / M = 0 identically
    dM = {(e,): Fr(C[0][e] ** 2) for e in range(m)}
    from vkern import fwedge
    print("   H=%d polygon %d: dM ^ dM = %s (so -1/2 dM M^-1 dM = 0)" % (H, m, fwedge(dM, dM)))
print("total %.1fs" % (time.time() - T0))
