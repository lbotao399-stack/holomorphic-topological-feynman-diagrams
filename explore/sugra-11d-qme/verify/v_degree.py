"""Claim a3, derivative count at m = 1: flow integrand (grouped evaluation) is homogeneous of degree 13 = (12-m)+2m
in the momenta (currents and external lam scale together); D_theta adds 10: total 23 = 22+m."""
import sys, os, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, PQR, H, kirchhoff, rq, rsimplex, D_theta
from v_integrand import grouped_integrand, flows_m1
rng = random.Random(5)
for pqr in PQR:
    nv, edges, paths = theta_graph(*pqr)
    fls = flows_m1(nv, edges)
    ok = True
    for fl, typ in (fls[0], fls[-1]):
        lam = {v: [rq(rng) for _ in range(H)] for v in range(1, nv - 1)}; lam[0] = [Fr(0)] * H
        u = {v: [rq(rng) for _ in range(H)] for v in range(nv)}; vg = {v: [rq(rng) for _ in range(H)] for v in range(nv)}
        x = [rsimplex(rng, k) for k in pqr]; w = rsimplex(rng, 3)
        t = [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
        vals = []
        for s in (Fr(1), Fr(3)):
            ls = {v: [s * c for c in lam[v]] for v in lam}
            la = dict(ls); la[nv - 1] = [-sum(ls[v][i] for v in range(nv - 1)) for i in range(H)]
            ext = {v: dict(u=u[v], vg=vg[v], lam=la[v]) for v in range(nv)}
            y = kirchhoff(nv, edges, t, ls)
            vals.append(grouped_integrand(nv, edges, paths, fl, typ, y, ext) * D_theta(ls, paths))
        ok &= (vals[0] != 0 and vals[1] == Fr(3) ** 23 * vals[0])
    print(pqr, 'integrand x D_theta homogeneous of degree 23 in lambda (2 flows):', ok)

# torus weight (b3 premise): lam_i -> s_i lam_i, y_i -> s_i y_i, gamma_i -> s_i gamma_i, mu^i -> mu^i / s_i
# => integrand x D_theta -> (s_1...s_5)^3 x (integrand x D_theta)
for pqr in PQR:
    nv, edges, paths = theta_graph(*pqr)
    fls = flows_m1(nv, edges)
    ok = True
    for fl, typ in (fls[1], fls[-2]):
        lam = {v: [rq(rng) for _ in range(H)] for v in range(1, nv - 1)}; lam[0] = [Fr(0)] * H
        u = {v: [rq(rng) for _ in range(H)] for v in range(nv)}; vg = {v: [rq(rng) for _ in range(H)] for v in range(nv)}
        x = [rsimplex(rng, k) for k in pqr]; w = rsimplex(rng, 3)
        t = [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
        sc = [Fr(2), Fr(3), Fr(5), Fr(7), Fr(11, 13)]
        vals = []
        for S in ([Fr(1)] * H, sc):
            ls = {v: [S[i] * lam[v][i] for i in range(H)] for v in lam}
            la = dict(ls); la[nv - 1] = [-sum(ls[v][i] for v in range(nv - 1)) for i in range(H)]
            ext = {v: dict(u=[u[v][i] / S[i] for i in range(H)], vg=[S[i] * vg[v][i] for i in range(H)], lam=la[v]) for v in range(nv)}
            y = kirchhoff(nv, edges, t, ls)
            vals.append(grouped_integrand(nv, edges, paths, fl, typ, y, ext) * D_theta(ls, paths))
        P5 = Fr(1)
        for s in sc: P5 *= s
        ok &= (vals[0] != 0 and vals[1] == P5 ** 3 * vals[0])
    print(pqr, 'torus weight (3,3,3,3,3) of integrand x D_theta (2 flows):', ok)
