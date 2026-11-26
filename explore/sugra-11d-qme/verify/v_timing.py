"""Timing of the grouped evaluation (v_integrand.grouped_integrand) over all 200 m=1 flows at one random point each."""
import sys, os, time, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, PQR, H, kirchhoff, rq, rsimplex
from v_integrand import grouped_integrand, flows_m1
rng = random.Random(99)
t0 = time.time(); n = 0; nz = 0
for pqr in PQR:
    nv, edges, paths = theta_graph(*pqr)
    for fl, typ in flows_m1(nv, edges):
        lam = {v: [rq(rng) for _ in range(H)] for v in range(1, nv - 1)}; lam[0] = [Fr(0)] * H
        lam_all = dict(lam); lam_all[nv - 1] = [-sum(lam[v][i] for v in range(nv - 1)) for i in range(H)]
        ext = {v: dict(u=[rq(rng) for _ in range(H)], vg=[rq(rng) for _ in range(H)], lam=lam_all[v]) for v in range(nv)}
        x = [rsimplex(rng, k) for k in pqr]; w = rsimplex(rng, 3)
        t = [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
        y = kirchhoff(nv, edges, t, lam)
        val = grouped_integrand(nv, edges, paths, fl, typ, y, ext)
        n += 1; nz += (val != 0)
print('flows evaluated', n, '; nonzero', nz, '; total sec %.2f ; per flow %.3f s' % (time.time() - t0, (time.time() - t0) / n))
