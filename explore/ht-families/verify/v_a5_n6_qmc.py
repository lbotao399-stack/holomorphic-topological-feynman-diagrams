"""a5 at n = 6, (H,T) = (1,1), own QMC (no engine input): banana (2,2,2,2,2,2,1) of the b21.py construction,
and the explorer's n=6 graph 1 with its printed sources [3,1,-3,1,2,-3,1]; compared with
F = (5!!/2^3) D_Gamma prod_P Gamma(3/2+1-a_P)/(a_P-1)! / Gamma(1+7/2) = D_Gamma pi^3/7."""
import sys, math, time, random
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import vcore as vc
import vsp
from v_b2 import banana, D_banana


def qmc(nv, edges, lamf, m, nrep, seed):
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    rc = vsp.RatioCoords(vsp.sp_decompose(nv, edges), net.N)
    f = lambda t, J, full: vc.density_batch(nv, edges, C, t, lamf) * J
    ex = vsp.probe_exponents(rc, f, random.Random(seed), npts=2)
    alpha = {b: max(0.5, ex[b] + 1) for b in rc.names}
    return vsp.qmc_integrate(rc, f, alpha, m=m, nrep=nrep, seed=seed, batch=2048)


nv, edges, ints = banana((2, 2, 2, 2, 2, 2, 1))
for lam_i in ([2, -1, 3, 1, -2, 1, 2], [1, 3, -2, 2, 1, -1, 3]):
    lam = [[Fr(x)] for x in lam_i]
    t0 = time.time()
    mean, err, vals, nex = qmc(nv, edges, np.array([[float(x)] for x in lam_i]), 10, 8, 3)
    Dg = D_banana(lam, ints, 1)
    F = float(Dg) * math.pi ** 3 / 7
    print("banana(2,2,2,2,2,2,1) lambda %s: own QMC K/pi^3 = %.8f +- %.1e ; F/pi^3 = D_Gamma/7 = %.8f  (%.0fs)" % (
        lam_i, mean / math.pi ** 3, err / math.pi ** 3, F / math.pi ** 3, time.time() - t0), flush=True)
edges1 = [(0, 5), (0, 6), (1, 5), (1, 6), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
lam_i = [3, 1, -3, 1, 2, -3, 1]
t0 = time.time()
mean, err, vals, nex = qmc(8, edges1, np.array([[float(x)] for x in lam_i]), 10, 8, 5)
inter = lam_i[:5] + [-sum(lam_i)]
print("explorer n=6 graph 1, lambda %s: own QMC K/pi^3 = %.8f +- %.1e ; |prod interior sources|/7 = %.8f  (%.0fs)" % (
    lam_i, mean / math.pi ** 3, err / math.pi ** 3, abs(math.prod(inter)) / 7, time.time() - t0), flush=True)
