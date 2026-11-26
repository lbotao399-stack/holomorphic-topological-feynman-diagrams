"""b5: C^2 n=3 non-SP graph [(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,4)], direct log-lattice (vk density, no reduction).
Sources lam3 = e1, lam0 = e2, lam2 = (a, b), lam1 = 0, lam4 = -(sum): p = <l3,l0> = 1, q = <l3,l2> = b, s = <l0,l2> = -a.
16 source configurations (a, b) in {-1,0,1,2}^2 integrated simultaneously; cubic form in (p,q,s) fitted by least squares.
usage: python3 run_c2n3.py <h> <tol>"""
import sys, time, itertools
import numpy as np
from hepp import K4Int
from lattice2 import lattice

h = float(sys.argv[1]); tol = float(sys.argv[2])
E = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
grid = [(a, b) for a in (-1, 0, 1, 2) for b in (-1, 0, 1, 2)]
lams = [[[0, 1], [0, 0], [a, b], [1, 0]] for (a, b) in grid]
I = K4Int(5, E, lams, m=2)
def Fc(Ex):
    R = np.ones((len(Ex), 6)); R[:, :5] = Ex
    return I.f_R(R) * np.prod(R, 1)[:, None]
t0 = time.time()
v0, b0, n0 = lattice(Fc, 5, h, tol, 0.0)
v1, b1, n1 = lattice(Fc, 5, h, tol, h / 2)
mons = [(i, j, k) for i in range(4) for j in range(4) for k in range(4) if i + j + k == 3]
def basis(a, b):
    p, q, s = 1.0, b, -a
    return [p ** i * q ** j * s ** k for (i, j, k) in mons]
A = np.array([basis(a, b) for (a, b) in grid])
claim = {(3, 0, 0): -1, (2, 1, 0): -1, (2, 0, 1): 3, (1, 2, 0): 5, (1, 1, 1): 12, (1, 0, 2): 3}
print("h=%.2f tol=%.0e points %d  %.0fs  boundary-layer mass %.1e" % (h, tol, n0 + n1, time.time() - t0, b0))
for lab, v in (("shift 0", v0), ("shift h/2", v1)):
    co, res, rk, sv = np.linalg.lstsq(A, v, rcond=None)
    resid = np.max(np.abs(A @ co - v))
    print(lab, "max fit residual %.1e" % resid)
    for m, c in zip(mons, co):
        print("   p^%d q^%d s^%d : 576*coef = %+.9f   claim %+d   diff %.1e" % (m + (576 * c, claim.get(m, 0), 576 * c - claim.get(m, 0))))
# the two sample points of the explorer
for (a, b), lab in (((1 / 3, 2 / 5), "lam2=(1/3,2/5)"),):
    pass
print("K at (a,b) = (0,0) [sources at 3,0 and 4 only]: shift0 %.12e  shift h/2 %.12e  (-1/576 = %.12e)" % (v0[grid.index((0, 0))], v1[grid.index((0, 0))], -1 / 576))
