"""C^3 K4 subdivisions: direct trapezoidal log-lattice for K(lam;0), sources e1,e2,e3,-(sum) at d1<d2<d3<d4.
usage: python3 run_c3.py <graph index 0..5> <h> <tol>"""
import sys, time
import numpy as np
import mpmath as mp
from hepp import K4Int
from lattice2 import lattice
sys.path.insert(0, '..')
from graphs import K4SUB

z3 = float(mp.zeta(3))
CLAIM = {0: (39 - 32 * z3) / 67200, 1: -1 / 77760, 2: (1296 * z3 - 1597) / 5443200,
         3: (34 - 27 * z3) / 453600, 4: 0.0, 5: (39 - 32 * z3) / 67200}

w = int(sys.argv[1]); h = float(sys.argv[2]); tol = float(sys.argv[3])
name, E = K4SUB[w]
deg = [0] * 8
for a, b in E:
    deg[a] += 1; deg[b] += 1
d2 = [v for v in range(8) if deg[v] == 2]
assert d2[3] == 7
lam = np.zeros((7, 3)); lam[d2[0]] = [1, 0, 0]; lam[d2[1]] = [0, 1, 0]; lam[d2[2]] = [0, 0, 1]
I = K4Int(8, E, [lam], m=2)
def Fc(Ex):
    R = np.ones((len(Ex), 6)); R[:, :5] = Ex
    return I.f_R(R) * np.prod(R, 1)[:, None]
t0 = time.time()
v0, b0, n0 = lattice(Fc, 5, h, tol, 0.0)
v1, b1, n1 = lattice(Fc, 5, h, tol, h / 2)
c = CLAIM[w]
print("%s h=%.2f tol=%.0e  K(shift0)=%.12e  K(shift h/2)=%.12e  claim=%.12e  rel dev %.2e %.2e  boundary-layer mass %.1e  pts %d  %.0fs"
      % (name, h, tol, v0[0], v1[0], c, (v0[0] - c) / c if c else v0[0], (v1[0] - c) / c if c else v1[0], b0, n0 + n1, time.time() - t0), flush=True)
