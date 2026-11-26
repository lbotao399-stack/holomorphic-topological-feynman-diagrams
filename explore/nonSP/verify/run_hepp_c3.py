"""Second quadrature for C^3 graphs: all 720 Hepp sectors in the 6 chain resistances, tensor Gauss-Legendre n^5.
usage: python3 run_hepp_c3.py <graph index> <n>"""
import sys, time
import numpy as np
from hepp import K4Int
sys.path.insert(0, '..')
from graphs import K4SUB
w, n = int(sys.argv[1]), int(sys.argv[2])
name, E = K4SUB[w]
deg = [0] * 8
for a, b in E:
    deg[a] += 1; deg[b] += 1
d2 = [v for v in range(8) if deg[v] == 2]
lam = np.zeros((7, 3)); lam[d2[0]] = [1, 0, 0]; lam[d2[1]] = [0, 1, 0]; lam[d2[2]] = [0, 0, 1]
I = K4Int(8, E, [lam], m=2)
t0 = time.time()
v, per = I.sector_integral(n, chunk=8000)
per = np.array([p[0] for p in per])
print("%s Hepp n=%d: K = %.10e   sum|sector| = %.4e  sectors >0: %d <0: %d  %.0fs" % (name, n, v[0], np.abs(per).sum(), (per > 0).sum(), (per < 0).sum(), time.time() - t0), flush=True)
