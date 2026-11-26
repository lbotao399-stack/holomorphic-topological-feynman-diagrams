import time, numpy as np
from hepp import K4Int
from lattice2 import lattice
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
lam = [[0, 1], [0, 0], [-1, -1], [1, 0]]
I = K4Int(5, C2N3, [lam], m=2)
def Fc(E):
    R = np.ones((len(E), 6)); R[:, :5] = E
    return I.f_R(R) * np.prod(R, 1)[:, None]
ex = -1/576
for h, tol in [(0.8, 1e-9), (0.7, 1e-10)]:
    t0 = time.time()
    v0, b0, n0 = lattice(Fc, 5, h, tol, 0.0)
    v1, b1, n1 = lattice(Fc, 5, h, tol, h/2)
    print("h=%.2f tol=%.0e  K0=%.13e K1=%.13e rel0=%.2e rel1=%.2e below=%.1e pts=%d %.0fs" % (h, tol, v0[0], v1[0], abs(v0[0]/ex-1), abs(v1[0]/ex-1), b0, n0+n1, time.time()-t0), flush=True)
