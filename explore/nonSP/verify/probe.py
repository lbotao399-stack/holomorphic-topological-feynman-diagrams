import numpy as np
from hepp import K4Int
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
lam = [[0, 1], [0, 0], [-1, -1], [1, 0]]
I = K4Int(5, C2N3, [lam], m=3)
sg = (1, 0, 2, 5, 3, 4)
def F(u):
    u = np.atleast_2d(u)
    Rk = np.ones((len(u), 6))
    for k in range(1, 6):
        Rk[:, k] = Rk[:, k-1] * u[:, k-1]
    R = np.empty_like(Rk); R[:, list(sg)] = Rk
    J = u[:, 0]**4 * u[:, 1]**3 * u[:, 2]**2 * u[:, 3]
    return I.f_R(R)[:, 0] * J
base = np.array([0.5, 0.4, 0.6, 0.3, 0.7])
for k in range(5):
    vals = []
    for s in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        u = base.copy(); u[k] = s
        vals.append(F(u)[0])
    print(k, ["%.6e" % v for v in vals])
