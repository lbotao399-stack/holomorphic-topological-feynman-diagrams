import time, numpy as np, itertools
from hepp import K4Int
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
lam = [[0, 1], [0, 0], [-1, -1], [1, 0]]
I = K4Int(5, C2N3, [lam], m=3)
res = {}
for n in (4, 6, 8):
    v, per = I.sector_integral(n)
    res[n] = np.array([p[0] for p in per])
    print(n, v[0], np.abs(res[n]).max(), np.abs(res[n]).sum(), flush=True)
d = np.abs(res[8]-res[6])
idx = np.argsort(-d)[:10]
perms = list(itertools.permutations(range(6)))
print("chains", I.chains)
for i in idx:
    print(perms[i], res[4][i], res[6][i], res[8][i])
