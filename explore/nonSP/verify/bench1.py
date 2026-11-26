import time, numpy as np
from hepp import K4Int
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
lam = [[0, 1], [0, 0], [-1, -1], [1, 0]]
for m in (3, 4):
    I = K4Int(5, C2N3, [lam], m=m)
    for n in (3, 4, 5, 6):
        t0 = time.time()
        v, per = I.sector_integral(n)
        print("m=%d n=%d K=%.12f  -1/576=%.12f  rel=%.2e  %.1fs" % (m, n, v[0], -1/576, abs(v[0]+1/576)*576, time.time()-t0), flush=True)
