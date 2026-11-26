import numpy as np
from k33 import K33Int
lam = [[[1, 0], [0, 1], [1/3, 2/5], [-4/3, -7/5], [0, 0]]]
rng = np.random.default_rng(5)
Ex = np.exp(rng.normal(size=(6, 5)) * 1.5)
ref = K33Int(lam, 32, 3)(Ex)[:, 0]
ref2 = K33Int(lam, 28, 4)(Ex)[:, 0]
print("ref(m=32,p=3) vs ref(m=28,p=4): max rel %.1e" % np.max(np.abs(ref - ref2) / np.abs(ref).max()))
for p in (2, 3, 4):
    for m in (8, 10, 12, 16):
        v = K33Int(lam, m, p)(Ex)[:, 0]
        print("p=%d m=%d max rel err %.1e" % (p, m, np.max(np.abs(v - ref) / np.abs(ref).max())))
