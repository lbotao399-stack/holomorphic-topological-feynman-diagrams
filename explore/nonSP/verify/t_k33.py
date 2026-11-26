import numpy as np
from k33 import K33Int
lam = [[[1, 0], [0, 1], [1/3, 2/5], [-4/3, -7/5], [0, 0]]]
rng = np.random.default_rng(5)
Ex = np.exp(rng.normal(size=(4, 5)))
for p in (1, 2, 3):
    vals = {m: K33Int(lam, m, p)(Ex)[:, 0] for m in (4, 8, 16, 24)}
    print("p=%d" % p, ["m=%d: %s" % (m, np.array2string(vals[m], precision=10)) for m in vals])
