"""x-degree test: compare the x-quadrature with m = 2, 3, 4, 6 at random R (exactness of the polynomial rule)."""
import numpy as np, sys
sys.path.insert(0, '..')
from graphs import K4SUB
from hepp import K4Int
rng = np.random.default_rng(1)
for name, E in K4SUB:
    deg = [0]*8
    for a, b in E: deg[a] += 1; deg[b] += 1
    d2 = [v for v in range(8) if deg[v] == 2]
    lam = np.zeros((7, 3)); lam[d2[0]] = [1, 0, 0]; lam[d2[1]] = [0, 1, 0]; lam[d2[2]] = [0, 0, 1]
    R = np.exp(rng.normal(size=(4, 6)))
    vals = {m: K4Int(8, E, [lam], m=m).f_R(R)[:, 0] for m in (2, 3, 4, 6)}
    print(name, "a =", [len(es) for _, es in K4Int(8, E, [lam], 2).chains],
          "rel(m=2 vs 6) %.1e  rel(m=3 vs 6) %.1e  rel(m=4 vs 6) %.1e" % tuple(np.max(np.abs(vals[m]-vals[6])/np.abs(vals[6])) for m in (2, 3, 4)))
