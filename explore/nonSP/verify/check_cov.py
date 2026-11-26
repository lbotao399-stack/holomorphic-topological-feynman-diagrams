"""Pointwise GL(H) covariance D(t, g lam) = det(g)^n D(t, lam) (fresh solver), used for K = c Delta^3 (C^3) and
for the SL(2)-invariance of the C^2 kernels (polynomial in the brackets)."""
import sys, random
from fractions import Fraction as Fr
from vk import G, density_exact, fdet
sys.path.insert(0, '..')
from graphs import K4SUB, PRISM, K33
rng = random.Random(5)
for name, nv, E, H in [(K4SUB[0][0], 8, K4SUB[0][1], 3), (K4SUB[3][0], 8, K4SUB[3][1], 3), ("prism", 6, PRISM, 2), ("K33", 6, K33, 2)]:
    g = G(nv, E)
    out = []
    for _ in range(2):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in E]
        lam = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(H)] for _ in range(nv - 1)]
        gm = [[Fr(rng.randint(-3, 3), rng.randint(1, 2)) for _ in range(H)] for _ in range(H)]
        glam = [[sum(gm[i][j] * l[j] for j in range(H)) for i in range(H)] for l in lam]
        out.append(density_exact(g, t, glam) / density_exact(g, t, lam) == fdet(gm) ** g.n)
    print(name, "D(t, g lam) = det(g)^n D(t, lam):", out)
