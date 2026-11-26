"""Convention check: my density_exact vs htk.brute.density (the claims are stated in htk's K normalization)."""
import sys, os, random
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from htk.graph import Graph
from htk import brute
from vk import G, density_exact, NG
sys.path.insert(0, '..')
from graphs import K4SUB, PRISM, K33   # edge lists only
rng = random.Random(11)
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
cases = [(name, 8, E, 3) for name, E in K4SUB] + [("prism", 6, PRISM, 2), ("K33", 6, K33, 2), ("c2n3", 5, C2N3, 2)]
for name, nv, E, H in cases:
    g = G(nv, E); hg = Graph(nv, E)
    ok = []
    for _ in range(3):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in E]
        lam = [[Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(H)] for _ in range(nv - 1)]
        a = density_exact(g, t, lam)
        b = brute.density(hg, t, lam).rational()
        nd = NG(g, [lam]).D(np.array([[float(x) for x in t]]))[0, 0]
        ok.append((a == b, abs(nd - float(a)) / abs(float(a))))
    print(name, ok)
