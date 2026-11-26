"""Q_1 = 0 at T = 1: the one-loop graph of Q_1 (box for H = 2, D = 3; (H+2)-gon in general) has
rho_1 = 0, so its T = 1 density vanishes identically (brute force, rational t, random sources)."""
import random
from fractions import Fraction as Fr
import htk_path  # noqa
import htk
from htk import brute
rng = random.Random(3)
for H in (1, 2, 3):
    G = htk.polygon(H + 2)
    vals = []
    for _ in range(5):
        t = [Fr(rng.randint(1, 20), rng.randint(1, 7)) for _ in range(G.ne)]
        lam = [[Fr(rng.randint(-5, 5)) for _ in range(H)] for _ in range(G.nv - 1)]
        vals.append(brute.density(G, t, lam, T=1).coeff)
    print("H=%d, T=1, polygon %d (n=1, |I| = H+2): density at 5 rational points %s" % (H, H + 2, vals))
