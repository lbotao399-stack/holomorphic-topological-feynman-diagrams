"""Controls for the vanishing checks of v7: the same engine calls on graphs where the result must be nonzero."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
random.seed(5)
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, _ = htk.theta(*pqr)
    res = []
    for _ in range(3):
        lam = [[Fr(random.randint(-4, 4)), Fr(random.randint(-4, 4))] for _ in range(G.nv - 1)]
        t = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(G.ne)]
        res.append(htk.brute.density(G, t, lam, T=1).coeff != 0)
    print("control theta%s (two loops, T=1): brute density coefficient nonzero at 3 random points: %s" % (pqr, res))
G = htk.polygon(4)
lam3 = [[Fr(1), Fr(2), Fr(3)], [Fr(-1), Fr(3), Fr(1)], [Fr(2), Fr(-5, 3), Fr(1, 2)]]
print("control polygon(4) at (H,T) = (3,0): K =", htk.kernel_exact(G, lam3, T=0))
