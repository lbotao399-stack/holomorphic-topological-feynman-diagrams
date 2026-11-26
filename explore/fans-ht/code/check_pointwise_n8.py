"""P2, P3 at n = 8 (one rational point each), own brute force only."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import fanlib as F
rng = random.Random(88)
n = 8
t0 = time.time()
u, v = F.random_r(n, rng)
t = F.t_of_r(n, u, v)
J = F.jacobian_r(n, u, v)
lam = [Fr(rng.randint(-9, 9)) for _ in range(n + 1)]
Dc, U = F.density_brute(n, t, [lam], T=1)
rat, rad = F.density_11(n, u, v, lam)
assert (rat > 0) - (rat < 0) == (Dc > 0) - (Dc < 0) and rat ** 2 * rad == (Dc * abs(J)) ** 2 / U
lam_e = [Fr(3)] + [Fr(0)] * (n - 1) + [Fr(-5)]
Dc, U = F.density_brute(n, t, [lam_e], T=1)
rat, rad = F.density_11_endpoint(n, u, v, 3, -5)
assert (rat > 0) - (rat < 0) == (Dc > 0) - (Dc < 0) and rat ** 2 * rad == (Dc * abs(J)) ** 2 / U
print("n=8: P2 (general sources) and P3 (endpoint) hold at one rational point; sign(det[t;dt/dr]) = %+d; %.0f s" % (1 if J > 0 else -1, time.time() - t0))
