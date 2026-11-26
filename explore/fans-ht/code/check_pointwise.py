"""Pointwise identities at random rational points r of the S1 domain (exact arithmetic):

 (P0) r_of_t(t_of_r(r)) = r.
 (P1) own brute-force D(t) (Kirchhoff + Pfaffian, fanlib.density_brute) = htk.brute.density.
 (P2) (1,1), general sources:  g(r) := 2^{-n/2} Pi(r) Haf[L_i M_j]  =  D(t(r)) |det[t; dt/dr]|  (sign and square).
 (P3) (1,1), endpoint sources: density_11_endpoint = D(t(r)) |det[t; dt/dr]|.
 (P4) (2,0): (-1)^{n(n+1)/2} x S1 eq-1-41 = D(t(r)) |det|  (endpoint),  potential form = D |det| (general sources).
"""
import sys
import random
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import fanlib as F
from htk import fan, brute

rng = random.Random(20261121)
for n in (2, 4, 6):
    G = fan(n)
    npts = 3 if n < 6 else 2
    for p in range(npts):
        u, v = F.random_r(n, rng)
        t = F.t_of_r(n, u, v)
        assert all(x > 0 for x in t)
        uu, vv = F.r_of_t(n, [7 * x for x in t])
        assert uu == u and vv == v, "P0"
        J = F.jacobian_r(n, u, v)
        # (1,1) general sources
        lam = [Fr(rng.randint(-9, 9)) for _ in range(n + 1)]
        Dc, U = F.density_brute(n, t, [lam], T=1)
        hb = brute.density(G, t, [[x] for x in lam], T=1)
        assert hb.sign() == ((Dc > 0) - (Dc < 0)) and hb.square() == Dc ** 2 / U, "P1"
        lhs_sq, lhs_sg = (Dc * abs(J)) ** 2 / U, (Dc > 0) - (Dc < 0)
        rat, rad = F.density_11(n, u, v, lam)
        assert (rat > 0) - (rat < 0) == lhs_sg and rat ** 2 * rad == lhs_sq, ("P2", n)
        # (1,1) endpoint
        lam_e = [Fr(0)] * (n + 1)
        lam_e[0], lam_e[n] = Fr(rng.randint(1, 9)), Fr(rng.randint(-9, -1))
        Dc, U = F.density_brute(n, t, [lam_e], T=1)
        rat, rad = F.density_11_endpoint(n, u, v, lam_e[0], lam_e[n])
        assert (rat > 0) - (rat < 0) == (Dc > 0) - (Dc < 0) and rat ** 2 * rad == (Dc * abs(J)) ** 2 / U, ("P3", n)
        rat2, rad2 = F.density_11(n, u, v, lam_e)
        assert rat2 == rat and rad2 == rad
        # (2,0)
        l1 = [Fr(rng.randint(-9, 9)) for _ in range(2)]
        l2 = [Fr(rng.randint(-9, 9)) for _ in range(2)]
        lam1 = [l1[0]] + [0] * (n - 1) + [l2[0]]
        lam2 = [l1[1]] + [0] * (n - 1) + [l2[1]]
        D20, _ = F.density_brute(n, t, [lam1, lam2], T=0)
        hb = brute.density(G, t, [[x, y] for x, y in zip(lam1, lam2)], T=0).rational()
        assert hb == D20, "P1 (2,0)"
        d = l1[0] * l2[1] - l1[1] * l2[0]
        assert F.density_20_endpoint(n, u, v, d) == D20 * abs(J), ("P4 endpoint", n)
        lg1 = [Fr(rng.randint(-9, 9)) for _ in range(n + 1)]
        lg2 = [Fr(rng.randint(-9, 9)) for _ in range(n + 1)]
        D20g, _ = F.density_brute(n, t, [lg1, lg2], T=0)
        assert F.density_20_potential(n, u, v, lg1, lg2) == D20g * abs(J), ("P4 potential", n)
        print("n=%d point %d: sign(det[t;dt/dr]) = %+d; P0-P4 hold; g_11 = %s * sqrt(%s)" % (
            n, p, 1 if J > 0 else -1, rat, rad))
print("all pointwise identities hold")
