"""b1 as polynomial identities in the sources: exact K (engine density at numeric sources, own pointwise
check at 2 rational points, own Dirichlet integration) == claimed formula at 20 random integer source vectors
with entries in [-100, 100] per graph (Schwartz-Zippel: a nonzero quartic vanishes at a uniform random point
of {-100..100}^5 with probability <= 4/201)."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import sympy as sp
from v_b1 import GR, FORM, L, engine_exact
rng = random.Random(4242)
for k in range(5):
    t0 = time.time()
    oks, pts = 0, 0
    for _ in range(20):
        lam_i = [rng.randint(-100, 100) for _ in range(5)]
        K, ok, dens = engine_exact(6, GR[k], [[Fr(x)] for x in lam_i], rng, npts=2)
        Ke = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())
        claim = FORM[k].subs({s: x for s, x in zip(L, lam_i)})
        oks += int(ok and sp.simplify(Ke - claim) == 0)
        pts += 1
    print("G%d: exact K == claimed formula (and pointwise check passed) at %d/%d random source vectors  (%.1fs)" % (k, oks, pts, time.time() - t0), flush=True)
