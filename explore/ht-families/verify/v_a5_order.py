"""a5 closed form versus exact K for bananas with the chains in different orders (v_b2.banana construction):
H = 1, n = 4: (1,2,2,2,2), (2,1,2,2,2), (2,2,2,2,1); H = 2, n = 4: (3,3,1,3,3), (3,3,3,3,1); H = 1, n = 2: (2,2,1), (2,1,2)."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import sympy as sp
from v_b1 import engine_exact
from v_b2 import banana, banana_formula
rng = random.Random(23)
for lens, H, n in [((1, 2, 2, 2, 2), 1, 4), ((2, 1, 2, 2, 2), 1, 4), ((2, 2, 2, 2, 1), 1, 4),
                   ((3, 3, 1, 3, 3), 2, 4), ((3, 3, 3, 3, 1), 2, 4), ((2, 2, 1), 1, 2), ((2, 1, 2), 1, 2)]:
    nv, edges, ints = banana(lens)
    out = []
    for _ in range(2):
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)]
        K, ok, dens = engine_exact(nv, edges, lam, rng, npts=2)
        Ke = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())
        f = banana_formula(lam, ints, lens, H, n)
        out.append((ok, sp.simplify(Ke / f) if f != 0 else "f=0,K=%s" % Ke))
    print("banana%s H=%d n=%d: [pointwise ok, K / closed form] = %s" % (lens, H, n, out), flush=True)
