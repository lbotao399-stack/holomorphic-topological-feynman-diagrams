"""a5 at n = 2: banana closed form (1/2) D_Gamma prod_P Gamma(H+3/2-a_P)/(a_P-1)! / Gamma(H+3/2) versus the exact
kernel (engine density, own pointwise check, own Dirichlet integration) for theta graphs built by v_b2.banana
(p = 0, q = root, chains in order), H = 1, 2, 3; also banana (2,2,2,2,1) at H = 1, n = 4."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from v_b1 import engine_exact
from v_b2 import banana, banana_formula
import htk
rng = random.Random(21)
for lens, H, n in [((1, 2, 2), 1, 2), ((1, 3, 3), 2, 2), ((2, 2, 3), 2, 2), ((1, 4, 4), 3, 2), ((2, 3, 4), 3, 2), ((3, 3, 3), 3, 2),
                   ((2, 2, 2, 2, 1), 1, 4)]:
    nv, edges, ints = banana(lens)
    Gt = htk.theta(*lens)[0] if n == 2 else None
    same = (Gt is not None and Gt.nv == nv and list(Gt.edges) == list(edges))
    out = []
    for _ in range(2):
        lam = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)]
        K, ok, dens = engine_exact(nv, edges, lam, rng, npts=2)
        Ke = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())
        f = banana_formula(lam, ints, lens, H, n)
        out.append((ok, sp.simplify(Ke / f) if f != 0 else "f=0,K=%s" % Ke))
    print("banana%s H=%d n=%d (same graph as htk.theta: %s): [pointwise ok, K / closed form] = %s" % (lens, H, n, same, out))

# comparison with the theta formula of CONVENTIONS.md section 3 (two-loop thetas at T=1):
# K = s (D_theta/2) prod_P Gamma(H+3/2-n_P) / (Gamma(H+3/2) prod_P (n_P-1)!), s = (-1)^{ceil(H/2)+(H+1)(q+1)},
# D_theta = det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]] (htk.brute.theta_D)
from htk import brute
import math
for lens, H in [((1, 2, 2), 1), ((1, 3, 3), 2), ((2, 2, 3), 2), ((1, 4, 4), 3), ((2, 3, 4), 3), ((3, 3, 3), 3)]:
    G, interior = htk.theta(*lens)
    p, q, r = lens
    lam = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(G.nv - 1)]
    K, ok, dens = engine_exact(G.nv, list(G.edges), lam, rng, npts=2)
    Ke = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** pp for pp, v in K.items())
    s = (-1) ** (math.ceil(H / 2) + (H + 1) * (q + 1))
    Dt = brute.theta_D(lam, interior)
    g = sp.Integer(1)
    for a in lens:
        g *= sp.gamma(sp.Rational(2 * H + 3 - 2 * a, 2)) / sp.factorial(a - 1)
    g /= sp.gamma(sp.Rational(2 * H + 3, 2))
    fconv = s * sp.Rational(Dt.numerator, Dt.denominator) / 2 * g
    nv, edges, ints = banana(lens)
    f5 = banana_formula(lam, ints, lens, H, 2)
    print("theta%s H=%d: K / CONVENTIONS formula = %s ; K / a5 closed form = %s ; s = %d" % (
        lens, H, sp.simplify(Ke / fconv), sp.simplify(Ke / f5), s))
