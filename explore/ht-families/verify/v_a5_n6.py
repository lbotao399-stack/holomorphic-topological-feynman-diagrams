"""a5 at n = 6, (H,T) = (1,1): banana (2,2,2,2,2,2,1) in the explorer's n=6 list (graph 1: poles 5, 6,
root 7 interior).  Engine density (numeric sources) checked pointwise against the own exact brute force at 2
rational points; own Dirichlet integration; compared with (5!!/2^3) D_Gamma prod_P Gamma(H+3/2-a_P)/(a_P-1)! /
Gamma(H+7/2) = D_Gamma pi^3/7 up to sign (|D_Gamma| = |product of the six interior sources|)."""
import sys, random, time, math
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from v_b1 import engine_exact
edges = [(0, 5), (0, 6), (1, 5), (1, 6), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
rng = random.Random(66)
for lam_i in ([3, 1, -3, 1, 2, -3, 1], [2, -1, 3, 1, -2, 1, 2]):
    lam = [[Fr(x)] for x in lam_i]
    t0 = time.time()
    K, ok, dens = engine_exact(8, edges, lam, rng, npts=2)
    Ke = sum(sp.Rational(v.numerator, v.denominator) * sp.pi ** p for p, v in K.items())
    inter = [lam_i[v] for v in range(5)] + [-sum(lam_i)]
    Dg = math.prod(inter)
    print("lambda %s: pointwise %s ; K = %s ; |D_Gamma| pi^3/7 = %s ; K/(D_Gamma pi^3/7) = %s  (%.1fs)" % (
        lam_i, ok, Ke, abs(Dg) * sp.pi ** 3 / 7, sp.simplify(Ke / (Dg * sp.pi ** 3 / 7)), time.time() - t0), flush=True)
