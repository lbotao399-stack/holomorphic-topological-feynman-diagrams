"""Symbolic T = 1, H = 1 density of fan(2) (= bitriangle) with sources (l0, l1, l2) at the path
vertices (polybrute: polynomial numerators, exact), compared with htk.brute at rational points.
Output: bitriangle_density.log"""
import random
from fractions import Fraction as Fr
import sympy as sp
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "code"))
from polybrute import density_poly
from htk import graph, brute

G = graph.fan(2)
for lam in ([["l0"], [0], ["l2"]], [["l0"], ["l1"], ["l2"]]):
    R, gens, num, upow, U = density_poly(G, lam, 1, 1, extra=("l0", "l1", "l2"))
    rng = random.Random(3)
    ok = True
    for _ in range(3):
        t = [Fr(rng.randint(1, 20), rng.randint(1, 7)) for _ in range(5)]
        ls = [Fr(rng.randint(-5, 5)) for _ in range(3)]
        lamn = [[ls[0] if lam[0][0] != 0 else 0], [ls[1] if lam[1][0] != 0 else 0], [ls[2] if lam[2][0] != 0 else 0]]
        Db = brute.density(G, t, lamn, T=1)
        vals = t + ls
        nv, Uv = num(*vals), U(*vals)
        coeff = Fr(int(nv.numerator), int(nv.denominator)) / t[4] / Fr(int(Uv.numerator), int(Uv.denominator)) ** upow
        ok &= brute.Density(coeff, Fr(int(Uv.numerator), int(Uv.denominator))).same_as(Db)
    D = sp.factor(num.as_expr() / (sp.Symbol("t4") * U.as_expr() ** upow))
    print("sources", lam, ": D * U^(1/2) =", D, " (brute match at 3 points: %s)" % ok)
