"""Assembly of O_1 for every T = 0 twist.

omega_1 := coefficient of theta_0..theta_H in O_1[A = sum_j theta_j x_j e^{p_j z}]
         = N^(H) sum_Gamma |Aut Gamma|^-1 sum_{s in S_{H+1}} sgn(s) K_Gamma(p_s) W_Gamma(x_s)
directed (H+1)-gon: |Aut| = H+1, W = Str_V(x_s0 .. x_sH);  undirected (self-paired): |Aut| = 2(H+1).
sgn(s) K(p_s) = c_H Delta (polygon.log (c))  =>  omega_1 = N^(H) c_H H! Delta S_H(x),  c_H H! = (-1)^H.
Cross-checks: hCS on C^3 (anomalies REPORT: c_box = -1/6, |Aut| = 8, T_1 = Str x^4/48, omega_1 = 3 Delta T_1/pi^3,
N_1 = 1/(8 pi^3)); hCS on C^5 (type-i README: -det Tr_ad X^6 / (2 (2 pi i)^5), wheel +det/(120 (2 pi i)^5)).
GW characteristic polynomial: P_{H+1}(X) = S_H(X) / ((H+1)! (2 pi i)^H)  (directed even beta-gamma: Tr_V X^{H+1}/...).
"""
from fractions import Fraction as F
from math import factorial
import sympy as sp

from common import Log

log = Log("assemble.log")
pi, I = sp.pi, sp.I
Str, Delta = sp.symbols("Str_ad Delta")


def cH(H):
    return sp.Rational((-1) ** H, factorial(H))


log("(a) combinatorics: sum_s K W / |Aut| = c_H Delta (H+1)! S / (H+1) = c_H H! Delta S, c_H H! =",
    [cH(H) * factorial(H) for H in range(1, 9)])

# hCS on C^3: S_3 = 1/2 Str_{g[1]} x^4 = -1/2 Str_ad x^4
N1 = 1 / (8 * pi ** 3)
omega_box = N1 * cH(3) * factorial(3) * Delta * (-Str / 2)
T1 = sp.Rational(-1, 6) / 8 * (-Str)          # c_box / |Aut| * Sym W, Sym W = -Str (one ribbon, orientation rule)
log("(b) H=3 hCS: omega_1 = N_1 c_3 3! Delta S_3 =", sp.simplify(omega_box), "; anomalies REPORT 3 Delta T_1/pi^3 with"
    " T_1 = Str/48:", sp.simplify(3 * Delta * (Str / 48) / pi ** 3), "; equal:",
    sp.simplify(omega_box - 3 * Delta * (Str / 48) / pi ** 3) == 0, "; T_1 from c_box/|Aut|:", T1)

N5 = -1 / (2 * pi * I) ** 5
omega_hex = N5 * cH(5) * factorial(5) * Delta * (-Str / 2)
log("(c) H=5 hCS: N^(5) c_5 =", sp.simplify(N5 * cH(5)), "(type-i: +1/(120 (2 pi i)^5));  omega_1 =",
    sp.simplify(omega_hex), "; type-i README -Delta Tr/(2 (2 pi i)^5):",
    sp.simplify(omega_hex + Delta * Str / (2 * (2 * pi * I) ** 5)) == 0)

log("(d) per twist: c_H, S_H (supertrace.log), P_{H+1} = S_H/((H+1)!(2 pi i)^H), physical I_phys = -S_H")
ROWS = [
    ("2d (2,2), chirals R", 1, "0"),
    ("2d (0,2), chirals R, Fermis F", 1, "-Tr_adj + Tr_R - Tr_F"),
    ("2d (4,0), hypers U", 1, "-2 Tr_adj + Tr_U"),
    ("2d (4,4), hypers U", 1, "0"),
    ("2d (8,8)", 1, "0"),
    ("2d (N+,0)", 1, "-(N+/2) Tr_adj"),
    ("4d N=1, chirals R", 2, "Tr_R - Tr_adj  (= Tr_R)"),
    ("4d N=2, hypers U", 2, "Tr_U (= 0)"),
    ("4d N=4", 2, "0"),
    ("6d (1,0), U//g", 3, "-Tr_adj + 1/2 Tr_U"),
    ("6d (1,0), BF R/g", 3, "-Tr_adj + Tr_R"),
    ("6d (1,1)", 3, "0"),
    ("8d N=1", 4, "-Tr_adj (= 0)"),
    ("10d hCS C^5", 5, "-1/2 Tr_adj"),
]
for name, H, S in ROWS:
    log("  %-32s H=%d  c_H = %-7s c_H H! = %-3s S_%d = %-26s P_%d = [%s] / (%d! (2 pi i)^%d)"
        % (name, H, cH(H), cH(H) * factorial(H), H, S, H + 1, S, H + 1, H))
log.close()
