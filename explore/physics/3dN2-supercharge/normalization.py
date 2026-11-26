"""Propagator normalization: Higher.tex (3.15) super-propagator -> position-space P; closure and flux;
CDG (5.7) as printed; GKW (3.59) bitriangle vs htk.  Log: normalization.log"""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from fractions import Fraction as Fr

tau, u, x, yv = sp.symbols("tau u x y", real=True, positive=False)
z, zb = sp.symbols("z zb")
tp = sp.Symbol("tau", positive=True)
r2 = u**2 + z * zb
# (3.15): P(x;tau) = pi^{-1/2} exp(-s^2 - z*Y) dY ds,  Y = zb/tau, s = u/tau^{1/2}
# forms in (tau, zb, u): dY = dzb/tau - zb dtau/tau^2, ds = du/tau^{1/2} - u dtau/(2 tau^{3/2})
pref = sp.pi**sp.Rational(-1, 2) * sp.exp(-r2 / tp)
# dY ^ ds = (1/tau^{3/2}) dzb^du - (u/(2 tau^{5/2})) dzb^dtau - (zb/tau^{5/2}) dtau^du
C_zu = pref / tp**sp.Rational(3, 2)              # dzb ^ du
A_tz = pref * u / (2 * tp**sp.Rational(5, 2))     # dtau ^ dzb  (= -(u/2tau^{5/2}) dzb^dtau)
B_tu = -pref * zb / tp**sp.Rational(5, 2)         # dtau ^ du
# closure in (tau, zb, u) modulo dz:  d(A dtau dzb + B dtau du + C dzb du) = 0
clos = sp.simplify(sp.diff(A_tz, u) * (-1) + sp.diff(B_tu, zb) * (+1) + sp.diff(C_zu, tp))
# coefficient of dtau^dzb^du: d(A dtau dzb) = dA/du du^dtau^dzb = dA/du dtau^dzb^du ;
# d(B dtau du) = dB/dzb dzb^dtau^du = -dB/dzb dtau^dzb^du ; d(C dzb du) = dC/dtau dtau^dzb^du
clos = sp.simplify(sp.diff(A_tz, u) - sp.diff(B_tu, zb) + sp.diff(C_zu, tp))
print("(3.15) super-propagator closed in (tau, zbar, u) mod dz:", clos == 0)
# d'(A dzb + B du) = (dB/dzb - dA/du) dzb^du = dC/dtau dzb^du  =>  P := -int_0^oo (A dzb + B du) dtau has d'P = C(0+) = delta
r = sp.Symbol("r", positive=True)
Pz = -sp.integrate(A_tz.subs(r2, r**2), (tp, 0, sp.oo))
Pu = -sp.integrate(B_tu.subs(r2, r**2), (tp, 0, sp.oo))
Pz = sp.simplify(sp.integrate((-A_tz).subs(z * zb, r**2 - u**2), (tp, 0, sp.oo), conds="none"))
Pu = sp.simplify(sp.integrate((-B_tu).subs(z * zb, r**2 - u**2), (tp, 0, sp.oo), conds="none"))
print("P_GKW = Pzb dzb + Pu du with")
print("   Pzb =", Pz, "   Pu =", Pu, "   (r^2 = u^2 + |z|^2)")
# closure of P away from 0, with r = sqrt(u^2 + z zb)
R = sp.sqrt(u**2 + z * zb)
PzR, PuR = Pz.subs(r, R), Pu.subs(r, R)
print("d'P_GKW = 0 for x != 0:", sp.simplify(sp.diff(PuR, zb) - sp.diff(PzR, u)) == 0)
# flux: d'P = delta with int dVol dzb du delta = 1 (Higher (3.2) orientation), dVol = dz/(2 pi i).
# Stokes on the ball: int_ball dVol ^ d'P = -int_S2 dVol ^ P (dVol closed 1-form).  Use real coords
# z = a + i b; dVol ^ dzb ^ du = (1/2pi i)(-2i) da^db^du = -(1/pi) da db du, so (3.2) orientation is
# the one with -(1/pi) da db du negative-oriented... compute the (3.2) integral sign directly:
a, b = sp.symbols("a b", real=True)
# int pi^{-1/2} e^{-u^2-a^2-b^2} * (-(1/pi)) da db du over R^3 with standard orientation = -1
# so (3.2) = +1 requires orientation o = -1 relative to da^db^du.
o = -1
# flux through unit sphere of the 2-form F = dVol ^ P in standard orientation, times o, times -1 (Stokes)
th, ph = sp.symbols("theta phi", real=True)
U_ = sp.cos(th); A_ = sp.sin(th) * sp.cos(ph); B_ = sp.sin(th) * sp.sin(ph)
Zs = A_ + sp.I * B_; Zbs = A_ - sp.I * B_
dz_ = [sp.diff(Zs, th), sp.diff(Zs, ph)]
dzb_ = [sp.diff(Zbs, th), sp.diff(Zbs, ph)]
du_ = [sp.diff(U_, th), sp.diff(U_, ph)]
subs = {u: U_, z: Zs, zb: Zbs}
cz, cu = PzR.subs(subs), PuR.subs(subs)
# F = (1/2pi i) dz ^ (cz dzb + cu du); pull back: coefficient of dth^dph
def w2(p, q):
    return p[0] * q[1] - p[1] * q[0]
Fpull = (cz * w2(dz_, dzb_) + cu * w2(dz_, du_)) / (2 * sp.pi * sp.I)
Fpull = sp.simplify(Fpull.subs(sp.sqrt(sp.sin(th)**2 + sp.cos(th)**2), 1))
flux = sp.simplify(sp.integrate(sp.integrate(sp.simplify(sp.expand_complex(Fpull)), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
# sphere orientation (th, ph) is outward for the standard orientation of (u, a, b) with u = cos th?
# outward normal check: (u,a,b) frame: sign of det[n, d/dth, d/dph] at th = pi/2, ph = 0
n = sp.Matrix([U_, A_, B_]); J = sp.Matrix.hstack(n, sp.Matrix([sp.diff(U_, th), sp.diff(A_, th), sp.diff(B_, th)]),
                                               sp.Matrix([sp.diff(U_, ph), sp.diff(A_, ph), sp.diff(B_, ph)]))
s_or = sp.sign(sp.simplify(J.det().subs({th: sp.pi / 2, ph: 0})))
# standard orientation da^db^du = (orientation of (u,a,b)) * sign of perm (a,b,u)->(u,a,b) (+1, cyclic)
total = sp.simplify(-o * s_or * flux)
print("int dVol ^ d'P_GKW over a ball containing 0, (3.2) orientation =", total)

# CDG (5.7) as printed: P = 3/(16 pi i) |x|^{-3} [zb dt - t dzb] dz
cP_t, cP_zb = zb / R**3, -u / R**3
print("CDG (5.7) printed: d' of (zb dt - t dzb)/|x|^3 = ", sp.factor(sp.simplify(sp.diff(cP_zb, u) - sp.diff(cP_t, zb))), " dt^dzb  (nonzero: not closed)")
for bb in [1, 2, 4]:
    Rb = sp.sqrt(u**2 + bb * z * zb)
    print("   with |x|^2 = t^2 + %d|z|^2:" % bb, sp.factor(sp.simplify(sp.diff(-u / Rb**3, u) - sp.diff(zb / Rb**3, zb))))
print("closed numerator: d'[(2 zb dt - t dzb)/|x|^3] =", sp.simplify(sp.diff(-u / R**3, u) - sp.diff(2 * zb / R**3, zb)))
print("P_GKW ^ dVol in CDG form:  (1/2pi i) * (Pu du + Pzb dzb) ^ dz  -> coefficient of (zb dt - ...):",
      sp.simplify(PuR * R**3 / zb / (2 * sp.pi * sp.I)), ",", sp.simplify(PzR * R**3 / u / (2 * sp.pi * sp.I)))

# ---------------- GKW (3.59) vs htk
from htk import graph, brute
# GKW labels v1..v4 -> 0..3, root v4 = 3; edges in GKW order (12, 23, 24, 31, 43)
G = graph.Graph(4, [(0, 1), (1, 2), (1, 3), (2, 0), (3, 2)], "bitriangle-GKW")
import random
rng = random.Random(3)
print("GKW (3.59) vs htk, edge order (t12,t23,t24,t31,t43):")
for trial in range(4):
    t = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(5)]
    lam = [[Fr(rng.randint(-5, 5))] for _ in range(3)]
    D = brute.density(G, t, lam, T=1)
    t12, t23, t24, t31, t43 = t
    U = (t12*t23 + t12*t24 + t23*t24 + t23*t31 + t24*t31 + t12*t43 + t23*t43 + t31*t43)
    l1, l2, l3 = lam[0][0], lam[1][0], lam[2][0]
    # omega(3.59) = -(pi^2/4) U^{-5/2} (prod t) l1 (l1+l2+l3) dVol_t ;  prod t * dVol_t = -iota_E prod dt
    gk = Fr(1, 4) * l1 * (l1 + l2 + l3)            # coefficient of pi^2 U^{-5/2} iota_E prod dt
    # R := omega(3.59) / [(4 pi)^{-1} Omega_htk] = pi^3 * l1 (l1+l2+l3) U^{-5/2} / D ;  R^2/pi^6 rational
    if D.coeff == 0:
        print("   t=%s lam=%s : D = 0, l1(l1+l2+l3) = %s" % ([str(v) for v in t], [str(v[0]) for v in lam], l1 * (l1 + l2 + l3)))
        continue
    num = l1 * (l1 + l2 + l3)
    R2 = num**2 / (D.square() * U**5)
    sg = (1 if num > 0 else -1) * D.sign()
    print("   t=%s lam=%s : (R/pi^3)^2 = %s, sign(R) = %+d ; D^2 U^5 / (l1 l4)^2 = %s"
          % ([str(v) for v in t], [str(v[0]) for v in lam], R2, sg, D.square() * U**5 / num**2))
