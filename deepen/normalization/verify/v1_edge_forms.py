"""V1: heat kernel, edge form, propagator on R x C^2 (verifier's own sympy code, Wirtinger coordinates).
Generators: 0 dtau, 1 dz1, 2 dz2, 3 dzb1, 4 dzb2, 5 dt.  Coordinates tau, z1, z2, zb1, zb2 independent.
D = dtau d_tau + dzb_i d_zb_i ; D* = -iota_tau d_tau - 4 iota_zb_i d_z_i."""
import sys
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/verify')
from vext import wedge, gen, interior, add, scale, clean, wedge_all

tau, z1, z2, zb1, zb2 = sp.symbols('tau z1 z2 zb1 zb2')
t = sp.Symbol('t', positive=True)
X = {0: tau, 3: zb1, 4: zb2}           # D directions
Z = {3: z1, 4: z2}                     # holomorphic partner of zb_i
r2 = tau**2 + z1 * zb1 + z2 * zb2


def simp(f):
    return sp.simplify(f)


def D(a):
    out = {}
    for k, v in a.items():
        for g, x in X.items():
            out = add(out, wedge(gen(g), {k: sp.diff(v, x)}))
    return clean(out, simp)


def dt_(a):
    out = {}
    for k, v in a.items():
        out = add(out, wedge(gen(5), {k: sp.diff(v, t)}))
    return clean(out, simp)


def Dstar(a):
    out = {}
    for k, v in a.items():
        out = add(out, scale(interior(0, {k: 1}), -sp.diff(v, tau)))
        for g, zz in Z.items():
            out = add(out, scale(interior(g, {k: 1}), -4 * sp.diff(v, zz)))
    return clean(out, simp)


def iota_E(a):
    out = {}
    for k, v in a.items():
        out = add(out, scale(interior(0, {k: 1}), v * tau / 2))
        out = add(out, scale(interior(3, {k: 1}), v * zb1))
        out = add(out, scale(interior(4, {k: 1}), v * zb2))
    return clean(out, simp)


def eq(a, b):
    d = add(a, scale(b, -1))
    return all(sp.simplify(v) == 0 for v in d.values())


print("V1  heat kernel, edge form, propagator on R x C^2")
gauss = sp.exp(-r2 / (4 * t))
K = {(0, 3, 4): sp.Rational(1, 4) * (4 * sp.pi * t) ** sp.Rational(-5, 2) * gauss}
# (a) H on functions
f = sp.Function('f')(tau, z1, z2, zb1, zb2)
Hf = add(D(Dstar({(): f})), Dstar(D({(): f})))
print("(a) [D,D*] f = -f_tautau - 4 sum f_{z zb}:",
      sp.simplify(Hf.get((), 0) - (-sp.diff(f, tau, 2) - 4 * sp.diff(f, z1, zb1) - 4 * sp.diff(f, z2, zb2))) == 0)
# (b) heat equation and D* K = t^-1 iota_E K
HK = add(D(Dstar(K)), Dstar(D(K)))
print("(b) [D,D*] K_t = -d_t K_t:", eq(HK, scale(K, -1) if False else {k: -sp.diff(v, t) for k, v in K.items()}))
print("    D* K_t = t^-1 iota_E K_t, E = (tau/2) d_tau + zb^i d_zb^i:", eq(Dstar(K), scale(iota_E(K), 1 / t)))
print("    D* K_t components:", {k: sp.factor(v / gauss) for k, v in Dstar(K).items()})
# (c) edge form
eta = tau / (2 * sp.sqrt(t))
y1, y2 = zb1 / (4 * t), zb2 / (4 * t)


def total_d(fn):
    """d fn on (tau, zb1, zb2, t) only (fn depends on these)."""
    return clean({(0,): sp.diff(fn, tau), (3,): sp.diff(fn, zb1), (4,): sp.diff(fn, zb2), (5,): sp.diff(fn, t)}, simp)


Theta_rhs = scale(wedge_all([total_d(eta), total_d(y1), total_d(y2)]),
                  1 / (4 * sp.pi ** sp.Rational(5, 2)) * sp.exp(-eta**2 - y1 * z1 - y2 * z2))
Theta_lhs = add(K, scale(wedge(gen(5), Dstar(K)), -1))
print("(c) K_t - dt ^ D*K_t = (4 pi^{5/2})^-1 e^{-eta^2 - y.z} deta ^ dy1 ^ dy2:", eq(Theta_lhs, Theta_rhs))
print("    (d_t + D) Theta = 0:", eq(add(dt_(Theta_lhs), D(Theta_lhs)), {}))
# (d) propagator
a = sp.Symbol('a', positive=True)
P = {}
for k, v in Dstar(K).items():
    integrand = sp.simplify(v.subs(r2, 4 * a * t) if False else v)
    # write integrand = c(tau,zb) t^p e^{-r2/4t}
    c = sp.simplify(integrand / gauss)
    cc = sp.simplify(c.subs(t, 1))
    p = sp.simplify(sp.log(sp.simplify(c / cc)) / sp.log(t))
    I = sp.integrate(t ** p * sp.exp(-a / t), (t, 0, sp.oo))
    P[k] = sp.simplify(cc * I.subs(a, r2 / 4))
rr = sp.sqrt(r2)
P_claim = {(3, 4): 3 / (16 * sp.pi**2 * rr**5) * tau / 2,
           (0, 4): -3 / (16 * sp.pi**2 * rr**5) * zb1,
           (0, 3): 3 / (16 * sp.pi**2 * rr**5) * zb2}
print("(d) int_0^oo D*K_t dt =", {k: sp.simplify(v * rr**5) for k, v in P.items()}, "/ r^5")
print("    = (3/(16 pi^2 r^5)) [ (tau/2) dzb1 dzb2 - dtau ^ (zb1 dzb2 - zb2 dzb1) ]:", eq(P, P_claim))
print("    D P = 0:", eq(D(P), {}), "; D* P = 0:", eq(Dstar(P), {}))
# IJZ (2.12) with (x, z, t) -> (z1, z2, tau): (3/16pi^2)((1/2) t dxb^dzb + xb dzb^dt - zb dxb^dt)/r^5
IJZ = add(add(scale(wedge(gen(3), gen(4)), tau / 2), scale(wedge(gen(4), gen(0)), zb1)), scale(wedge(gen(3), gen(0)), -zb2))
print("    IJZ (2.12) numerator, (x,z,t)->(z1,z2,tau), eps1 -> 1:", eq(P, scale(IJZ, 3 / (16 * sp.pi**2 * rr**5))))

# (e) real coordinates: generators 0 dtau, 1 dx1, 2 dy1, 3 dx2, 4 dy2
I = sp.I
cplx_to_real = {0: {(0,): 1}, 1: {(1,): 1, (2,): I}, 2: {(3,): 1, (4,): I}, 3: {(1,): 1, (2,): -I}, 4: {(3,): 1, (4,): -I}}


def to_real(a):
    out = {}
    for k, v in a.items():
        out = add(out, scale(wedge_all([cplx_to_real[g] for g in k]), v))
    return clean(out, sp.expand)


Omega = wedge(gen(1), gen(2))
print("(e) Omega ^ dtau ^ dzb1 ^ dzb2 in real coordinates:", to_real(wedge(Omega, {(0, 3, 4): 1})))
x1, yy1, x2, yy2 = sp.symbols('x1 y1 x2 y2', real=True)
real_sub = {z1: x1 + I * yy1, zb1: x1 - I * yy1, z2: x2 + I * yy2, zb2: x2 - I * yy2}
OK_ = to_real(wedge(Omega, K))
g = sp.simplify(OK_[(0, 1, 2, 3, 4)].subs(real_sub))
norm = sp.integrate(g, (tau, -sp.oo, sp.oo), (x1, -sp.oo, sp.oo), (yy1, -sp.oo, sp.oo), (x2, -sp.oo, sp.oo), (yy2, -sp.oo, sp.oo))
print("    int Omega ^ [K_t] =", sp.simplify(norm))
# flux of Omega ^ P through S^4
OP = to_real(wedge(Omega, P))
coords = [tau, x1, yy1, x2, yy2]
V = [0] * 5
for k, v in OP.items():
    missing = [i for i in range(5) if i not in k][0]
    V[missing] = (-1) ** missing * sp.simplify(v.subs(real_sub))
R2 = sum(c**2 for c in coords)
xdotV = sp.simplify(sum(c * Vi for c, Vi in zip(coords, V)))
divV = sp.simplify(sum(sp.diff(Vi, c) for c, Vi in zip(coords, V)))
print("    Omega ^ [P] = iota_V d^5x ; x.V =", sp.simplify(xdotV * R2 ** sp.Rational(3, 2)), "/ r^3 ; div V =", divV)
area = 2 * sp.pi ** sp.Rational(5, 2) / sp.gamma(sp.Rational(5, 2))
print("    flux = (x.V/r) r^4 Area(S^4) =", sp.simplify(xdotV * R2 ** sp.Rational(3, 2) * area), "; Area(S^4) =", area)

# (f) WW gauge (21),(22) literal, pulled back by x = sqrt2 tau
s = sp.Symbol('s', positive=True)
HWW = sp.exp(-(2 * (z1 * zb1 + z2 * zb2) + 2 * tau**2) / (4 * t)) / (2**3 * (sp.pi * t) ** sp.Rational(5, 2))
# form d^2zb d x = dzb1 dzb2 dx, dx = sqrt2 dtau  -> sqrt2 dtau dzb1 dzb2 (dtau moved past two 1-forms: +)
HWW_form = {(0, 3, 4): sp.sqrt(2) * HWW}
# D*_WW = -2 d_z iota_zb - d_x iota_x ; d_x = d_tau/sqrt2, iota_x dtau = 1/sqrt2 -> -(1/2) d_tau iota_tau
def DstarWW(a):
    out = {}
    for k, v in a.items():
        out = add(out, scale(interior(0, {k: 1}), -sp.Rational(1, 2) * sp.diff(v, tau)))
        for gg, zz in Z.items():
            out = add(out, scale(interior(gg, {k: 1}), -2 * sp.diff(v, zz)))
    return clean(out, simp)


ratioK = sp.simplify(HWW_form[(0, 3, 4)] / K[(0, 3, 4)].subs(t, t / 2))
print("(f) WW (22) pulled back by x = sqrt2 tau: H_WW(t) / K_{t/2} =", ratioK,
      "; int Omega ^ H_WW =", sp.simplify(ratioK))
PWW = {}
for k, v in DstarWW(HWW_form).items():
    c = sp.simplify(v / sp.exp(-r2 / (2 * t)))
    cc = sp.simplify(c.subs(t, 1))
    p = sp.simplify(sp.log(sp.simplify(c / cc)) / sp.log(t))
    Iw = sp.integrate(t ** p * sp.exp(-a / t), (t, 0, sp.oo))
    PWW[k] = sp.simplify(cc * Iw.subs(a, r2 / 2))
print("    P_WW (literal (22)) / P per component:", {k: sp.simplify(PWW[k] / P[k]) for k in P})
