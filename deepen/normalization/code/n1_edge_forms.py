"""(a)/(b) Edge forms of commutative CS on R x C^2 in heat-kernel gauge, exact (sympy) and one numerical flux.

Generators: 0 = dtau, 1 = dx1, 2 = dy1, 3 = dx2, 4 = dy2, 5 = dt (Schwinger time).
dz_i = dx_i + i dy_i, dzb_i = dx_i - i dy_i, Omega = dz1 dz2.
D = dtau d_tau + sum dzb_i d_{zb_i},  D* = -iota_{d_tau} d_tau - 4 sum iota_{d_zb_i} d_{z_i},
[K_t]_Omega = (1/4) e^{-(tau^2+|z|^2)/4t} (4 pi t)^{-5/2} dtau dzb1 dzb2.
"""
import sys
import sympy as sp
import mpmath as mp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/code')
import extalg as ea

I = sp.I
tau, x1, y1, x2, y2 = sp.symbols('tau x1 y1 x2 y2', real=True)
t = sp.symbols('t', positive=True)
X = [tau, x1, y1, x2, y2]
z = [x1 + I * y1, x2 + I * y2]
zb = [x1 - I * y1, x2 - I * y2]
r2 = tau ** 2 + x1 ** 2 + y1 ** 2 + x2 ** 2 + y2 ** 2
dtau = ea.one_form({0: 1})
dz = [ea.one_form({1: 1, 2: I}), ea.one_form({3: 1, 4: I})]
dzb = [ea.one_form({1: 1, 2: -I}), ea.one_form({3: 1, 4: -I})]
dt = ea.one_form({5: 1})
Omega = ea.wedge(dz[0], dz[1])
vol5 = {(0, 1, 2, 3, 4): 1}

print("(1) Omega ^ dtau ^ dzb1 ^ dzb2 =", ea.wedge_all([Omega, dtau, dzb[0], dzb[1]]), " (key (0,1,2,3,4) = dtau d^4x)")

# vector fields as contraction data
iota_tau = {0: 1}
iota_zb = [{1: sp.Rational(1, 2), 2: I / 2}, {3: sp.Rational(1, 2), 4: I / 2}]
print("    iota_{d_zb_i} dzb_j, dz_j:", [[ea.interior(iota_zb[i], dzb[j]).get((), 0) for j in range(2)] for i in range(2)],
      [[ea.interior(iota_zb[i], dz[j]).get((), 0) for j in range(2)] for i in range(2)])


def d_z(f, i):
    xs, ys = (x1, y1) if i == 0 else (x2, y2)
    return (sp.diff(f, xs) - I * sp.diff(f, ys)) / 2


def d_zb(f, i):
    xs, ys = (x1, y1) if i == 0 else (x2, y2)
    return (sp.diff(f, xs) + I * sp.diff(f, ys)) / 2


def apply_coeff(a, op):
    return {k: op(v) for k, v in a.items() if op(v) != 0}


def Dstar(a):
    out = ea.scale(ea.interior(iota_tau, apply_coeff(a, lambda f: sp.diff(f, tau))), -1)
    for i in range(2):
        out = ea.add(out, ea.scale(ea.interior(iota_zb[i], apply_coeff(a, lambda f, i=i: d_z(f, i))), -4))
    return out


def Dop(a):
    """D = dtau d_tau + sum_i dzb_i d_{zb_i} acting on forms with function coefficients."""
    out = ea.wedge(dtau, apply_coeff(a, lambda f: sp.diff(f, tau)))
    for i in range(2):
        out = ea.add(out, ea.wedge(dzb[i], apply_coeff(a, lambda f, i=i: d_zb(f, i))))
    return out


def dt_op(a):
    return ea.wedge(dt, apply_coeff(a, lambda f: sp.diff(f, t)))


g = sp.exp(-r2 / (4 * t)) / (4 * sp.pi * t) ** sp.Rational(5, 2)
K = ea.scale(ea.wedge_all([dtau, dzb[0], dzb[1]]), g / 4)

# (2) normalization
OK = ea.wedge(Omega, K)
dens = sp.simplify(OK[(0, 1, 2, 3, 4)])
val = sp.simplify(dens.subs({x1: 0, y1: 0, x2: 0, y2: 0, tau: 0}) * 0 + 1)
integral = sp.integrate(sp.exp(-tau ** 2 / (4 * t)), (tau, -sp.oo, sp.oo)) ** 5 / (4 * sp.pi * t) ** sp.Rational(5, 2)
print("(2) Omega ^ [K_t] = g d^5x with g/(e^{-r^2/4t}(4 pi t)^{-5/2}) =",
      sp.simplify(dens / (sp.exp(-r2 / (4 * t)) / (4 * sp.pi * t) ** sp.Rational(5, 2))),
      "; int g d^5x =", sp.simplify(integral))

# (3) D* K = (1/t) iota_E K
iota_E = {0: tau / 2}
for i in range(2):
    for gnr, c in iota_zb[i].items():
        iota_E[gnr] = iota_E.get(gnr, 0) + zb[i] * c
DsK = Dstar(K)
rhs = ea.scale(ea.interior(iota_E, K), 1 / t)
print("(3) D*K_t = t^{-1} iota_E K_t, E = (tau/2) d_tau + zb^i d_{zb^i}:", ea.equal(DsK, rhs))
print("    iota_E(dtau dzb1 dzb2) = (tau/2) dzb1 dzb2 - dtau (zb1 dzb2 - zb2 dzb1):",
      ea.equal(ea.interior(iota_E, ea.wedge_all([dtau, dzb[0], dzb[1]])),
               ea.add(ea.scale(ea.wedge(dzb[0], dzb[1]), tau / 2),
                      ea.wedge(dtau, ea.add(ea.scale(dzb[1], zb[0]), ea.scale(dzb[0], zb[1]), -1)), -1)))

# (4) Theta = K - dt ^ D*K = (1/(4 pi^{5/2})) e^{-eta^2 - y.z} deta ^ dy1 ^ dy2
Theta = ea.add(K, ea.wedge(dt, DsK), -1)
eta = tau / (2 * sp.sqrt(t))
ys = [zb[0] / (4 * t), zb[1] / (4 * t)]


def dfun(f):
    return ea.one_form({0: sp.diff(f, tau), 1: sp.diff(f, x1), 2: sp.diff(f, y1), 3: sp.diff(f, x2),
                        4: sp.diff(f, y2), 5: sp.diff(f, t)})


expo = sp.exp(-eta ** 2 - ys[0] * z[0] - ys[1] * z[1])
Theta2 = ea.scale(ea.wedge_all([dfun(eta), dfun(ys[0]), dfun(ys[1])]), expo / (4 * sp.pi ** sp.Rational(5, 2)))
print("(4) K_t - dt ^ D*K_t = (4 pi^{5/2})^{-1} e^{-eta^2 - y.z} deta ^ d^2y, eta = tau/(2 sqrt t), y = zb/(4t):",
      ea.equal(Theta, Theta2))

# (5) (d_t + D) Theta = 0 ; H K = D D* K = -d_t K
print("(5) (d_t + D) Theta = 0:", ea.equal(ea.add(dt_op(Theta), Dop(Theta)), {}),
      "; D D* K_t = -d_t K_t:", ea.equal(Dop(DsK), ea.scale(apply_coeff(K, lambda f: sp.diff(f, t)), -1)),
      "; D* D K_t = 0:", ea.equal(Dstar(Dop(K)), {}))

# (6) propagator
r = sp.symbols('r', positive=True)
Pform = {}
for k, v in DsK.items():
    # v = (polynomial in coords) * e^{-r^2/4t} t^{-7/2} * const ; integrate over t in (0, oo)
    vv = sp.simplify(v * sp.exp(r2 / (4 * t)))
    vv = sp.simplify(vv.subs(r2, r ** 2))
    integ = sp.integrate(sp.expand(vv) * sp.exp(-r ** 2 / (4 * t)), (t, 0, sp.oo))
    Pform[k] = sp.simplify(integ)
Pclosed = ea.scale(ea.add(ea.scale(ea.wedge(dzb[0], dzb[1]), tau / 2),
                          ea.wedge(dtau, ea.add(ea.scale(dzb[1], zb[0]), ea.scale(dzb[0], zb[1]), -1)), -1),
                   3 / (16 * sp.pi ** 2 * r ** 5))
print("(6) [P]_Omega = int_0^oo dt D*K_t = 3/(16 pi^2 r^5) [ (tau/2) dzb1 dzb2 - dtau (zb1 dzb2 - zb2 dzb1) ]:",
      ea.equal(Pform, Pclosed))
# IJZ (2.12) numerator with (x, z) -> (z1, z2), t -> tau:  (1/2) t dxb dzb + xb dzb dt - zb dxb dt
ijz = ea.add(ea.add(ea.scale(ea.wedge(dzb[0], dzb[1]), tau / 2), ea.scale(ea.wedge(dzb[1], dtau), zb[0])),
             ea.scale(ea.wedge(dzb[0], dtau), zb[1]), -1)
print("    IJZ (2.12) numerator equals ours:", ea.equal(ea.scale(ijz, 3 / (16 * sp.pi ** 2 * r ** 5)), Pclosed))
Pr = {k: v.subs(r, sp.sqrt(r2)) for k, v in Pclosed.items()}
print("(7) D[P]_Omega = 0 for r > 0:", ea.equal(Dop(Pr), {}),
      "; D*[P]_Omega = 0 (gauge, IJZ (2.11)):", ea.equal(Dstar(Pr), {}))
# parity under x -> -x
sub = {tau: -tau, x1: -x1, y1: -y1, x2: -x2, y2: -y2}
print("    [P](-x) = -[P](x) (coefficients):", all(sp.simplify(v.subs(sub, simultaneous=True) + v) == 0 for v in Pr.values()))

# (8) flux of Omega ^ P through S^4 (Stokes: = int_B Omega ^ D P = int_B Omega ^ K_0)
w = ea.wedge(Omega, Pr)
V = []
for i in range(5):
    key = tuple(j for j in range(5) if j != i)
    V.append(sp.simplify((-1) ** i * w.get(key, 0)))
radial = [sp.simplify(V[i] * r2 ** sp.Rational(5, 2)) for i in range(5)]
print("(8) Omega ^ [P] = iota_V d^5x, r^5 V =", radial)
print("    Re(r^5 V) = (3/(8 pi^2)) x:", all(sp.simplify(sp.re(radial[i]) - 3 * X[i] / (8 * sp.pi ** 2)) == 0 for i in range(5)),
      "; Im(r^5 V) = (3/(8 pi^2)) (0, -y1, x1, -y2, x2) (tangent to S^4):",
      [sp.simplify(sp.im(radial[i])) for i in range(5)])
import numpy as np
Vf = sp.lambdify(X, V, 'numpy')


def flux(R, m=24):
    xs, ws = np.polynomial.legendre.leggauss(m)
    a = 0.5 * np.pi * (xs + 1); wa = 0.5 * np.pi * ws
    b = np.pi * (xs + 1); wb = np.pi * ws
    A1, A2, A3, A4 = np.meshgrid(a, a, a, b, indexing='ij')
    W = np.einsum('i,j,k,l->ijkl', wa, wa, wa, wb)
    s1, s2, s3 = np.sin(A1), np.sin(A2), np.sin(A3)
    n = [np.cos(A1), s1 * np.cos(A2), s1 * s2 * np.cos(A3), s1 * s2 * s3 * np.cos(A4), s1 * s2 * s3 * np.sin(A4)]
    v = Vf(*[R * c for c in n])
    vn = sum(v[i] * n[i] for i in range(5))
    return np.sum(W * vn * R ** 4 * s1 ** 3 * s2 ** 2 * s3)


for R in (0.7, 1.3):
    F = flux(R)
    print("    numerical flux int_{S^4_R} Omega ^ [P] (R = %.1f, 24^4-point Gauss-Legendre): %.15f %+.1e i ; |F - 1| = %.1e"
          % (R, F.real, F.imag, abs(F - 1)))
print("    exact: 8 pi^2/3 * 3/(8 pi^2) =", sp.simplify(sp.Rational(8, 3) * sp.pi ** 2 * 3 / (8 * sp.pi ** 2)))
