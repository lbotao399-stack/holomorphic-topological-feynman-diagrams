"""Independent checks of the claims 'twist-5d', 'action-hbar' and 'propagator-T1' (sympy, exact).

Grassmann elements: dict {sorted tuple of odd generator ids: sympy coefficient}.
(1) On sigma-independent forms on C^{H+1} (sigma = Im x^{H+1}, tau = Re x^{H+1}):
    Phi(dxbar^{H+1}) = 2 dtau:  Phi dbar = (dbar_H + d_tau) Phi,  Phi dbar* = D* Phi,
    D* = -2 sum_{a<=H} iota_abar d_a - (1/2) iota_tau d_tau,  [D, D*] = -(1/2) Laplacian_{R^{2H+1}};
    with the pull-back dxbar^{H+1} -> dtau instead: Phi' dbar = (dbar_H + (1/2) d_tau) Phi'.
(2) Omega ^ dxbar^1..dxbar^H (^ dtau) in terms of d^D x_E.
(3) chi = sum_{i<j} e^i e^j chi_ij (10d), dbar chi = sum e^i e^j e^k d_k chi_ij,
    chi ^ dbar chi = e^1..e^5 L_hol, L_hol = -1/4 eps^{ijklm} chi_jk d_i chi_lm;
    reduced: int Omega ^ [int d^k theta (1/2) C2 D C2] = c * int d^D x_E L_hol^red,  C2 = -(i/4) Phi(chi).
(4) heat kernel, edge form, propagator for T = 1 (claim 'propagator-T1').
"""
import itertools
import sympy as sp

I = sp.I


def gmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            if set(ka) & set(kb):
                continue
            seq = list(ka + kb)
            s = 1
            for i in range(len(seq)):
                for j in range(i + 1, len(seq)):
                    if seq[i] > seq[j]:
                        s = -s
            key = tuple(sorted(seq))
            out[key] = out.get(key, 0) + s * va * vb
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def gadd(*xs):
    out = {}
    for x in xs:
        for k, v in x.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def gsc(a, c):
    return {k: c * v for k, v in a.items()}


def g(i, c=1):
    return {(i,): sp.sympify(c)}


def capply(a, op):
    return {k: op(v) for k, v in a.items()}


def iota(a, gen):
    """interior product with the vector dual to generator gen (acts from the left)."""
    out = {}
    for k, v in a.items():
        if gen in k:
            p = k.index(gen)
            key = k[:p] + k[p + 1:]
            out[key] = out.get(key, 0) + (-1) ** p * v
    return out


def is_zero(a):
    return all(sp.simplify(v) == 0 for v in a.values())


print("(1) reduction of dbar, dbar* along sigma = Im x^{H+1}:")
for H in (1, 2, 3):
    X = sp.symbols('X1:%d' % (2 * H + 1), real=True)        # x_E^1..x_E^{2H}
    tau = sp.Symbol('tau', real=True)
    db = [lambda f, a=a: (sp.diff(f, X[2 * a]) + I * sp.diff(f, X[2 * a + 1])) / 2 for a in range(H)]
    dh = [lambda f, a=a: (sp.diff(f, X[2 * a]) - I * sp.diff(f, X[2 * a + 1])) / 2 for a in range(H)]
    # on sigma-independent functions: d/dxbar^{H+1} = d/dx^{H+1} = (1/2) d_tau
    db6 = db + [lambda f: sp.diff(f, tau) / 2]
    dh6 = dh + [lambda f: sp.diff(f, tau) / 2]
    w = {}
    for r in range(H + 2):
        for S in itertools.combinations(range(H + 1), r):
            w[S] = sp.Function('w' + ''.join(map(str, S)))(*X, tau)

    def Phi(a, c=2):
        return {k: v * (c if H in k else 1) for k, v in a.items()}

    def dbar6(a):
        return gadd(*[gmul(g(i), capply(a, db6[i])) for i in range(H + 1)])

    def dbs6(a):
        return gadd(*[gsc(iota(capply(a, dh6[i]), i), -2) for i in range(H + 1)])

    def Dop(a, ctau=1):
        return gadd(*[gmul(g(i), capply(a, db[i])) for i in range(H)], gmul(g(H), capply(a, lambda f: ctau * sp.diff(f, tau))))

    def Dst(a):
        return gadd(*[gsc(iota(capply(a, dh[i]), i), -2) for i in range(H)],
                    gsc(iota(capply(a, lambda f: sp.diff(f, tau)), H), sp.Rational(-1, 2)))

    r1 = is_zero(gadd(Phi(dbar6(w)), gsc(Dop(Phi(w)), -1)))
    r2 = is_zero(gadd(Phi(dbs6(w)), gsc(Dst(Phi(w)), -1)))
    lap = lambda f: sum(sp.diff(f, v, 2) for v in X) + sp.diff(f, tau, 2)
    r3 = is_zero(gadd(Dop(Dst(w)), Dst(Dop(w)), capply(w, lambda f: lap(f) / 2)))
    r4 = is_zero(gadd(Phi(dbar6(w), 1), gsc(Dop(Phi(w, 1), sp.Rational(1, 2)), -1)))
    print("   H=%d: Phi dbar = (dbar+d_tau) Phi: %s; Phi dbar* = D* Phi: %s; [D,D*] = -Lap/2: %s; pull-back: Phi' dbar = (dbar + d_tau/2) Phi': %s"
          % (H, r1, r2, r3, r4))

print("(2) measures:")
for H, T, cO in [(3, 0, I), (1, 1, 2 * I), (2, 1, 1), (3, 1, I / 2)]:
    # real generators dx_E^1..dx_E^{2H+T}: ids 0..; dx^a = e_{2a} + i e_{2a+1}
    dx = [gadd(g(2 * a), g(2 * a + 1, I)) for a in range(H)]
    dxb = [gadd(g(2 * a), g(2 * a + 1, -I)) for a in range(H)]
    f = {(): cO}
    for a in range(H):
        f = gmul(f, dx[a])
    for a in range(H):
        f = gmul(f, dxb[a])
    if T:
        f = gmul(f, g(2 * H))
    print("   H=%d T=%d Omega = %s dx^1..dx^H : Omega ^ dxbar^1..dxbar^H %s= %s d^%dx_E"
          % (H, T, cO, "^ dtau " if T else "", sp.simplify(f[tuple(range(2 * H + T))]), 2 * H + T))

print("(3) holomorphic term:")
# generators: e^1..e^5 -> 0..4, chi_ij (i<j) -> 10 + 5i + j, d_k chi_ij -> 100 + 25 k + 5 i + j
pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
CH = lambda i, j: 10 + 5 * i + j
DCH = lambda k, i, j: 100 + 25 * k + 5 * i + j


def chi_comp(i, j):
    if i == j:
        return {}
    return g(CH(i, j)) if i < j else g(CH(j, i), -1)


def dchi_comp(k, i, j):
    if i == j:
        return {}
    return g(DCH(k, i, j)) if i < j else g(DCH(k, j, i), -1)


chi = gadd(*[gmul(gmul(g(i), g(j)), g(CH(i, j))) for (i, j) in pairs])
dchi = gadd(*[gmul(gmul(gmul(g(i), g(j)), g(k)), g(DCH(k, i, j))) for (i, j) in pairs for k in range(5)])
lhs = gmul(chi, dchi)
Lhol = {}
for perm in itertools.permutations(range(5)):
    i, j, k, l, m = perm
    sgn = sp.combinatorics.Permutation(list(perm)).signature()
    Lhol = gadd(Lhol, gsc(gmul(chi_comp(j, k), dchi_comp(i, l, m)), sp.Rational(-1, 4) * sgn))
rhs = gmul({(0, 1, 2, 3, 4): sp.Integer(1)}, Lhol)
print("   10d: chi ^ dbar chi = e^1..e^5 L_hol: %s" % is_zero(gadd(lhs, gsc(rhs, -1))))


def reduce_(a, H, T):
    """d_k chi with k < H kept; k = H -> (1/2) d_tau (T = 1, same symbol times 1/2); k > H(+T-1) -> 0."""
    out = {}
    for key, v in a.items():
        fac, keep = sp.Integer(1), True
        for x in key:
            if x >= 100:
                k = (x - 100) // 25
                if k < H:
                    continue
                if T == 1 and k == H:
                    fac /= 2
                    continue
                keep = False
        if keep:
            out[key] = out.get(key, 0) + fac * v
    return {k: v for k, v in out.items() if v != 0}


for H, T, cmeas in [(3, 0, 8), (1, 1, 4), (2, 1, 4), (3, 1, 4)]:
    k = 5 - H - T
    phic = 2 if T else 1
    # Phi(e^{H+1}) = 2 dtau when T = 1: overall factor 2 on the top form; e^{H+T+1..5} = theta_1..theta_k
    act = gsc(reduce_(gmul(chi, dchi), H, T), sp.Rational(1, 2) * (-I / 4) ** 2 * phic)
    # top component e^1..e^5 = (e-form part) (theta block); int d^k theta from the left: move theta block
    # (k odd generators) to the left of the H+T e-generators: sign (-1)^{k (H+T)}
    top = {key[5:]: v for key, v in act.items() if key[:5] == (0, 1, 2, 3, 4)}
    sgn = (-1) ** (k * (H + T))
    val = {kk: sgn * cmeas * v for kk, v in top.items()}
    Lred = reduce_(Lhol, H, T)
    ratios = {sp.simplify(val.get(kk, 0) / Lred[kk]) for kk in Lred}
    print("   (H,T,k)=(%d,%d,%d): int Omega ^ [int d^k theta (1/2) C2 D C2] / int d^Dx L_hol^red = %s (same support: %s)"
          % (H, T, k, ratios, set(val) == set(Lred)))

print("(4) T = 1 heat kernel, edge form, propagator:")
for H in (1, 2, 3):
    x = sp.symbols('x1:%d' % (H + 1)); xb = sp.symbols('xb1:%d' % (H + 1))
    tau = sp.Symbol('tau', real=True); t = sp.Symbol('t', positive=True)
    GT, Gt = H, H + 1               # dtau, dt generators; dxbar^a -> a
    r2 = sum(x[a] * xb[a] for a in range(H)) + tau ** 2
    Kc = sp.exp(-r2 / (2 * t)) / (4 * (2 * sp.pi * t) ** (H + sp.Rational(1, 2)))
    K = {tuple(range(H + 1)): Kc}

    def D(a):
        return gadd(*[gmul(g(i), capply(a, lambda f, i=i: sp.diff(f, xb[i]))) for i in range(H)],
                    gmul(g(GT), capply(a, lambda f: sp.diff(f, tau))))

    def Ds(a):
        return gadd(*[gsc(iota(capply(a, lambda f, i=i: sp.diff(f, x[i])), i), -2) for i in range(H)],
                    gsc(iota(capply(a, lambda f: sp.diff(f, tau)), GT), sp.Rational(-1, 2)))

    # normalization: Omega ^ dxbar^1..dxbar^H dtau = 4 d^{2H+1}x_E and |x|^2 = x.xbar + tau^2
    u = sp.Symbol('u', real=True)
    norm = sp.simplify(4 * Kc.subs({x[a]: 0 for a in range(H)}).subs(tau, 0) * sp.integrate(sp.exp(-u ** 2 / (2 * t)), (u, -sp.oo, sp.oo)) ** (2 * H + 1))
    heat = is_zero(gadd(capply(K, lambda f: sp.diff(f, t)), D(Ds(K)), Ds(D(K))))
    iE = gadd(*[gsc(iota(K, i), xb[i]) for i in range(H)], gsc(iota(K, GT), tau / 2))
    okDs = is_zero(gadd(Ds(K), gsc(iE, -1 / t)))
    Theta = gadd(K, gsc(gmul(g(Gt), Ds(K)), -1))
    y = [xb[a] / (2 * t) for a in range(H)]
    eta = tau / sp.sqrt(2 * t)

    def dd(f):
        return gadd(*[g(i, sp.diff(f, xb[i])) for i in range(H)], g(GT, sp.diff(f, tau)), g(Gt, sp.diff(f, t)))
    form = {(): sp.exp(-eta ** 2 - sum(y[a] * x[a] for a in range(H))) / (4 * sp.pi ** (H + sp.Rational(1, 2)))}
    for a in range(H):
        form = gmul(form, dd(y[a]))
    form = gmul(form, dd(eta))
    okTheta = is_zero(gadd(Theta, gsc(form, -1)))
    closed = is_zero(gadd(gmul(g(Gt), capply(Theta, lambda f: sp.diff(f, t))), D(Theta)))
    rr = sp.Symbol('r', positive=True)
    cP = sp.simplify(sp.integrate(Kc.subs({x[a]: 0 for a in range(H)}).subs(tau, rr) / t, (t, 0, sp.oo)) * rr ** (2 * H + 1))
    cPc = sp.gamma(H + sp.Rational(1, 2)) / (4 * sp.pi ** (H + sp.Rational(1, 2)))
    Pform = {}
    vol = {tuple(range(H + 1)): sp.Integer(1)}
    Pform = gadd(*[gsc(iota(vol, i), xb[i]) for i in range(H)], gsc(iota(vol, GT), tau / 2))
    Pform = gsc(Pform, cPc / sp.sqrt(r2) ** (2 * H + 1))
    okDP = is_zero(D(Pform))
    print("   H=%d: int Omega^K = %s; (d_t + [D,D*])K = 0: %s; D*K = iota_E K/t: %s; K - dt^D*K = (4 pi^(H+1/2))^-1 e^(-eta^2-y.x) d^Hy^deta: %s;"
          " (d_t + D)Theta = 0: %s; P coefficient %s (claim Gamma(H+1/2)/(4 pi^(H+1/2)) = %s: %s); D P = 0: %s; P degree %d"
          % (H, norm, heat, okDs, okTheta, closed, cP, cPc, sp.simplify(cP - cPc) == 0, okDP, H))
