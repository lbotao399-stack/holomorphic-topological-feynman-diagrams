"""Propagator, edge form and position-space normalization for R x C^H, T = 1 (S3 Section 4.1 and
Theorem thm-master, with the conventions of twist_check.py):

  D = dbar + d_tau,  D* = -2 sum_a iota_{abar} d_a - (1/2) iota_tau d_tau,  [D, D*] = -(1/2) Lap,
  Omega ^ dxbar^1..dxbar^H ^ dtau = 4 d^{2H+1}x_E,
  [K_t]_Omega = exp(-(x.xbar + tau^2)/2t) / (4 (2 pi t)^{H+1/2}) dxbar^1..dxbar^H dtau,
  [D* K_t]_Omega = (1/t) iota_E [K_t]_Omega,  E = xbar^a d/dxbar^a + (tau/2) d/dtau,
  Theta_e = [K_t] - dt ^ [D* K_t] = (4 pi^{H+1/2})^{-1} exp(-eta^2 - y.x) deta ^ d^H y,
      y^a = xbar^a / 2t,  eta = tau / sqrt(2t),
  [P]_Omega = Gamma(H+1/2) / (4 pi^{H+1/2} r^{2H+1}) iota_E(dxbar^1..dxbar^H dtau),  D[P] = 0 (r > 0).

Raw graph integral  G = lim int (wedge_v Omega_v) e^{sum lam.x} D wedge_e P_e |_UV  =  s_n (8 pi^{H+1/2})^{-n} K_Gamma.
"""
import itertools
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import brute, conventions, linalg
from forms import gen, one, wedge, wedge_all, add, scale, clean, interior, coeff_apply, is_zero

I = sp.I
T0 = time.time()

# ------------------------------------------------------------------ (a) heat kernel, edge form, propagator
print("(a) heat kernel, edge form, propagator (Wirtinger variables x^a, xbar^a independent):")
for H in (1, 2, 3):
    x = sp.symbols('x1:%d' % (H + 1))
    xb = sp.symbols('xb1:%d' % (H + 1))
    tau, t = sp.symbols('tau t', positive=True)
    # generators: dxbar^a -> a (0..H-1), dtau -> H, dt -> H+1
    G_TAU, G_T = H, H + 1
    r2 = sum(x[a] * xb[a] for a in range(H)) + tau ** 2
    N = 1 / (4 * (2 * sp.pi * t) ** (H + sp.Rational(1, 2)))
    K = {tuple(range(H + 1)): N * sp.exp(-r2 / (2 * t))}

    def D(p):
        out = {}
        for a in range(H):
            out = add(out, wedge(gen(a), coeff_apply(p, lambda f, a=a: sp.diff(f, xb[a]))))
        return add(out, wedge(gen(G_TAU), coeff_apply(p, lambda f: sp.diff(f, tau))))

    def Dstar(p):
        out = {}
        for a in range(H):
            out = add(out, scale(interior(coeff_apply(p, lambda f, a=a: sp.diff(f, x[a])), {a: 1}), -2))
        return add(out, scale(interior(coeff_apply(p, lambda f: sp.diff(f, tau)), {G_TAU: 1}), sp.Rational(-1, 2)))

    # normalization: int Omega ^ K = 4 * int d^Dx_E exp(-|x|^2/2t) * N; |x|^2 = x.xbar + tau^2
    u = sp.symbols('u', real=True)
    g1 = sp.integrate(sp.exp(-u ** 2 / (2 * t)), (u, -sp.oo, sp.oo))
    normK = sp.simplify(4 * N * g1 ** (2 * H + 1))
    DK = D(K)
    heat = add(coeff_apply(K, lambda f: sp.diff(f, t)), D(Dstar(K)), Dstar(D(K)))
    E = {a: xb[a] for a in range(H)}
    E[G_TAU] = tau / 2
    iE = scale(interior(K, E), 1 / t)
    ok_DsK = is_zero(add(Dstar(K), scale(iE, -1)))
    # edge form
    Theta = add(K, scale(wedge(gen(G_T), Dstar(K)), -1))
    y = [xb[a] / (2 * t) for a in range(H)]
    eta = tau / sp.sqrt(2 * t)

    def dfun(f):
        p = {}
        for a in range(H):
            p = add(p, gen(a, sp.diff(f, xb[a])))
        p = add(p, gen(G_TAU, sp.diff(f, tau)))
        return add(p, gen(G_T, sp.diff(f, t)))

    cTheta = 1 / (4 * sp.pi ** (H + sp.Rational(1, 2)))
    rhs = scale(wedge_all([dfun(ya) for ya in y] + [dfun(eta)]),
                cTheta * sp.exp(-eta ** 2 - sum(y[a] * x[a] for a in range(H))))
    ok_Theta = is_zero(add(Theta, scale(rhs, -1)))
    rhs2 = scale(wedge_all([dfun(eta)] + [dfun(ya) for ya in y]),
                 cTheta * sp.exp(-eta ** 2 - sum(y[a] * x[a] for a in range(H))))
    ok_Theta2 = is_zero(add(Theta, scale(rhs2, -(-1) ** H)))
    closed = add(coeff_apply(Theta, lambda f: 0), wedge(gen(G_T), coeff_apply(Theta, lambda f: sp.diff(f, t))), D(Theta))
    ok_closed = is_zero(closed)
    # propagator
    rr = sp.symbols('rr', positive=True)
    cP = sp.simplify(sp.integrate((1 / t) * N * sp.exp(-rr ** 2 / (2 * t)), (t, 0, sp.oo)) * rr ** (2 * H + 1))
    cP_pred = sp.gamma(H + sp.Rational(1, 2)) / (4 * sp.pi ** (H + sp.Rational(1, 2)))
    P = clean({k: cP_pred * v / sp.sqrt(r2) ** (2 * H + 1) for k, v in interior({tuple(range(H + 1)): 1}, E).items()})
    ok_DP = is_zero(D(P))
    print("   H=%d: int Omega^K = %s; DK = 0: %s; (d_t + [D,D*])K = 0: %s; D*K = iota_E K / t: %s;"
          " Theta = c exp(..) d^Hy^deta, c = 1/(4 pi^(H+1/2)): %s; = (-1)^H c exp(..) deta^d^Hy: %s; (d_t + D) Theta = 0: %s;"
          " P coefficient %s (predicted %s, equal %s); D P = 0: %s"
          % (H, normK, DK == {}, is_zero(heat), ok_DsK, ok_Theta, ok_Theta2, ok_closed, cP, cP_pred,
             sp.simplify(cP - cP_pred) == 0, ok_DP))

# ------------------------------------------------------------------ (b), (c) Gaussian integrals at rational t
print("(b) tau-Gaussian e^{-sum_e tau_e^2/(2 t_e)}: <w_e w_f> = (1/4)[delta/t - C^T M^{-1} C],  w_e = tau_e/(2 t_e);")
print("(c) holomorphic Gaussian: xbar = 2 Lambda^{-1} lambda gives y = xbar_e/(2 t_e) = Kirchhoff currents:")
rng = random.Random(11)
for name, G in [("A", htk.graph_A()), ("B", htk.graph_B()), ("theta(2,3,4)", htk.theta(2, 3, 4)[0])]:
    t = [Fr(rng.randint(1, 30), rng.randint(1, 9)) for _ in range(G.ne)]
    Bm = G.incidence()
    L = brute.laplacian(G, t)
    Linv = linalg.inverse(L)
    g = G.nv - 1
    C = G.cycle_matrix()
    M = brute.period_matrix(G, t, C)
    Minv = linalg.inverse(M)
    # quadratic form sum_e tau_e^2 / (2 t_e) = (1/2) tau^T Lambda tau  -> covariance Lambda^{-1}
    ok_q = all(sum(Bm[e][u] * Bm[e][v] / (2 * t[e]) for e in range(G.ne)) == L[u][v] / 2 for u in range(g) for v in range(g))
    ok_w = True
    for e in range(G.ne):
        for f in range(G.ne):
            cov = sum(Bm[e][u] * Linv[u][v] * Bm[f][v] for u in range(g) for v in range(g)) / (4 * t[e] * t[f])
            pred = (Fr(1) / t[e] if e == f else 0) - sum(C[a][e] * Minv[a][b] * C[b][f] for a in range(len(C)) for b in range(len(C)))
            ok_w &= cov == pred / 4
    H = 2
    lam = [[Fr(rng.randint(-5, 5)) for _ in range(H)] for _ in range(g)]
    xbar = [[2 * sum(Linv[v][u] * lam[u][i] for u in range(g)) for i in range(H)] for v in range(g)] + [[Fr(0)] * H]
    ycur = [[(xbar[a][i] - xbar[b][i]) / (2 * t[e]) for i in range(H)] for e, (a, b) in enumerate(G.edges)]
    print("   %s: quadratic form %s, <w w> %s, currents %s" % (name, ok_q, ok_w, ycur == brute.currents(G, t, lam)))

# ------------------------------------------------------------------ (d) extraction sign
print("(d) extraction sign (S3 eq-extraction, CONVENTIONS Section 4), edge factor theta^1..theta^H theta^tau,"
      " theta^a = dxbar^a_e - 2 y^a_e dt_e, theta^tau = dtau_e - w_e dt_e:")
EXT = {}
cases = [("A", htk.graph_A(), 2), ("B", htk.graph_B(), 2)]
for pqr in [(1, 4, 4), (2, 3, 4), (3, 3, 3)]:
    cases.append(("theta%s" % (pqr,), htk.theta(*pqr)[0], 3))
cases.append(("theta(1,2,2)", htk.theta(1, 2, 2)[0], 1))
for name, G, H in cases:
    Dd = H + 1
    val = conventions.extraction_sign(G, Dd, seed=3, c=[2] * H + [1], k=[2] * Dd)
    frm = conventions.extraction_sign_formula(G, Dd)
    EXT[name] = val
    print("   %s (H=%d, D=%d, |I|=%d, |V|=%d): exact %+d, formula %+d" % (name, H, Dd, G.ne, G.nv, val, frm))


# ------------------------------------------------------------------ (e) bookkeeping of signs and constants
def perm_sign(src, dst):
    pos = {s: i for i, s in enumerate(dst)}
    p = [pos[s] for s in src]
    return linalg.perm_sign(p)


def sign_pos(H, T, Kp, tau_first=False):
    """sign of (wedge_v Omega_v) ^ [prod_v dX^{a1}_v ... prod_v dX^{aD}_v] -> prod_v (Omega_v ^ F_v),
    Omega_v = dx^1_v .. dx^H_v;  tau_first = False: a = (1..H, tau), F_v = dxbar^1_v .. dxbar^H_v dtau_v
    (this note);  tau_first = True: a = (tau, 1..H), F_v = dtau_v dxbar^1_v .. dxbar^H_v (S3)."""
    src = []
    for v in range(Kp):
        src += [('x', a, v) for a in range(H)]
    order = list(range(H)) + (['tau'] if T else [])
    if tau_first:
        order = (['tau'] if T else []) + list(range(H))
    for a in order:
        for v in range(Kp):
            src.append(('xb', a, v) if a != 'tau' else ('tau', 0, v))
    dst = []
    for v in range(Kp):
        F = [('xb', a, v) for a in range(H)]
        if T:
            F = [('tau', 0, v)] + F if tau_first else F + [('tau', 0, v)]
        dst += [('x', a, v) for a in range(H)] + F
    return perm_sign(src, dst)


def total_sign(H, T, n, ext):
    Dd = H + T
    m, Kp = Dd * n + 1, (Dd - 1) * n + 1
    P_even = (Dd - 1) % 2 == 0
    s_uv = 1 if P_even else (-1) ** (m * (m - 1) // 2)
    s_fib = (-1) ** ((Dd * Kp) * (m - 1))
    s_tau = (-1) ** (T * H * n * n)          # <prod psi^tau> moved in front of prod_i prod_r psi^i
    return s_uv * ext * s_fib * sign_pos(H, T, Kp) * s_tau


def constant(H, T, n, cOmega):
    """c_Theta^m 2^{Hn} c_Omega^{K'} (2 pi)^{H K'} (2 pi)^{T K'/2} 2^{-(H + T/2) m} (1/2)^{T n/2};
    the products of t_e cancel against (prod t / U)^{H + T/2} up to U^{-H} (U^{-1/2} rho_n)^T."""
    Dd = H + T
    m, Kp = Dd * n + 1, (Dd - 1) * n + 1
    hh = H + sp.Rational(T, 2)
    cTheta = 1 / (cOmega * sp.pi ** hh)
    return sp.simplify(cTheta ** m * 2 ** (H * n) * cOmega ** Kp * (2 * sp.pi) ** (hh * Kp)
                       * 2 ** (-hh * m) * sp.Rational(1, 2) ** (sp.Rational(T * n, 2)))


print("(e) constant and sign of G = sign * const * K:")
for H, T, cO in [(3, 0, 8), (0, 1, 1), (1, 1, 4), (2, 1, 4), (3, 1, 4)]:
    row = []
    for n in (1, 2, 3, 4):
        Dd = H + T
        m, V = Dd * n + 1, (Dd - 1) * n + 2
        cD = Dd * (Dd - 1) // 2
        ext = (-1) ** ((m * (m - 1) // 2) * cD + n * Dd + (V - 1) * n * cD)
        s = total_sign(H, T, n, ext)
        row.append("n=%d: %s * %+d" % (n, constant(H, T, n, cO), s))
    print("   (H,T)=(%d,%d), c_Omega=%d: %s" % (H, T, cO, "; ".join(row)))
print("   S2 check (H,T)=(3,0): sign should be (-1)^n (-1)^{n(n+1)/2}:",
      [(-1) ** n * conventions.sign_n(n) for n in (1, 2, 3, 4)])
# S3 check: (5,1) with Omega odd and even edge forms: sign = (-1)^n (-1)^{(5n+1)5n/2} ext
s3 = []
for n in (1, 2, 3):
    m, Kp = 6 * n + 1, 5 * n + 1
    s3.append((sign_pos(5, 1, Kp, tau_first=True) * (-1) ** (m * (m - 1) // 2), (-1) ** n * (-1) ** ((5 * n + 1) * 5 * n // 2)))
print("   S3 check (H,T)=(5,1): (s_uv * s_pos, (-1)^n (-1)^{(5n+1)5n/2}) =", s3)

# ------------------------------------------------------------------ (f) H = 0 check (S3 remark after Corollary cor-odd)
Kb = htk.kernel_exact(htk.banana(3), None, T=1)
print("(f) H=0: K(banana_3) =", Kb, "; G = (2 sqrt(pi))^{-2} K =", sp.nsimplify(Kb.sympy() / (4 * sp.pi)),
      "(S3: 3 int P^2 dP = 1/4)")
print("total time %.1fs" % (time.time() - T0))
