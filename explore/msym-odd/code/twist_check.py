"""Twist of maximal SYM reduced to R x C^H with k = 4 - H odd coordinates (H + 1 + k = 5).

Reduction of S2 Section 3: 10d complex directions i = 1..5 (S2 eq-cx); fields depend on
x^1..x^H, on tau = x_E^{2H+1} = Re x^{H+1}, and not on Im x^{H+1}, x^{H+2..5};
theta_j = dxbar^{H+1+j}.  Phi: dxbar^{H+1} -> 2 dtau on sigma-independent forms.

Checks (exact, sympy):
 (1) Omega_H ^ dxbar^1..dxbar^H ^ dtau = 4 d^{2H+1}x_E,  Omega_1 = 2i dx, Omega_2 = dx^1 dx^2,
     Omega_3 = (i/2) dx^1 dx^2 dx^3;  6d (S2 eq-measure): i dx^1dx^2dx^3 ^ dxbar^1..3 = 8 d^6x_E.
 (2) Phi o dbar_{H+1} = (dbar + d_tau) o Phi,  Phi o dbar*_{H+1} = D* o Phi,
     D* = -2 sum_a iota_{abar} d_a - (1/2) iota_tau d_tau,  [D, D*] = -(1/2) Laplacian.
 (3) chi ^ dbar chi = e^1..e^5 L_hol (S2 eq-Lhol-form) and
     int Omega ^ [int d^k theta (1/2) C2 D C2]_top = -(1/4) int d^{2H+1}x_E L_hol,  C2 = -(i/4) Phi(chi),
     hence S^phys = S_tw / hbar with hbar = g^2/4.
"""
import itertools
import time
import sympy as sp
from forms import gen, one, wedge, wedge_all, add, scale, clean, interior, coeff_apply, is_zero

I = sp.I
T0 = time.time()


def omega_const(H, T):
    if T == 0:
        return {3: I}[H]
    return {1: 2 * I, 2: sp.Integer(1), 3: I / 2}[H]


# ------------------------------------------------------------------ (1) measures
def measure(H, T):
    d = 2 * H + T
    dx = [add(gen(2 * a), gen(2 * a + 1, I)) for a in range(H)]
    dxb = [add(gen(2 * a), gen(2 * a + 1, -I)) for a in range(H)]
    form = scale(wedge_all(dx), omega_const(H, T))
    form = wedge(form, wedge_all(dxb))
    if T:
        form = wedge(form, gen(2 * H))
    return form.get(tuple(range(d)), 0), len(form)


print("(1) Omega ^ dxbar^1..dxbar^H (^ dtau) / d^Dx_E:")
for H, T in [(3, 0), (1, 1), (2, 1), (3, 1)]:
    c, nterms = measure(H, T)
    print("   H=%d T=%d  Omega = %s dx^1..dx^H :  %s  (terms: %d)" % (H, T, omega_const(H, T), c, nterms))

# ------------------------------------------------------------------ (2) reduction of dbar and dbar*
print("(2) Phi o dbar = D o Phi, Phi o dbar* = D* o Phi, [D, D*] = -(1/2) Laplacian:")
for H in (1, 2, 3):
    xs = sp.symbols('u0:%d' % (2 * H + 1), real=True)       # x_E^1..x_E^{2H+1}, tau = xs[2H]
    tau = xs[2 * H]
    db = [lambda f, a=a: (sp.diff(f, xs[2 * a]) + I * sp.diff(f, xs[2 * a + 1])) / 2 for a in range(H)]
    dh = [lambda f, a=a: (sp.diff(f, xs[2 * a]) - I * sp.diff(f, xs[2 * a + 1])) / 2 for a in range(H)]
    # 6d-type operators on (H+1)-complex forms, sigma-independent: d/dxbar^{H+1} = d/dx^{H+1} = (1/2) d_tau
    db6 = db + [lambda f: sp.diff(f, tau) / 2]
    dh6 = dh + [lambda f: sp.diff(f, tau) / 2]
    gens = list(range(H + 1))                                # dxbar^1..dxbar^{H+1}; gen H = dtau after Phi
    omega = {}
    for r in range(H + 2):
        for S in itertools.combinations(gens, r):
            omega[S] = sp.Function('f_' + ''.join(map(str, S)) if S else 'f_')(*xs)

    def Phi(p):
        return clean({k: v * (2 if H in k else 1) for k, v in p.items()})

    def dbar6(p):
        out = {}
        for a in range(H + 1):
            out = add(out, wedge(gen(a), coeff_apply(p, db6[a])))
        return out

    def D5(p):
        out = {}
        for a in range(H):
            out = add(out, wedge(gen(a), coeff_apply(p, db[a])))
        return add(out, wedge(gen(H), coeff_apply(p, lambda f: sp.diff(f, tau))))

    def dbarstar6(p):
        out = {}
        for a in range(H + 1):
            out = add(out, scale(interior(coeff_apply(p, dh6[a]), {a: 1}), -2))
        return out

    def Dstar5(p):
        out = {}
        for a in range(H):
            out = add(out, scale(interior(coeff_apply(p, dh[a]), {a: 1}), -2))
        return add(out, scale(interior(coeff_apply(p, lambda f: sp.diff(f, tau)), {H: 1}), sp.Rational(-1, 2)))

    ok1 = is_zero(add(Phi(dbar6(omega)), scale(D5(Phi(omega)), -1)))
    ok2 = is_zero(add(Phi(dbarstar6(omega)), scale(Dstar5(Phi(omega)), -1)))
    lap = lambda f: sum(sp.diff(f, x, 2) for x in xs)
    Hop = add(D5(Dstar5(omega)), Dstar5(D5(omega)))
    ok3 = is_zero(add(Hop, coeff_apply(omega, lambda f: lap(f) / 2)))
    print("   H=%d: Phi dbar = D Phi: %s;  Phi dbar* = D* Phi: %s;  [D,D*] = -Lap/2: %s" % (H, ok1, ok2, ok3))

# ------------------------------------------------------------------ (3) holomorphic term and hbar
print("(3) chi ^ dbar chi = e^1..e^5 L_hol;  int Omega ^ [int d^k theta (1/2) C2 D C2] = c * L_hol:")
# generators: e^1..e^5 -> 0..4;  chi_ij (i<j) -> 5..14;  chi_{ij;k} = dbar_k chi_ij -> 15..64
pairs = list(itertools.combinations(range(5), 2))
CHI = {p: 5 + n for n, p in enumerate(pairs)}
DCHI = {(p, k): 15 + 5 * n + k for n, p in enumerate(pairs) for k in range(5)}


def chi_c(i, j):
    """chi_ij as (generator, sign)."""
    if i == j:
        return None
    return (CHI[(i, j)], 1) if i < j else (CHI[(j, i)], -1)


def dchi_c(i, j, k):
    if i == j:
        return None
    return (DCHI[((i, j), k)], 1) if i < j else (DCHI[((j, i), k)], -1)


chi = {}
dbchi = {}
for (i, j) in pairs:
    chi = add(chi, wedge(wedge(gen(i), gen(j)), gen(CHI[(i, j)])))
    for k in range(5):
        dbchi = add(dbchi, wedge(wedge(wedge(gen(k), gen(i)), gen(j)), gen(DCHI[((i, j), k)])))
lhs = wedge(chi, dbchi)
Lhol = {}
for perm in itertools.permutations(range(5)):
    i, j, k, l, m = perm
    sgn = sp.combinatorics.Permutation(list(perm)).signature()
    a, b = chi_c(j, k), dchi_c(l, m, i)
    Lhol = add(Lhol, scale(wedge(gen(a[0]), gen(b[0])), sp.Rational(-1, 4) * sgn * a[1] * b[1]))
rhs = wedge(wedge_all([gen(i) for i in range(5)]), Lhol)
print("   10d identity chi ^ dbar chi = e^1..e^5 L_hol:", is_zero(add(lhs, scale(rhs, -1))),
      "(terms of L_hol: %d)" % len(Lhol))


def reduce_jets(p, H, T):
    """chi_{;k}: k < H kept;  k = H (tau) -> (1/2) chi_{;tau} (T = 1), stored on the same symbol
    times 1/2;  other k -> 0."""
    out = {}
    for key, c in p.items():
        fac, keep = sp.Integer(1), True
        for g in key:
            if g >= 15:
                k = (g - 15) % 5
                if k < H:
                    continue
                if T == 1 and k == H:
                    fac *= sp.Rational(1, 2)
                    continue
                keep = False
        if keep:
            out[key] = out.get(key, 0) + fac * c
    return clean(out)


for H, T in [(3, 0), (1, 1), (2, 1), (3, 1)]:
    k = 5 - H - T
    # Phi: e^{H+1} -> 2 dtau (T = 1); D acts on chi through dbar_k chi with the reduced jets
    phi_e = {g: (2 if (T == 1 and g == H) else 1) for g in range(5)}

    def Phi(p):
        out = {}
        for key, c in p.items():
            f = sp.Integer(1)
            for g in key:
                if g < 5:
                    f *= phi_e[g]
            out[key] = c * f
        return clean(out)

    C2 = scale(Phi(chi), -I / 4)
    # D Phi(chi) = Phi(dbar chi) with dbar_{H+1} -> (1/2) d_tau, dbar_{>H+1} -> 0 (reduce_jets)
    DC2 = scale(Phi(dbchi), -I / 4)
    act = scale(wedge(C2, DC2), sp.Rational(1, 2))
    act = reduce_jets(act, H, T)
    top = tuple(range(5))                       # e^1..e^{H+T} theta_1..theta_k: the theta part last
    # int d^k theta theta_1..theta_k = 1 (theta's to the right of the form generators, even block
    # moves for k = 2; for odd k the theta block is moved to the left of the integrand first)
    coeff = {}
    for key, c in act.items():
        if tuple(g for g in key if g < 5) == top:
            rest = tuple(g for g in key if g >= 5)
            coeff[rest] = coeff.get(rest, 0) + c
    # theta_1..theta_k sits at the right end of e^1..e^5; int d^k theta from the left:
    # e_form theta_block X = (-1)^{k * (H+T)} theta_block e_form X
    sgn_theta = (-1) ** (k * (H + T))
    cmeasure = 8 if T == 0 else 4
    val = scale(clean(coeff), sgn_theta * cmeasure)
    Lred = reduce_jets(Lhol, H, T)
    ratio = None
    for key, c in val.items():
        r = sp.simplify(c / Lred.get(key, sp.nan))
        ratio = r if ratio is None else ratio
        if r != ratio:
            ratio = 'not proportional'
            break
    same_support = set(val) == set(Lred)
    print("   H=%d T=%d k=%d: int Omega ^ [int d^k theta (1/2) C2 D C2] = %s * int d^Dx_E L_hol"
          " (same support: %s, terms %d)  =>  hbar = -ratio g^2 = %s g^2"
          % (H, T, k, ratio, same_support, len(val), -ratio))
print("total time %.1fs" % (time.time() - T0))
