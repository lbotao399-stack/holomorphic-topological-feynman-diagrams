"""V2: exact pointwise position-space identity (verifier's own derivation and code).

For a graph Gamma (root = last vertex at 0), edge forms Theta_e (edge order), and plane waves e^{sum_v lam_v . z_v}:
    Phi(t) := sum_j dt_0..^j..dt_{N-1} int_{X^{m}} (wedge_{v<root} Omega_v) e^{lam.z} a_j(t),
where wedge_e Theta_e = sum_j dt_0..^j..dt_{N-1} ^ a_j + (other dt degrees)  (dt's moved to the front).
By the fiber-integration identity (re-derived in VERIFY.txt), the raw UV integral is
    G = lim int (wedge Omega_v) e^{lam.z} D wedge_e P_e |_UV = int_{RP^{N-1}_{>0}, iota_E dt^N} Phi .
Claim tested: Phi = c * Omega_Gamma^{htk},  Omega_Gamma^{htk} = D(t) iota_E(dt_0..dt_{N-1}), D from htk.brute.density.

Every term is computed by explicit expansion of wedge_e prod_a (dX^a_e - c^a_e dt_e) (no extraction formula):
  T-type   : deta_e = (2 sqrt t_e)^-1 (dtau_e - w_e dt_e),  w_e = tau_e/(2 t_e),  exponent -tau_e^2/(4 t_e)
  hol. i   : dy^i_e = (alpha t_e)^-1 (dzb^i_e - u^i_e dt_e), u^i_e = zb^i_e/t_e, exponent -zb_e z_e/(alpha t_e)
Gaussian integrals over real Lebesgue measure (exact):
  int d^m tau e^{-tau Lam tau/4} P(tau) = (4 pi)^{m/2} det(Lam)^{-1/2} E[P],  E[tau_v tau_w] = 2 Lam^{-1}_{vw}
  int d^{2m}x e^{-zb Lam z/alpha + lam z} f(zb) = (alpha pi)^m det(Lam)^{-1} f(alpha Lam^{-1} lam)
The constant mu: Omega_v ^ (wedge_a dX^a_v) = mu d^D x_v, computed from dz = dx + i dy, dzb = dx - i dy.
"""
import sys
import itertools
import random
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/verify')
import htk
from htk import Graph
from vext import perm_sign, wedge_all, add, scale, clean

rng = random.Random(77)


def mat_det(A):
    A = [row[:] for row in A]
    n = len(A)
    d = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            A[i], A[p] = A[p], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            if f:
                for c in range(i, n):
                    A[r][c] -= f * A[i][c]
    return d


def mat_inv(A):
    n = len(A)
    M = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [x / piv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [a - f * b for a, b in zip(M[r], M[i])]
    return [row[n:] for row in M]


def matchings(lst):
    if not lst:
        yield []
        return
    a = lst[0]
    for k in range(1, len(lst)):
        rest = lst[1:k] + lst[k + 1:]
        for m in matchings(rest):
            yield [(a, lst[k])] + m


# ---------------------------------------------------------------- theory data
I = sp.I
THEORIES = {
    # name: (T, H, alpha, Omega letters, Omega factor, K-coefficient(t) without exponential, sympy)
    '(2,1)': (1, 2, 4, ['z1', 'z2'], 1, lambda t: sp.Rational(1, 4) * (4 * sp.pi * t) ** sp.Rational(-5, 2)),
    '(3,0)': (0, 3, 2, ['z1', 'z2', 'z3'], I, lambda t: 1 / (8 * (2 * sp.pi * t) ** 3)),
    '(0,1)': (1, 0, None, [], 1, lambda t: (4 * sp.pi * t) ** sp.Rational(-1, 2)),
    '(1,0)': (0, 1, 4, ['z1'], 1, lambda t: 1 / ((-2 * I) * 4 * sp.pi * t)),
}


def mu_constant(T, H, Om_letters, Om_factor):
    """Omega ^ (dtau) ^ dzb^1 .. dzb^H = mu * (dtau dx1 dy1 .. dxH dyH)."""
    # real generators: 0 = dtau (if T), then (dx_i, dy_i)
    base = 1 if T else 0
    def dz(i, conj):
        return {(base + 2 * (i - 1),): 1, (base + 2 * (i - 1) + 1,): (-I if conj else I)}
    forms = [dz(int(l[1:]), False) for l in Om_letters]
    if T:
        forms.append({(0,): 1})
    forms += [dz(i, True) for i in range(1, H + 1)]
    w = wedge_all(forms)
    w = clean({k: sp.expand(v) for k, v in w.items()})
    assert len(w) == 1
    (k, v), = w.items()
    assert k == tuple(range(base + 2 * H))
    return sp.simplify(Om_factor * v)


def heat_edge_c0(T, H, alpha, Kcoef):
    t = sp.Symbol('t', positive=True)
    norm = (1 / (2 * sp.sqrt(t))) ** T * (1 / (alpha * t)) ** H if H else (1 / (2 * sp.sqrt(t))) ** T
    c0 = sp.simplify(Kcoef(t) / norm)
    assert sp.diff(c0, t) == 0
    return c0


def pointwise(G, theory, tq, lam):
    """tq: list of positive rationals q_e, t_e = q_e^2.  lam: list over non-root vertices of H-vectors (Fractions).
    Returns (list over j of Phi_j * sqrt(U) as sympy number), U, htk density coeff D*sqrt(U), N."""
    T, H, alpha, Om, Omf, Kc = THEORIES[theory]
    N, nv = G.ne, G.nv
    m = nv - 1
    n = N - m
    D = T + H
    assert D * n == N - 1, (D, n, N)
    t = [q * q for q in tq]
    B = [[Fr(int(v == a) - int(v == b)) for v in range(m)] for (a, b) in G.edges]   # B_ev, v < root
    Lam = [[sum(B[e][v] * B[e][w] / t[e] for e in range(N)) for w in range(m)] for v in range(m)]
    detLam = mat_det(Lam)
    Li = mat_inv(Lam)
    prod_t = Fr(1)
    for x in t:
        prod_t *= x
    U = detLam * prod_t
    # currents y = t^-1 B Lam^-1 lam  (per holomorphic component)
    yv = []
    for i in range(H):
        phi = [sum(Li[v][w] * lam[w][i] for w in range(m)) for v in range(m)]
        yv.append([sum(B[e][v] * phi[v] for v in range(m)) / t[e] for e in range(N)])
    # tau covariance of w_e = tau_e/(2 t_e):  E[w_e w_f] = 2 (B Li B^T)_ef / (4 t_e t_f)
    BLB = [[sum(B[e][v] * Li[v][w] * B[f][w] for v in range(m) for w in range(m)) for f in range(N)] for e in range(N)]
    Eww = [[2 * BLB[e][f] / (4 * t[e] * t[f]) for f in range(N)] for e in range(N)]
    types = (['T'] if T else []) + ['h%d' % i for i in range(H)]
    # minors det B[R, :]
    minor = {}
    for R in itertools.combinations(range(N), m):
        minor[R] = mat_det([B[e] for e in R])
    # s2: (wedge_v Omega_v) ^ (wedge_a wedge_v dX^a_v)  ->  prod_v (Omega_v ^ wedge_a dX^a_v)
    letters_now = [('O', v, l) for v in range(m) for l in range(len(Om))] + [('X', a, v) for a in range(D) for v in range(m)]
    letters_ref = [x for v in range(m) for x in ([('O', v, l) for l in range(len(Om))] + [('X', a, v) for a in range(D)])]
    pos = {x: i for i, x in enumerate(letters_ref)}
    s2 = perm_sign([pos[x] for x in letters_now])
    Phi = []
    for j in range(N):
        others = [e for e in range(N) if e != j]
        tot = Fr(0)
        # assignment: which type each other edge gives its dt to; n per type
        for combo in itertools.product(range(D), repeat=len(others)):
            if any(combo.count(a) != n for a in range(D)):
                continue
            ass = dict(zip(others, combo))
            word = []
            for e in range(N):
                for a in range(D):
                    if e != j and ass[e] == a:
                        word.append(('t', e))
                    else:
                        word.append(('X', a, e))
            dts = [x for x in word if x[0] == 't']
            Xs = sorted([x for x in word if x[0] == 'X'], key=lambda x: (x[1], x[2]))
            target = dts + Xs
            ps = {x: i for i, x in enumerate(target)}
            s1 = perm_sign([ps[x] for x in word])
            coef = Fr(s1)
            for a in range(D):
                R = tuple(e for e in range(N) if not (e != j and ass[e] == a))
                coef *= minor[R]
            if coef == 0:
                continue
            coef *= (-1) ** (N - 1)
            # moments
            for a in range(D):
                Sa = [e for e in others if ass[e] == a]
                if types[a] == 'T':
                    mom = Fr(0)
                    for mt in matchings(Sa):
                        pr = Fr(1)
                        for (e, f) in mt:
                            pr *= Eww[e][f]
                        mom += pr
                    coef *= mom
                else:
                    i = int(types[a][1:])
                    for e in Sa:
                        coef *= alpha * yv[i][e]      # u = zb_e/t_e -> alpha y_e
            tot += coef * s2
        Phi.append(tot)
    # constants: prod_e (c0 * norm(t_e)) * mu^m * Gaussian normalisations ; sqrt(U) factored out for T = 1
    c0 = heat_edge_c0(T, H, alpha, Kc)
    mu = mu_constant(T, H, Om, Omf)
    const = c0 ** N * mu ** m
    rat = Fr(1)
    for e in range(N):
        if T:
            rat *= Fr(1, 2) / tq[e]                    # 1/(2 sqrt t_e)
        rat *= Fr(1) / (alpha * t[e]) ** H if H else 1
    if T:
        const *= (4 * sp.pi) ** sp.Rational(m, 2)      # tau Gaussian ; det(Lam)^{-1/2} prod 1/sqrt(t) = U^{-1/2}
        rat *= 1                                        # det(Lam)^{-1/2} prod_e t_e^{-1/2} = U^{-1/2}: handled by caller
        # the factor prod_e 1/sqrt(t_e) is included in 1/(2 sqrt t_e) above, so divide it out and use U^{-1/2}:
        pq = Fr(1)
        for q in tq:
            pq *= q
        rat *= pq                                       # undo prod 1/sqrt(t_e)
        # remaining: det(Lam)^{-1/2} prod_e t_e^{-1/2} = U^{-1/2}
    if H:
        const *= ((alpha * sp.pi) ** m) ** H
        rat *= Fr(1) / detLam ** H
    # Phi_j * sqrt(U)^T
    vals = [const * sp.Rational(rat.numerator, rat.denominator) * sp.Rational(p.numerator, p.denominator) for p in Phi]
    dens = htk.brute.density(G, t, lam if H else None, T=T)
    assert dens.expo == 0
    if T:
        assert dens.radicand == U, (dens.radicand, U)
    return vals, U, dens.coeff, N


def run(G, theory, label, npts=2, H=None):
    T, Hh, alpha, Om, Omf, Kc = THEORIES[theory]
    out = []
    for _ in range(npts):
        tq = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(G.ne)]
        lam = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(Hh)] for _ in range(G.nv - 1)]
        vals, U, Dc, N = pointwise(G, theory, tq, lam)
        t = [q * q for q in tq]
        ratios = []
        for j in range(N):
            ref = Dc * (-1) ** j * t[j]
            ratios.append(sp.simplify(vals[j] / sp.Rational(ref.numerator, ref.denominator)) if ref != 0 else None)
        out.append(ratios)
    allr = set(str(r) for rr in out for r in rr)
    print("  %-22s %s : Phi_j / (D(t) (-1)^j t_j), all j, %d random (t, lam): %s" % (label, theory, npts, sorted(allr)))
    return out


if __name__ == '__main__':
    print("V2  pointwise position-space identity  Phi = c * D(t) iota_E vol  (exact)")
    for th in THEORIES:
        T, H, alpha, Om, Omf, Kc = THEORIES[th]
        print("  theory %s: mu = %s, c0 = %s" % (th, mu_constant(T, H, Om, Omf), heat_edge_c0(T, H, alpha, Kc)))
    print("(2,1) thetas, engine labelling:")
    res = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, inter = htk.theta(*pqr)
        res[pqr] = run(G, '(2,1)', 'theta%s' % (pqr,), npts=3)
    print("    claim b4: c = (-1)^{n(n-1)/2} (8 pi^{5/2})^{-n} at n=2 =", sp.Rational(-1, 64) / sp.pi**5,
          "; c2's varsigma=+1 value =", sp.Rational(1, 64) / sp.pi**5)
    print("(2,1) S2 / msym-odd graphs:")
    run(htk.graph_A(), '(2,1)', 'graph_A', npts=2)
    run(htk.graph_B(), '(2,1)', 'graph_B', npts=2)
    print("(3,0) S2 graphs:")
    run(htk.box(), '(3,0)', 'box (n=1)', npts=3)
    run(htk.graph_A(), '(3,0)', 'graph_A (n=2)', npts=2)
    run(htk.graph_B(), '(3,0)', 'graph_B (n=2)', npts=2)
    print("    I^S2 = (-1)^{n(n+1)/2} K: n=1: I=-K, n=2: I=-K.  S2 eq-B-norm: G = (-1/8pi^3)^n I  ->  c_K = +1/(8pi^3) (n=1), -1/(64pi^6) (n=2)")
    print("    explorer c2: G = +(8pi^3)^{-n} I  ->  c_K = -1/(8pi^3) (n=1), -1/(64pi^6) (n=2)")
    print("calibration graphs:")
    run(htk.banana(3), '(0,1)', 'banana(3) (n=2)', npts=2)
    run(htk.banana(5), '(0,1)', 'banana(5) (n=4)', npts=1)
    run(htk.banana(2), '(1,0)', 'banana(2) (n=1)', npts=2)
    run(htk.banana(3), '(1,0)', 'banana(3) (n=2)', npts=2)
