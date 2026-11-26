"""Part (c), kernels: theta(1,3,3), theta(2,2,3) at (H, T) = (2, 1).

 (1) engine K(lambda; 0) with symbolic sources at the four interior vertices (lambda_N = 0);
 (2) independent pointwise check (mybrute.py, own Kirchhoff solution, own rho_2 from the Pfaffian
     definition) of the closed-form density in the series-parallel coordinates
        D(t(xi)) |det[t; d_xi t]| = s * D_theta/2 * prod_P w_P^{H + 1/2 - n_P},
     xi = (w_B, w_C, x^A_2.., x^B_2.., x^C_2..), R_P = 1/w_P, t^P_j = R_P x^P_j,
     at random rational points, with random sources at ALL vertices (sum zero), sign included;
 (3) exact integral of the closed form: s * D/2 * prod_P Gamma(H+3/2-n_P) / (Gamma(H+3/2) prod_P (n_P-1)!);
 (4) engine oracle (direct lattice integration of the brute-force density, no reduction).
"""
import sys, random, time
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import theta
import mybrute as mb

H, T = 2, 1
rng = random.Random(7)


def br(a, b):
    return a[0] * b[1] - a[1] * b[0]


def Dtheta_H2(lam, interior):
    """det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]] (S3 eq-Dtheta) for H = 2."""
    A, B, C = [[lam[v] for v in P] for P in interior]
    cols = len(A) + len(B) + len(C)
    assert cols == 2 * H
    M = sp.zeros(4, 4)
    j = 0
    for v in A:
        M[0, j], M[1, j] = v[0], v[1]; j += 1
    for v in B:
        M[0, j], M[1, j] = v[0], v[1]; M[2, j], M[3, j] = v[0], v[1]; j += 1
    for v in C:
        M[2, j], M[3, j] = v[0], v[1]; j += 1
    return M.det()


for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, interior = theta(*pqr)
    nP = list(pqr)
    print("=== theta%s  edges %s  interior %s ===" % (pqr, G.edges, interior))
    # (1) engine, symbolic
    L = {v: list(sp.symbols('l%d_1:3' % v)) for P in interior for v in P}
    L[0] = [0, 0]
    t0 = time.time()
    K = sp.factor(htk.kernel_exact(G, L, T=T))
    print("(1) engine K =", K, "  (%.2fs)" % (time.time() - t0))
    Dsym = sp.expand(Dtheta_H2(L, interior))
    print("    D_theta =", sp.factor(Dsym), ";  K / (pi D_theta) =", sp.simplify(K / (sp.pi * Dsym)))
    # closed form constant
    s = (-1) ** ((H + 1) // 2 + (H + 1) * (pqr[1] + 1))  # ceil(H/2) = (H+1)//2
    cst = sp.Rational(1, 2) * sp.prod([sp.gamma(sp.Rational(2 * H + 3, 2) - k) for k in nP]) / (
        sp.gamma(sp.Rational(2 * H + 3, 2)) * sp.prod([sp.factorial(k - 1) for k in nP]))
    print("(3) closed form: s =", s, ", s * prod Gamma(H+3/2-n_P)/(2 Gamma(H+3/2) prod (n_P-1)!) =", sp.nsimplify(s * cst),
          "; equals K/D_theta:", sp.simplify(s * cst - K / Dsym) == 0)
    # (2) pointwise check with sources at all vertices
    xs = []
    names = ['wB', 'wC']
    for Pn, k in zip('ABC', nP):
        names += ['x%s%d' % (Pn, j) for j in range(2, k + 1)]
    xi = sp.symbols(' '.join(names))
    wB, wC = xi[0], xi[1]
    wA = 1 - wB - wC
    w = {'A': wA, 'B': wB, 'C': wC}
    idx = 2
    tlist = []
    for Pn, k in zip('ABC', nP):
        xr = list(xi[idx: idx + k - 1]); idx += k - 1
        x1 = 1 - sum(xr)
        for xj in [x1] + xr:
            tlist.append(xj / w[Pn])
    J = sp.Matrix([tlist] + [[sp.diff(tt, v) for tt in tlist] for v in xi])
    ok = True
    signs = []
    for trial in range(4):
        # random point in the domain
        while True:
            vals = {}
            a = sorted([Fr(rng.randint(1, 999), 1000) for _ in range(2)])
            vals[wB], vals[wC] = a[0], a[1] - a[0]
            if 1 - vals[wB] - vals[wC] <= 0:
                continue
            idx = 2
            good = True
            for k in nP:
                cuts = sorted(Fr(rng.randint(1, 999), 1000) for _ in range(k - 1))
                parts = [b - c for b, c in zip(cuts + [Fr(1)], [Fr(0)] + cuts)]
                for j in range(k - 1):
                    vals[xi[idx + j]] = parts[j + 1]
                idx += k - 1
                if any(pp <= 0 for pp in parts):
                    good = False
            if good:
                break
        sub = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in vals.items()}
        t = [Fr(str(sp.nsimplify(tt.subs(sub)))) for tt in tlist]
        detJ = Fr(str(J.subs(sub).det()))
        lam = {v: [rng.randint(-6, 6), rng.randint(-6, 6)] for v in range(G.nv - 1)}
        q, U = mb.density(G.nv, G.edges, t, lam, T=1)       # D = q / sqrt(U)
        lhs_q, lhs_rad = q * abs(detJ), U                    # lhs = lhs_q / sqrt(lhs_rad); orientation iota_E vol > 0
        signs.append((detJ > 0) - (detJ < 0))
        lamfull = dict(lam)
        lamfull[G.nv - 1] = [-sum(lam[v][i] for v in lam) for i in range(2)]
        Dn = Dtheta_H2({v: lamfull[v] for v in range(G.nv)}, interior)
        wv = {'A': 1 - vals[wB] - vals[wC], 'B': vals[wB], 'C': vals[wC]}
        # rhs = s D/2 prod w^{H+1/2-n_P} = r_q * prod w^{1/2} ;  compare squares and signs
        rq = Fr(s) * Fr(int(Dn), 2)
        for Pn, k in zip('ABC', nP):
            rq *= wv[Pn] ** (H - k)
        rhs_sq = rq * rq * wv['A'] * wv['B'] * wv['C']
        lhs_sq = lhs_q * lhs_q / lhs_rad
        sign_ok = ((lhs_q > 0) - (lhs_q < 0)) == ((rq > 0) - (rq < 0))
        ok &= (lhs_sq == rhs_sq) and sign_ok and rq != 0
    print("(2) pointwise identity D(t(xi)) |det[t; d_xi t]| = s D_theta/2 prod w^(H+1/2-n_P) (4 random points, sources at all vertices, sign included):", ok, "; sgn det[t; d_xi t] =", signs)
    # (4) oracle
    if '--no-oracle' in sys.argv:
        continue
    lamn = [[rng.randint(-5, 5), rng.randint(-5, 5)] for v in range(G.nv - 1)]
    Kex = htk.kernel_exact(G, lamn, T=T)
    lamfull = {v: lamn[v] for v in range(G.nv - 1)}
    lamfull[G.nv - 1] = [-sum(l[i] for l in lamn) for i in range(2)]
    Dn = Dtheta_H2(lamfull, interior)
    t0 = time.time()
    val, err = htk.kernel_numeric(G, lamn, T=T)
    exact_f = float(sp.N(s * cst * Dn, 20))
    print("(4) sources %s: D_theta = %s, exact %s = %.12f, engine exact %s; oracle %.12f +- %.1e; "
          "rel. dev. %.2e (%.0fs)" % (lamn, Dn, sp.nsimplify(s * cst * Dn), exact_f, Kex, val, err,
                                       abs(val - exact_f) / abs(exact_f), time.time() - t0), flush=True)
