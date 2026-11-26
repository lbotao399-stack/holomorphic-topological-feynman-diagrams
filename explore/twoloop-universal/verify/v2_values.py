"""z = 0 kernel values (C5, C6, C7 sign of K_htk, C8 interior root) by own integration of the pointwise-verified
density, compared with the claimed tables/examples and with htk.kernel_exact (oracle).

Own integration chain (Cheng-Wu gauge R_C = 1, t_{P_j} = R_P x^P_j, d^{n}t_P = R_P^{n-1} dR_P d^{n-1}x):
 T=1: K[C3] = s1 (D_theta/2) prod_P 1/(n_P-1)! * J1,  J1 = int_{(0,inf)^2} R_A^{p-1} R_B^{q-1} (R_A R_B + R_A + R_B)^{-(H+3/2)}
 T=0: K[C3] = s0 prod_P 1/(n_P-1)! * J0 * E_x[Q(L(x))],  J0 = int R_A^p R_B^q (R_A R_B + R_A + R_B)^{-(H+2)}
 with alpha = 1/R:  int alpha^{a-1} beta^{b-1} (1+alpha+beta)^{-(a+b+c)} = Gamma(a)Gamma(b)Gamma(c)/Gamma(a+b+c);
 J1, J0 also integrated numerically with mpmath (independent of that identity) for D <= 4.
 E_x[Q(L)] from exact Dirichlet moments of x (own), with the full second moments (no lemma used).
"""
import sys, time, random, math
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.dont_write_bytecode = True
from vcore import *
import htk
import mpmath as mp

rng = random.Random(777)
t0 = time.time()


def bilinear(th, lam, H, u, v):
    z = [0] * H
    cols = [col(u, z, H), col(z, v, H)] + [col(lam[w], z, H) for w in th.interior["A"]] + \
           [col(lam[w], lam[w], H) for w in th.interior["B"]] + [col(z, lam[w], H) for w in th.interior["C"]]
    return det([[c[i] for c in cols] for i in range(2 * H)])


def EQ_full(th, lam, H):
    """E_x[Q(L(x))] with exact first and second moments of all sigma's (no lemma)."""
    terms = [(("N",), lam[th.N])]
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            terms.append(((P, k), lam[u]))
    mom = {P: sigma_moments(th.n[P]) for P in "ABC"}

    def E1(a):
        return Fr(1) if a == ("N",) else mom[a[0]][0][a[1] - 1]

    def E2(a, b):
        if a == ("N",) or b == ("N",):
            return E1(a) * E1(b)
        if a[0] != b[0]:
            return E1(a) * E1(b)
        return mom[a[0]][1][a[1] - 1][b[1] - 1]
    tot = Fr(0)
    for a, la in terms:
        for b, lb in terms:
            tot += E2(a, b) * bilinear(th, lam, H, la, lb)
    Lbar = [sum(E1(a) * Fr(la[i]) for a, la in terms) for i in range(H)]
    return tot, Lbar


def mu_claim(th):
    mu = {th.N: Fr(1), th.Np: Fr(0)}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            mu[u] = Fr(th.n[P] - k, th.n[P])
    return mu


claimed_c1 = {(0, (1, 1, 1)): Fr(1), (1, (1, 2, 2)): Fr(1, 3), (2, (1, 3, 3)): Fr(1, 20), (2, (2, 2, 3)): Fr(1, 30),
              (3, (1, 4, 4)): Fr(1, 252), (3, (2, 3, 4)): Fr(1, 420), (3, (3, 3, 3)): Fr(1, 840),
              (4, (1, 5, 5)): Fr(1, 5184), (4, (2, 4, 5)): Fr(1, 9072), (4, (3, 3, 5)): Fr(1, 10080),
              (4, (3, 4, 4)): Fr(1, 22680), (5, (1, 6, 6)): Fr(1, 158400), (5, (2, 5, 6)): Fr(1, 285120),
              (5, (3, 4, 6)): Fr(1, 332640), (5, (3, 5, 5)): Fr(1, 798336), (5, (4, 4, 5)): Fr(1, 997920)}
claimed_sign1 = {0: [1], 1: [-1], 2: [-1, 1], 3: [1, 1, 1], 4: [1, -1, 1, -1], 5: [-1] * 5}
claimed_c0 = {1: [Fr(1, 2)], 2: [Fr(1, 6)], 3: [Fr(1, 48)] * 2, 4: [Fr(1, 720), Fr(1, 720), Fr(1, 960)],
              5: [Fr(1, 17280)] * 3 + [Fr(1, 25920)], 6: [Fr(1, 604800)] * 3 + [Fr(1, 967680), Fr(1, 1088640)]}
claimed_sign0 = {1: [1], 2: [1], 3: [-1, 1], 4: [-1] * 3, 5: [1, -1, 1, -1], 6: [1] * 5}

mp.mp.dps = 30


def J_numeric(a, b, c):
    """int_{(0,inf)^2} alpha^{a-1} beta^{b-1} (1+alpha+beta)^{-(a+b+c)} numerically (alpha = e^u, beta = e^v)."""
    f = lambda u, v: mp.e ** (a * u + b * v) * (1 + mp.e ** u + mp.e ** v) ** (-(a + b + c))
    return mp.quad(f, [-mp.inf, -5, 0, 5, mp.inf], [-mp.inf, -5, 0, 5, mp.inf])


ok_all = True
for T in (0, 1):
    for D in range(1, 7):
        H = D - T
        if T == 0 and H < 1:
            continue
        for idx, pqr in enumerate(theta_list(D)):
            p, q, r = pqr
            th = theta_std(*pqr)
            nv = th.nv
            fac = Fr(1)
            for n_ in pqr:
                fac /= math.factorial(n_ - 1)
            res = {}
            if T == 1:
                # J1 / pi : Gamma(H+3/2-n) = sqrt(pi) ghalf(H+1-n)
                J1_over_pi = ghalf(H + 1 - p) * ghalf(H + 1 - q) * ghalf(H + 1 - r) / ghalf(H + 1)
                c1 = fac * J1_over_pi / 2
                res["c1 = claim"] = (c1 == claimed_c1[(H, pqr)])
                if D <= 4:
                    a, b, c = [mp.mpf(2 * H + 3 - 2 * n_) / 2 for n_ in pqr]
                    Jn = J_numeric(a, b, c)
                    res["J1 numeric digits"] = round(float(-mp.log10(abs(Jn / (J1_over_pi * mp.pi) - 1) + mp.mpf(10) ** -40)), 1)
                s1 = sg(H * (H + 1) // 2)
                shtk = s1 * sg((H + 1) * (q + 1))
                res["sign_htk = claim"] = (shtk == claimed_sign1[H][idx])
                okh = True
                for _ in range(2):
                    lam = rand_sources(nv, th.Np, H, rng)
                    K3 = s1 * c1 * D_theta(th, lam, H)        # coefficient of pi
                    G, _ = htk.graph.theta(*pqr)
                    Kh = htk.kernel_exact(G, [lam[v] for v in range(nv - 1)], T=1) if H else htk.kernel_exact(G, None, T=1)
                    okh &= (Kh.coeff(1) == sg((q + 1) * D) * K3 and Kh.coeff(0) == 0)
                res["K_htk = (-1)^((q+1)D) K[C3] (htk.kernel_exact)"] = okh
            else:
                J0 = Fr(math.factorial(H - p) * math.factorial(H - q) * math.factorial(H - r), math.factorial(H + 1))
                c0 = fac * J0
                res["c0 = claim"] = (c0 == claimed_c0[H][idx])
                if D <= 4:
                    a, b, c = [mp.mpf(H + 1 - n_) for n_ in pqr]
                    Jn = J_numeric(a, b, c)
                    res["J0 numeric digits"] = round(float(-mp.log10(abs(Jn / mp.mpf(J0.numerator) * J0.denominator - 1) + mp.mpf(10) ** -40)), 1)
                s0 = sg(H * (H + 1) // 2 + 1)
                shtk = s0 * sg(H * (q + 1))
                res["sign_htk = claim"] = (shtk == claimed_sign0[H][idx])
                okh, okl, okmu = True, True, True
                for _ in range(2):
                    lam = rand_sources(nv, th.Np, H, rng)
                    EQ, Lbar = EQ_full(th, lam, H)
                    okl &= (EQ == Qform(th, lam, H, Lbar))       # covariances drop out
                    mu = mu_claim(th)
                    Lmu = [sum(mu[v] * Fr(lam[v][i]) for v in mu) for i in range(H)]
                    okmu &= (Lmu == Lbar)
                    K3 = s0 * c0 * EQ
                    G, _ = htk.graph.theta(*pqr)
                    Kh = htk.kernel_exact(G, [lam[v] for v in range(nv - 1)], T=0)
                    okh &= (Kh == sg((q + 1) * D) * K3)
                res["E[Q(L)] = Q(E[L])"] = okl
                res["Lbar = sum mu_v lam_v (claimed mu)"] = okmu
                res["K_htk = (-1)^((q+1)D) K[C3] (htk.kernel_exact)"] = okh
            ok_all &= all(v for v in res.values() if isinstance(v, bool))
            print("T=%d H=%d %s " % (T, H, pqr) + "  ".join("%s: %s" % kv for kv in res.items()), flush=True)

# ---------------- sigma second moments against sympy iterated integration (n <= 4)
import sympy as sp
okm = True
for n in (2, 3, 4):
    xs = sp.symbols("x1:%d" % (n + 1))
    # uniform on Delta^{n-1}: integrate over x_2..x_n with x_1 = 1 - sum
    def integ(expr):
        e = expr.subs(xs[0], 1 - sum(xs[1:]))
        lims = []
        for j in range(n - 1, 0, -1):
            e = sp.integrate(e, (xs[j], 0, 1 - sum(xs[1:j])))
        return e * math.factorial(n - 1)
    m1, m2 = sigma_moments(n)
    sig = [sum(xs[k:]) for k in range(1, n)]
    for k in range(n - 1):
        okm &= (sp.Rational(m1[k].numerator, m1[k].denominator) == integ(sig[k]))
        for l in range(n - 1):
            okm &= (sp.Rational(m2[k][l].numerator, m2[k][l].denominator) == sp.nsimplify(integ(sig[k] * sig[l])))
print("sigma moments vs sympy iterated integrals (n=2,3,4):", okm)

# ---------------- examples of C5/C6 (htk convention) against own closed values
def names(th):
    nm = {th.N: "N", th.Np: "Np"}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            nm[u] = P + str(k)
    return {v: k for k, v in nm.items()}


def br(lam, vs):
    return det([[Fr(lam[v][i]) for v in vs] for i in range(len(vs))])


oke = True
for _ in range(3):
    # T=1 H=1 (1,2,2): K_htk = -(pi/3) B1 C1
    th = theta_std(1, 2, 2); v = names(th); lam = rand_sources(th.nv, th.Np, 1, rng)
    K = sg(1) * sg(2 * 3) * Fr(1, 3) * D_theta(th, lam, 1)
    oke &= (K == -Fr(1, 3) * lam[v["B1"]][0] * lam[v["C1"]][0])
    # T=1 H=2 (1,3,3): -(pi/20)[B1B2][C1C2]
    th = theta_std(1, 3, 3); v = names(th); lam = rand_sources(th.nv, th.Np, 2, rng)
    K = sg(3) * sg(3 * 4) * Fr(1, 20) * D_theta(th, lam, 2)
    oke &= (K == -Fr(1, 20) * br(lam, [v["B1"], v["B2"]]) * br(lam, [v["C1"], v["C2"]]))
    # T=1 H=2 (2,2,3): +(pi/30)[A1B1][C1C2]
    th = theta_std(2, 2, 3); v = names(th); lam = rand_sources(th.nv, th.Np, 2, rng)
    K = sg(3) * sg(3 * 3) * Fr(1, 30) * D_theta(th, lam, 2)
    oke &= (K == Fr(1, 30) * br(lam, [v["A1"], v["B1"]]) * br(lam, [v["C1"], v["C2"]]))
    # T=1 H=3 (3,3,3): (pi/840)([A1A2B1][B2C1C2] - [A1A2B2][B1C1C2])
    th = theta_std(3, 3, 3); v = names(th); lam = rand_sources(th.nv, th.Np, 3, rng)
    K = sg(6) * sg(4 * 4) * Fr(1, 840) * D_theta(th, lam, 3)
    oke &= (K == Fr(1, 840) * (br(lam, [v["A1"], v["A2"], v["B1"]]) * br(lam, [v["B2"], v["C1"], v["C2"]])
                               - br(lam, [v["A1"], v["A2"], v["B2"]]) * br(lam, [v["B1"], v["C1"], v["C2"]])))
    # T=0 H=2 (1,2,2): -(1/6)[Lbar B1][Lbar C1], Lbar = N + B1/2 + C1/2
    th = theta_std(1, 2, 2); v = names(th); lam = rand_sources(th.nv, th.Np, 2, rng)
    EQ, Lb = EQ_full(th, lam, 2)
    K = sg(4) * sg(2 * 3) * Fr(1, 6) * EQ
    Lc = [lam[v["N"]][i] + lam[v["B1"]][i] / 2 + lam[v["C1"]][i] / 2 for i in range(2)]
    d2 = lambda a, b: a[0] * b[1] - a[1] * b[0]
    oke &= (K == -Fr(1, 6) * d2(Lc, lam[v["B1"]]) * d2(Lc, lam[v["C1"]]))
    # T=0 H=3 (1,3,3): -(1/48)[Lbar B1 B2][Lbar C1 C2]
    th = theta_std(1, 3, 3); v = names(th); lam = rand_sources(th.nv, th.Np, 3, rng)
    EQ, Lb = EQ_full(th, lam, 3)
    K = sg(7) * sg(3 * 4) * Fr(1, 48) * EQ
    Lc = [lam[v["N"]][i] + Fr(2, 3) * lam[v["B1"]][i] + Fr(1, 3) * lam[v["B2"]][i] + Fr(2, 3) * lam[v["C1"]][i]
          + Fr(1, 3) * lam[v["C2"]][i] for i in range(3)]
    d3 = lambda a, b, c: det([[a[i], b[i], c[i]] for i in range(3)])
    oke &= (K == -Fr(1, 48) * d3(Lc, lam[v["B1"]], lam[v["B2"]]) * d3(Lc, lam[v["C1"]], lam[v["C2"]]))
    # T=0 H=3 (2,2,3): +(1/48)[Lbar A1 B1][Lbar C1 C2]
    th = theta_std(2, 2, 3); v = names(th); lam = rand_sources(th.nv, th.Np, 3, rng)
    EQ, Lb = EQ_full(th, lam, 3)
    K = sg(7) * sg(3 * 3) * Fr(1, 48) * EQ
    Lc = [lam[v["N"]][i] + lam[v["A1"]][i] / 2 + lam[v["B1"]][i] / 2 + Fr(2, 3) * lam[v["C1"]][i]
          + Fr(1, 3) * lam[v["C2"]][i] for i in range(3)]
    oke &= (K == Fr(1, 48) * d3(Lc, lam[v["A1"]], lam[v["B1"]]) * d3(Lc, lam[v["C1"]], lam[v["C2"]]))
print("C5/C6 examples (htk convention) equal own closed values:", oke)

# ---------------- C8: interior root, htk relabelled graphs, Lbar with mu_v - mu_rho
ok8, n8 = True, 0
for T, H, pqr in [(0, 2, (1, 2, 2)), (0, 3, (1, 3, 3)), (0, 3, (2, 2, 3)), (0, 4, (2, 3, 4)), (1, 1, (1, 2, 2)),
                  (1, 2, (1, 3, 3)), (1, 2, (2, 2, 3)), (1, 3, (3, 3, 3)), (1, 3, (2, 3, 4))]:
    th = theta_std(*pqr)
    D = H + T
    q = pqr[1]
    for P in "ABC":
        for rho in th.interior[P]:
            perm = list(range(th.nv)); perm[rho], perm[th.Np] = th.Np, rho
            th3 = transform(th, perm)                      # rho now carries the last label
            lam = rand_sources(th.nv, rho, H, rng)          # full sources, lam_rho = -sum, lam_N' free
            lam3 = [None] * th.nv
            for v in range(th.nv):
                lam3[perm[v]] = lam[v]
            G3 = htk.Graph(th.nv, th3.edges)
            Kh3 = htk.kernel_exact(G3, [lam3[v] for v in range(th.nv - 1)], T=T)
            fac = Fr(1)
            for n_ in pqr:
                fac /= math.factorial(n_ - 1)
            if T == 1:
                c1 = fac * ghalf(H + 1 - pqr[0]) * ghalf(H + 1 - pqr[1]) * ghalf(H + 1 - pqr[2]) / ghalf(H + 1) / 2
                K3 = sg(H * (H + 1) // 2) * c1 * D_theta(th, lam, H)
                ok8 &= (Kh3.coeff(1) == sg(D) * sg((q + 1) * D) * K3)
            else:
                c0 = fac * Fr(math.factorial(H - pqr[0]) * math.factorial(H - pqr[1]) * math.factorial(H - pqr[2]),
                              math.factorial(H + 1))
                mu = mu_claim(th)
                Lr = [sum((mu[v] - mu[rho]) * Fr(lam[v][i]) for v in range(th.nv) if v != rho) for i in range(H)]
                K3 = sg(H * (H + 1) // 2 + 1) * c0 * Qform(th, lam, H, Lr)
                EQ, _ = EQ_full(th, lam, H)
                ok8 &= (Kh3 == sg(D) * sg((q + 1) * D) * K3) and (EQ == Qform(th, lam, H, Lr))
            n8 += 1
print("C8 interior roots (%d cases): K_htk[relabelled] = (-1)^D K_htk, Lbar = sum (mu_v - mu_rho) lam_v: %s" % (n8, ok8))
print("ALL VALUE CHECKS PASS:", ok_all and okm and oke and ok8, " time %.1f s" % (time.time() - t0))
