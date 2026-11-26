"""Exact checks of the closed two-loop formulas for all theta graphs with D = H + T <= 6, T in {0, 1}.

(1) t-space density: brute force (own Kirchhoff + forms, cycle matrix C3) == closed D(t), root N',
    sources at all non-root vertices (N and interior);
(1b) same with the root moved to every interior vertex (graph relabelled by the transposition
    rho <-> N'), and the htk sign rule;
(2) chart density: D(t(xi)) |det[t; dt/dxi]| == closed g(w, x), sign of det[t; dt/dxi] = (-1)^{q-1},
    |det| = prod_P w_P^{-(n_P+1)};
(3) brute force with C3 against htk.brute.density (positive basis): ratio = sgn det[B^T; C3]^D,
    and sgn det[B^T; C3] = (-1)^{q+1};
(4) z = 0 kernel: closed formula == own exact integral of the chart density (T = 0: full quadratic
    moments, covariance terms included) == htk.kernel_exact;
(5) the expansions into products of two H x H determinants (expand_T0, expand_T1) == Q(bar L), D_theta;
(6) displacement map: -sum_e y_e.z_e (brute-force currents) == sum_v lam_v.X_v, root N' and interior root.
Writes exact_results.json.
"""
import json
import random
import sys
import time

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import *  # noqa
import htk
from htk import brute as hbrute
from htk.simplex import PiNumber

rng = random.Random(20261120)


def rf():
    return Fr(rng.randint(-9, 9), rng.randint(1, 4))


def rt():
    return Fr(rng.randint(1, 40), rng.randint(1, 40))


def rsimplex(k):
    """random rational point of the open simplex with k barycentric coordinates"""
    a = [Fr(rng.randint(1, 50)) for _ in range(k)]
    s = sum(a)
    return [x / s for x in a]


def ext_ratio(b, c):
    """b, c: (coeff, radicand) -> True if equal as numbers coeff * radicand^{-1/2}"""
    return (b[0] > 0) == (c[0] > 0) and (b[0] == 0) == (c[0] == 0) and b[0] ** 2 / b[1] == c[0] ** 2 / c[1]


def relabel(nv, edges, rho):
    """transposition rho <-> nv-1"""
    lab = list(range(nv))
    lab[rho], lab[nv - 1] = nv - 1, rho
    return [(lab[a], lab[b]) for a, b in edges], lab


def main():
    out = {"cases": []}
    t0 = time.time()
    npts = 3
    allok = True
    for T in (0, 1):
        for D in range(1, 7):
            H = D - T
            if H < 0 or (T == 0 and H < 1):
                continue
            for pqr in theta_list(D):
                p, q, r = pqr
                nv, edges = theta_graph(*pqr)
                root = nv - 1
                th = Theta(nv, edges)
                rec = {"T": T, "H": H, "pqr": list(pqr)}
                # (1) t-space, root N'
                ok1 = True
                for _ in range(npts):
                    lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                    t = [rt() for _ in edges]
                    y = currents(nv, edges, root, t, lam) if H else None
                    if H:
                        ok1 &= check_kirchhoff(nv, edges, root, t, lam, y, th.C3())
                    b = brute_density(nv, edges, root, t, lam, H, T, th.C3())
                    c = closed_density_t(th, t, lam, H, T)
                    ok1 &= (b == c)
                rec["tspace_rootNp"] = ok1
                # (1b) interior roots
                ok1b, nroots = True, 0
                for P in "ABC":
                    for rho in th.interior[P]:
                        nroots += 1
                        e2, lab = relabel(nv, edges, rho)
                        th2 = Theta(nv, e2)  # trivalent vertices 0 and rho (old N')
                        assert th2.N == 0 and th2.Np == rho and th2.pqr == pqr
                        lam = full_sources(nv, nv - 1, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                        t = [rt() for _ in edges]
                        b = brute_density(nv, e2, nv - 1, t, lam, H, T, th2.C3())
                        c = closed_density_t(th2, t, lam, H, T)
                        ok1b &= (b == c)
                        # htk rule: K[relabelled] = sgn(transposition)^D K[original] for C_htk;
                        # check with the htk brute layer pointwise
                        G1 = htk.Graph(nv, edges)
                        G2 = htk.Graph(nv, e2)
                        lam_old = [None] * nv
                        for v in range(nv):
                            lam_old[v] = lam[lab[v]]
                        d1 = hbrute.density(G1, t, [list(x) for x in lam_old[:-1]] if H else None, T)
                        d2 = hbrute.density(G2, t, [list(x) for x in lam[:-1]] if H else None, T)
                        ok1b &= (d2.same_as(d1.scaled(sgn_power(D))))
                rec["tspace_interior_roots"] = [ok1b, nroots]
                # (2) chart
                ok2 = True
                for _ in range(npts):
                    lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                    w = rsimplex(3)
                    x = {P: rsimplex(th.n[P]) for P in "ABC"}
                    t = chart_t(th, w, x)
                    J = chart_jacobian(th, w, x)
                    mono = Fr(1)
                    for i, P in enumerate("ABC"):
                        mono *= w[i] ** (-(th.n[P] + 1))
                    ok2 &= (abs(J) == mono) and ((J > 0) == (sgn_power(q - 1) > 0))
                    b = brute_density(nv, edges, root, t, lam, H, T, th.C3())
                    gb = (b[0] * abs(J), b[1])
                    gc = closed_chart_density(th, w, x, lam, H, T)
                    if T == 0:
                        ok2 &= gb == gc
                    else:
                        # gb = coeff |J| U^{-1/2}, gc = coeff' (w_A w_B w_C)^{+1/2}
                        ok2 &= ((gb[0] > 0) == (gc[0] > 0)) and gb[0] ** 2 / gb[1] == gc[0] ** 2 * gc[1]
                rec["chart"] = ok2
                # (3) conversion to htk positive basis
                G = htk.Graph(nv, edges)
                s3 = sgn_detBC(nv, edges, root, th.C3())
                ok3 = s3 == sgn_power(q + 1)
                for _ in range(2):
                    lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                    t = [rt() for _ in edges]
                    b = brute_density(nv, edges, root, t, lam, H, T, th.C3())
                    dh = hbrute.density(G, t, [list(x) for x in lam[:-1]] if H else None, T)
                    ok3 &= (dh.coeff == s3 ** D * b[0]) and dh.radicand == b[1]
                rec["htk_basis_sign"] = [ok3, s3]
                # (4) kernels at z = 0
                lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                Kc3 = kernel_C3_closed(th, lam, H, T)
                Kh = kernel_htk_closed(th, lam, H, T)
                if T == 0:
                    # own exact integration of prod w^{H-n} Q(L(sigma)) with the full quadratic moments
                    Lb = Lbar(th, lam)
                    vecs = [(th.N, None)] + [(u, (P, k)) for P in "ABC" for k, u in enumerate(th.interior[P], 1)]
                    mean = {th.N: Fr(1)}
                    for P in "ABC":
                        nP = th.n[P]
                        for k, u in enumerate(th.interior[P], 1):
                            mean[u] = Fr(nP - k, nP)

                    def pol(a, b_):
                        ap = [x + y for x, y in zip(a, b_)]
                        am = [x - y for x, y in zip(a, b_)]
                        return (Q_form(th, lam, ap) - Q_form(th, lam, am)) / 4

                    def second(u, v):
                        if u == th.N or v == th.N:
                            return mean[u] * mean[v]
                        (P1, k1), (P2, k2) = dict(vecs)[u], dict(vecs)[v]
                        if P1 != P2:
                            return mean[u] * mean[v]
                        nP = th.n[P1]
                        return Fr((nP - max(k1, k2)) + (nP - k1) * (nP - k2), nP * (nP + 1))
                    tot = Fr(0)
                    cov = Fr(0)
                    for u, _ in vecs:
                        for v, _ in vecs:
                            qq = pol(lam[u], lam[v])
                            tot += second(u, v) * qq
                            cov += (second(u, v) - mean[u] * mean[v]) * qq
                    wint = Fr(1)
                    for n_ in pqr:
                        wint *= math.factorial(H - n_)
                    wint /= math.factorial(H + 1)
                    for n_ in pqr:
                        wint /= math.factorial(n_ - 1)
                    own = sgn_power(H * (H + 1) // 2 + 1) * wint * tot
                    ok4a = (own == Kc3) and cov == 0 and tot == Q_form(th, lam, Lb)
                else:
                    ok4a = True
                    own = None
                t1 = time.time()
                Ke = htk.kernel_exact(G, [list(x) for x in lam[:-1]] if H else None, T)
                tex = time.time() - t1
                if T == 0:
                    ok4b = (Fr(Ke) == Kh)
                else:
                    ok4b = (Ke == PiNumber({1: Kh}))
                rec["kernel_z0"] = {"closed_equals_own_moments": ok4a, "closed_equals_htk_exact": ok4b,
                                    "K_htk": str(Kh) + ("" if T == 0 else "*pi"), "htk_time_s": round(tex, 2)}
                # (5) expansions
                ok5 = True
                terms = expand_T0(th, H) if T == 0 else expand_T1(th, H)
                for _ in range(2):
                    lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                    ref = Q_form(th, lam, Lbar(th, lam)) if T == 0 else D_theta(th, lam)
                    ok5 &= (eval_terms(terms, th, lam) == ref)
                rec["expansion"] = [ok5, len(terms)]
                # (6) displacement map
                ok6 = True
                if H:
                    for _ in range(npts):
                        lam = full_sources(nv, root, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                        w = rsimplex(3)
                        x = {P: rsimplex(th.n[P]) for P in "ABC"}
                        t = chart_t(th, w, x)
                        y = currents(nv, edges, root, t, lam)
                        z = [[rf() for _ in range(H)] for _ in edges]
                        lhs = -sum(y[e][i] * z[e][i] for e in range(len(edges)) for i in range(H))
                        X = X_map(th, w, x, z)
                        rhs = sum(lam[v][i] * X[v][i] for v in range(nv) for i in range(H))
                        ok6 &= lhs == rhs
                        # interior root rho: X^(rho)_v = X_v - X_rho
                        for P in "ABC":
                            for rho in th.interior[P][:1]:
                                e2, lab = relabel(nv, edges, rho)
                                lam2 = full_sources(nv, nv - 1, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                                th2 = Theta(nv, e2)
                                y2 = currents(nv, e2, nv - 1, t, lam2)
                                lhs2 = -sum(y2[e][i] * z[e][i] for e in range(len(edges)) for i in range(H))
                                X2 = X_map(th2, w, x, z)
                                rhs2 = sum(lam2[v][i] * (X2[v][i] - X2[nv - 1][i]) for v in range(nv) for i in range(H))
                                ok6 &= lhs2 == rhs2
                rec["displacement"] = ok6
                good = ok1 and ok1b and ok2 and ok3 and ok4a and ok4b and ok5 and ok6
                allok &= good
                rec["all"] = good
                out["cases"].append(rec)
                print("T=%d H=%d %-10s tspace %s interior-roots %s chart %s htk-sign %s(%+d) z0 own %s htk %s  K=%s  expansion %s(%d terms) X %s" % (
                    T, H, pqr, ok1, ok1b, ok2, ok3, s3, ok4a, ok4b, rec["kernel_z0"]["K_htk"], ok5, len(terms), ok6), flush=True)
    out["all_pass"] = allok
    out["time_s"] = round(time.time() - t0, 1)
    print("all pass:", allok, " time %.1f s" % (time.time() - t0))
    with open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/exact_results.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
