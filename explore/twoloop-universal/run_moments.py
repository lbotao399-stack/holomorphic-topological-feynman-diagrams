"""z-dependence: moments int g(w,x) prod (y^i_e)^p over Delta^2 x prod_P Delta^{n_P-1} from the integral
representation (closed chart density g, polynomial currents y(w, x)), integrated exactly (Dirichlet),
compared with htk.kernel_exact(..., powers) (engine reduction).  Since -sum y.z = sum lam.X pointwise
(run_exact.py, check 6), these moments are the Taylor coefficients of K(lam; z) = int g e^{sum lam.X}.
Writes moments_results.json."""
import json
import random
import sys
import time

import sympy as sp

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import *  # noqa
import htk
from htk.simplex import PiNumber

rng = random.Random(99)


def rv(H):
    return [Fr(rng.randint(-6, 6), rng.randint(1, 3)) for _ in range(H)]


def g_half(k2):
    """Gamma(k2/2) / sqrt(pi)^{k2 odd} as a Fraction; k2 = 2x, x in (1/2) Z, x > 0."""
    if k2 % 2 == 0:
        return Fr(math.factorial(k2 // 2 - 1)), 0
    m = (k2 - 1) // 2          # x = m + 1/2
    return Fr(math.factorial(2 * m), 4 ** m * math.factorial(m)), 1


def moment(th, lamfull, H, T, powers):
    """exact int_chart g * prod (y^i_e)^p ; returns Fraction (T=0) or coefficient of pi (T=1)."""
    wA, wB, wC = sp.symbols("wA wB wC")
    W = {"A": wA, "B": wB, "C": wC}
    xs = {P: sp.symbols("x%s1:%d" % (P, th.n[P] + 1)) for P in "ABC"}
    lam = {v: [sp.Rational(c.numerator, c.denominator) for c in vec] for v, vec in enumerate(lamfull)}
    Qv = {}
    mhat = {}
    for P in "ABC":
        acc = [sp.Integer(0)] * H
        Qs = []
        for j in range(th.n[P]):
            Qs.append(list(acc))
            if j < th.n[P] - 1:
                u = th.interior[P][j]
                acc = [a + b for a, b in zip(acc, lam[u])]
        Qv[P] = Qs
        mhat[P] = [sum(xs[P][j] * Qs[j][i] for j in range(th.n[P])) for i in range(H)]
    L = [lam[th.N][i] + sum(mhat[P][i] for P in "ABC") for i in range(H)]
    y = {}
    for P in "ABC":
        for j, (e, s) in enumerate(th.paths[P]):
            y[e] = [s * (W[P] * L[i] - mhat[P][i] + Qv[P][j][i]) for i in range(H)]
    F = sp.Integer(1)
    for (e, i), pw in powers.items():
        F *= y[e][i] ** pw
    if T == 0:
        # Q(L) = L^T G L, G from polarization (exact)
        basis = [[Fr(int(a == b)) for b in range(H)] for a in range(H)]
        G = [[(Q_form(th, lamfull, [x + y_ for x, y_ in zip(basis[a], basis[b])]) -
               Q_form(th, lamfull, [x - y_ for x, y_ in zip(basis[a], basis[b])])) / 4 for b in range(H)] for a in range(H)]
        QL = sum(sp.Rational(G[a][b].numerator, G[a][b].denominator) * L[a] * L[b] for a in range(H) for b in range(H))
        pref = th.order_sign * sgn_power(H * (H + 1) // 2 + 1)
        integrand = sp.expand(pref * QL * F)
        aexp2 = {P: 2 * (H - th.n[P]) for P in "ABC"}      # twice the w exponents
    else:
        Dt = D_theta(th, lamfull)
        pref = th.order_sign * sgn_power(H * (H + 1) // 2) * sp.Rational(Dt.numerator, Dt.denominator) / 2
        integrand = sp.expand(pref * F)
        aexp2 = {P: 2 * H + 1 - 2 * th.n[P] for P in "ABC"}
    gens = [wA, wB, wC] + [x for P in "ABC" for x in xs[P]]
    poly = sp.Poly(integrand, *gens)
    total = Fr(0)
    pipow = None
    for mon, coef in poly.terms():
        c = Fr(int(sp.numer(coef)), int(sp.denom(coef)))
        # w part: prod Gamma(c_P + 1) / Gamma(sum c_P + 3)
        num, k = Fr(1), 0
        s2 = 0
        for idx, P in enumerate("ABC"):
            e2 = aexp2[P] + 2 * mon[idx]           # 2 c_P
            gv, kk = g_half(e2 + 2)
            num *= gv
            k += kk
            s2 += e2
        gv, kk = g_half(s2 + 6)
        val = num / gv
        k -= kk
        # x part
        off = 3
        for P in "ABC":
            ex = mon[off:off + th.n[P]]
            off += th.n[P]
            nn = 1
            for b_ in ex:
                nn *= math.factorial(b_)
            val *= Fr(nn, math.factorial(th.n[P] - 1 + sum(ex)))
        assert pipow is None or pipow == k
        pipow = k
        total += c * val
    return total, pipow


def main():
    cases = [(0, 2, (1, 2, 2)), (0, 3, (1, 3, 3)), (0, 3, (2, 2, 3)), (0, 4, (2, 3, 4)), (0, 4, (3, 3, 3)),
             (1, 1, (1, 2, 2)), (1, 2, (1, 3, 3)), (1, 2, (2, 2, 3)), (1, 3, (2, 3, 4)), (1, 3, (3, 3, 3)),
             (1, 5, (2, 5, 6))]
    out = []
    allok = True
    for T, H, pqr in cases:
        nv, edges = theta_graph(*pqr)
        th = Theta(nv, edges)
        G = htk.Graph(nv, edges)
        L = [rv(H) for _ in range(nv - 1)]
        lamfull = full_sources(nv, nv - 1, {v: L[v] for v in range(nv - 1)})
        s = sgn_detBC(nv, edges, nv - 1, th.C3()) ** (H + T)
        Ne = len(edges)
        plist = [{}]
        for _ in range(2):
            plist.append({(rng.randrange(Ne), rng.randrange(H)): 1})
        for _ in range(2):
            e1, e2 = rng.randrange(Ne), rng.randrange(Ne)
            i1, i2 = rng.randrange(H), rng.randrange(H)
            d = {}
            d[(e1, i1)] = d.get((e1, i1), 0) + 1
            d[(e2, i2)] = d.get((e2, i2), 0) + 1
            plist.append(d)
        for powers in plist:
            t0 = time.time()
            mine, k = moment(th, lamfull, H, T, powers)
            mine *= s
            t1 = time.time()
            ref = htk.kernel_exact(G, L, T, powers=powers or None)
            t2 = time.time()
            if T == 0:
                ok = Fr(ref) == mine and k == 0
                if not powers:
                    ok &= mine == kernel_htk_closed(th, lamfull, H, T)
            else:
                ok = ref == PiNumber({1: mine}) and k == 2   # k counts powers of sqrt(pi)
                if not powers:
                    ok &= mine == kernel_htk_closed(th, lamfull, H, T)
            allok &= ok
            rec = {"T": T, "H": H, "pqr": pqr, "powers": {"%d,%d" % kk: v for kk, v in powers.items()},
                   "value": str(mine) + ("*pi" if T else ""), "equal_htk_exact": ok,
                   "t_own_s": round(t1 - t0, 2), "t_htk_s": round(t2 - t1, 2)}
            out.append(rec)
            print(rec, flush=True)
    print("all moments agree:", allok)
    with open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/moments_results.json", "w") as f:
        json.dump({"all_pass": allok, "cases": out}, f, indent=1)


if __name__ == "__main__":
    main()
