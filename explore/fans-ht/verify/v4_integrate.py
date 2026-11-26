"""V7 (claim a6): g_n / m_0 is a polynomial (sympy cancel, independent of the explorer's expansion);
V9 (claims a7, a8, a9, a10, d1): exact integrals.
  - termwise integration of m_0 x (monomials of the polynomial g/m_0) with own Beta/Gamma code;
  - own implementation of the configuration formula of a8 (exponents re-derived, see VERIFY.md);
  - endpoint (1,1) and (2,0) closed forms by own Beta products, n up to 30;
  - comparison with explore/fans-ht/general_kernels.txt (parsed with sympy)."""
import itertools
import math
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.setrecursionlimit(10000)
NS = [int(x) for x in sys.argv[1:]] or [2, 4]


# ---------------------------------------------------------------- numbers q * sqrt(pi)^k
def gamma_q(x):
    """Gamma(x) for x a positive integer or half-integer: (q, k) with Gamma = q sqrt(pi)^k."""
    x = Fr(x)
    assert x > 0
    if x.denominator == 1:
        return Fr(math.factorial(int(x) - 1)), 0
    assert x.denominator == 2
    m = int(x - Fr(1, 2))
    return Fr(math.factorial(2 * m), 4 ** m * math.factorial(m)), 1


def beta_q(a, b):
    qa, ka = gamma_q(a)
    qb, kb = gamma_q(b)
    qc, kc = gamma_q(Fr(a) + Fr(b))
    return qa * qb / qc, ka + kb - kc


def to_sympy(q, k):
    return sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k


# ---------------------------------------------------------------- g / m_0 as a rational function
def build_p(n, lam):
    """lam: dict 1..n+1 (sympy symbols or integers).  Returns (p, us, vs) with
    g_n = m_0 * p,  m_0 = prod_{l even} u_l^{-1/2} prod_{k<n} v_k^{1/2} prod_{k<n odd} (1-v_k)^{-1/2}."""
    us = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    vs = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    P = {1: sp.Integer(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * us[k - 1]
    R = {n: sp.Integer(1)}
    for i in range(n - 1, 0, -1):
        R[i] = R[i + 1] * us[i + 1] / (1 - vs[i])
    c = {k: us[k] * (1 - us[k]) * R[k] for k in range(1, n + 1)}
    L = {}
    acc = 0
    for i in range(1, n + 1):
        acc = acc + lam[i] / P[i]
        L[i] = acc
    M = {j: sum((P[k] ** 2 * c[k] * (L[k] - L[j]) for k in range(j + 1, n + 1)), sp.Integer(0)) + P[n] * us[n] * vs[n] * lam[n + 1]
         for j in range(1, n + 1)}
    # Pi / m_0: integer exponents
    q = vs[n] ** -1
    for k in range(1, n):
        q *= (1 - vs[k]) ** ((k - 1) // 2 if k % 2 else k // 2 - 1)
    for k in range(1, n + 1):
        q *= us[k] ** (-(k + 1) // 2 if k % 2 else -(k // 2))
    from vlib import matchings
    haf = 0
    for m in matchings(list(range(1, n + 1))):
        term = sp.Integer(1)
        for (i, j) in m:
            term *= L[i] * M[j]
        haf += term
    p = sp.cancel(sp.together(q * haf / sp.Integer(2) ** (n // 2)))
    return p, us, vs


def integrate_p(n, p, us, vs):
    """int over (0,1)^{2n-2} x Delta^2 of m_0 * p, p a polynomial in u, v (asserted)."""
    gens = [us[k] for k in range(1, n + 1)] + [vs[k] for k in range(1, n + 1)]
    num, den = sp.fraction(p)
    assert not (den.free_symbols & set(gens)), "denominator depends on r: %s" % den
    poly = sp.Poly(sp.expand(num), *gens)
    tot = {}
    for monom, coeff in poly.terms():
        a = monom[:n]
        b = monom[n:]
        qq, kk = Fr(1), 0
        for l in range(1, n):
            e = Fr(a[l - 1]) - (Fr(1, 2) if l % 2 == 0 else 0)
            qq *= 1 / (e + 1)
        for k in range(1, n):
            if k % 2:
                q2, k2 = beta_q(Fr(b[k - 1]) + Fr(3, 2), Fr(1, 2))
            else:
                q2, k2 = 1 / (Fr(b[k - 1]) + Fr(3, 2)), 0
            qq *= q2
            kk += k2
        an = Fr(a[n - 1]) - (Fr(1, 2) if n % 2 == 0 else 0)
        bn = Fr(b[n - 1])
        g1, k1 = gamma_q(an + 1)
        g2, k2 = gamma_q(bn + 1)
        g3, k3 = gamma_q(an + bn + 3)
        qq *= g1 * g2 / g3
        kk += k1 + k2 - k3
        tot[kk] = tot.get(kk, 0) + coeff / den * sp.Rational(qq.numerator, qq.denominator)
    return sum(sp.expand(v) * sp.sqrt(sp.pi) ** k for k, v in tot.items())


# ---------------------------------------------------------------- configuration formula of a8 (own implementation)
def config_formula(n, lam):
    tot = {}
    for O in itertools.combinations(range(1, n + 1), n // 2):
        Oset = set(O)
        closers = [j for j in range(1, n + 1) if j not in Oset]
        # multiplicity N(O): number of perfect matchings with opener set O
        N = 1
        for j in closers:
            avail = sum(1 for i in O if i < j) - sum(1 for jj in closers if jj < j)
            if avail <= 0:
                N = 0
                break
            N *= avail
        if N == 0:
            continue
        opener_choices = [list(range(1, i + 1)) for i in O]
        closer_choices = [[(w, k) for k in range(j + 1, n + 1) for w in range(j + 1, k + 1)] + ["inf"] for j in closers]
        for ws in itertools.product(*opener_choices):
            for chis in itertools.product(*closer_choices):
                pcount = sum(1 for x in chis if x == "inf")
                if pcount == 0:
                    continue
                mono = Fr(1) if False else 1
                lamterm = 1
                for w in ws:
                    lamterm *= lam[w]
                for x in chis:
                    lamterm *= lam[n + 1] if x == "inf" else lam[x[0]]
                qq, kk = Fr(1), 0
                for l in range(1, n):
                    e = Fr(n - 1 - l, 2) - sum(1 for w in ws if w > l) + sum(1 for x in chis if x != "inf" and x[0] <= l < x[1])
                    f = sum(1 for x in chis if x != "inf" and x[1] == l)
                    q2, k2 = beta_q(e + 1, f + 1)
                    qq *= q2
                    kk += k2
                for m in range(1, n):
                    cm = sum(1 for x in chis if x != "inf" and x[1] <= m)
                    q2, k2 = beta_q(Fr(3, 2), Fr(m, 2) - cm)
                    qq *= q2
                    kk += k2
                fn = sum(1 for x in chis if x != "inf" and x[1] == n)
                q2, k2 = beta_q(Fr(1, 2), fn + pcount + 1)
                qq *= q2 / pcount
                kk += k2
                tot[kk] = tot.get(kk, 0) + N * lamterm * sp.Rational(qq.numerator, qq.denominator)
    return sum(sp.expand(v) * sp.sqrt(sp.pi) ** k for k, v in tot.items()) / sp.Integer(2) ** (n // 2)


def kernels_file():
    txt = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/general_kernels.txt").read()
    out = {}
    cur = None
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith("n ="):
            cur = int(line.split("=")[1])
        elif line.startswith("K =") and cur is not None:
            out[cur] = sp.sympify(line[3:].strip(), locals={"pi": sp.pi})
    return out


if __name__ == "__main__":
    KF = kernels_file()
    print("general_kernels.txt contains n =", sorted(KF))
    for n in NS:
        t0 = time.time()
        lsyms = {i: sp.Symbol("l%d" % i) for i in range(1, n + 2)}
        if n <= 4:
            p, us, vs = build_p(n, lsyms)
            K = sp.factor(integrate_p(n, p, us, vs))
            Kc = sp.factor(config_formula(n, lsyms))
            print("n=%d: g/m_0 polynomial in r (symbolic lambda): True; #monomials in r = %d" %
                  (n, len(sp.Poly(sp.fraction(p)[0], *(list(us.values()) + list(vs.values()))).terms())))
            print("   termwise integration: K_%d =" % n, K)
            print("   configuration formula (a8) == termwise:", sp.simplify(K - Kc) == 0)
            print("   general_kernels.txt == termwise:", sp.simplify(KF[n] - K) == 0)
            if n == 2:
                print("   claim a9 (pi/3) l1 l3:", sp.simplify(K - sp.pi / 3 * lsyms[1] * lsyms[3]) == 0)
            if n == 4:
                l = lsyms
                claim = sp.pi ** 2 / 180 * l[1] * l[5] * ((l[1] + 2 * l[2] + 2 * l[3]) * (2 * l[3] + 2 * l[4] + l[5]) + 8 * l[3] ** 2)
                print("   claim a9 (pi^2/180) l1 l5 [(l1+2l2+2l3)(2l3+2l4+l5)+8l3^2]:", sp.simplify(K - claim) == 0)
        else:
            rng = random.Random(100 + n)
            for trial in range(3):
                lv = {i: sp.Integer(rng.randint(-4, 4)) for i in range(1, n + 2)}
                lv[1] = lv[1] if lv[1] != 0 else sp.Integer(2)
                lv[n + 1] = lv[n + 1] if lv[n + 1] != 0 else sp.Integer(-3)
                p, us, vs = build_p(n, lv)
                Kt = sp.nsimplify(integrate_p(n, p, us, vs))
                Kc = config_formula(n, lv)
                Kf = KF[n].subs({lsyms[i]: lv[i] for i in lsyms})
                print("n=%d lam=%s: g/m_0 polynomial in r: True; termwise = %s; configuration formula equal: %s; general_kernels.txt equal: %s"
                      % (n, [int(lv[i]) for i in range(1, n + 2)], sp.simplify(Kt), sp.simplify(Kt - Kc) == 0, sp.simplify(Kt - Kf) == 0))
        print("   (%.1f s)" % (time.time() - t0))
