"""Trace identities for S_H in fundamental traces p_j = tr x^j; anomaly-free and Green-Schwarz conditions.

Generating functions (exact, verified on explicit matrices):
  V (x) V*:  Tr e^{sX} = tr e^{sX} tr e^{-sX}      -> gl/sl adjoint
  Lambda^2V: Tr e^{sX} = (tr(e^{sX})^2 - tr e^{2sX})/2  -> so adjoint
  Sym^2 V:   Tr e^{sX} = (tr(e^{sX})^2 + tr e^{2sX})/2  -> sp adjoint
"""
from fractions import Fraction as F
from math import comb
import random
import sympy as sp

from common import Log
from lie import basis_sl, basis_so, basis_sp, adjoint, sym2, alt2, dual, madd, mat, tr_pow

log = Log("casimirs.log")
rng = random.Random(11)
N = sp.symbols("N")
p = sp.symbols("p0:8")   # p[j] = tr x^j, p0 = N


def gen_VVd(k):
    return sum(comb(k, j) * (-1) ** (k - j) * p[j] * p[k - j] for j in range(k + 1))


def gen_L2(k):
    return sp.Rational(1, 2) * (sum(comb(k, j) * p[j] * p[k - j] for j in range(k + 1)) - 2 ** k * p[k])


def gen_S2(k):
    return sp.Rational(1, 2) * (sum(comb(k, j) * p[j] * p[k - j] for j in range(k + 1)) + 2 ** k * p[k])


def rand_elem(B):
    X = mat(len(B[0]))
    for b in B:
        X = madd(X, b, F(rng.randint(-6, 6), rng.randint(1, 4)))
    return X


def subs_p(expr, x, n):
    vals = {p[0]: n}
    for j in range(1, 8):
        vals[p[j]] = tr_pow(x, j)
    return sp.Rational(expr.subs(vals))


ok_all = True
log("(a) explicit-matrix check of the generating-function traces, k = 2..6")
for n in range(2, 6):
    B = basis_sl(n)
    x = rand_elem(B)
    ad = adjoint(B, x)
    for k in range(2, 7):
        ok = (sp.Rational(tr_pow(ad, k)) == subs_p(gen_VVd(k), x, n)
              and sp.Rational(tr_pow(sym2(x), k)) == subs_p(gen_S2(k), x, n)
              and sp.Rational(tr_pow(alt2(x), k)) == subs_p(gen_L2(k), x, n))
        ok_all &= ok
    log("  sl(%d): adj = V x V* - 1, Sym^2, Lambda^2: %s" % (n, ok_all))
for n in range(3, 9):
    B = basis_so(n)
    x = rand_elem(B)
    ad = adjoint(B, x)
    ok = all(sp.Rational(tr_pow(ad, k)) == subs_p(gen_L2(k), x, n) for k in range(2, 7))
    ok_all &= ok
    log("  so(%d): adj = Lambda^2 V: %s" % (n, ok))
for m in range(1, 4):
    B, J = basis_sp(2 * m)
    x = rand_elem(B)
    ad = adjoint(B, x)
    ok = all(sp.Rational(tr_pow(ad, k)) == subs_p(gen_S2(k), x, 2 * m) for k in range(2, 7))
    ok_all &= ok
    log("  sp(%d): adj = Sym^2 V: %s" % (2 * m, ok))


def show(name, expr, odd_zero=False, traceless=False):
    e = expr.subs(p[0], N)
    if odd_zero:
        e = e.subs({p[1]: 0, p[3]: 0, p[5]: 0, p[7]: 0})
    if traceless:
        e = e.subs(p[1], 0)
    log("  %-28s %s" % (name, sp.factor_terms(sp.expand(e))))
    return sp.expand(e)


log("(b) Tr_adj x^k in p_j = tr_fund x^j (sl: p1 = 0; so, sp: p_odd = 0)")
for k in range(2, 7):
    show("sl(N)  Tr_adj x^%d" % k, gen_VVd(k), traceless=True)
for k in (2, 4, 6):
    show("so(N)  Tr_adj x^%d" % k, gen_L2(k), odd_zero=True)
for k in (2, 4, 6):
    show("sp(N)  Tr_adj x^%d" % k, gen_S2(k), odd_zero=True)
for k in (2, 3, 4):
    show("sl(N)  Tr_Lambda2 x^%d" % k, gen_L2(k), traceless=True)
    show("sl(N)  Tr_Sym2 x^%d" % k, gen_S2(k), traceless=True)

log("(c) S_H for classical gauge algebras")
nf, nA, nS, nv, nF = sp.symbols("n_f n_A n_S n_v n_F")
# 2d (0,2): S_1 = -Tr_adj + Tr_R - Tr_F ; quadratic index: tr_R x^2 = T(R) p2
log("  2d (0,2), sl(N), chirals n_f(N+Nb), Fermis n_F N:  S_1 =",
    sp.factor(sp.expand(-gen_VVd(2).subs({p[0]: N, p[1]: 0}) + 2 * nf * p[2] - nF * p[2])))
log("  2d (N+,0) pure: S_1 = -(N+/2) Tr_adj x^2 = -(N+/2) h^vee... (sl(N): -(N+) N p2)")
# 4d: S_2 = Tr_R x^3 (Tr_adj x^3 = 0)
L23 = sp.expand(gen_L2(3).subs({p[0]: N, p[1]: 0}))
S23 = sp.expand(gen_S2(3).subs({p[0]: N, p[1]: 0}))
log("  4d N=1, sl(N): anomaly coefficients A(fund)=1, A(fund*)=-1, A(Lambda2) =", sp.factor(L23 / p[3]),
    ", A(Sym2) =", sp.factor(S23 / p[3]), ", A(adj) =", sp.expand(gen_VVd(3).subs({p[0]: N, p[1]: 0})))
# 6d: S_3 = -Tr_adj + Tr_R  (R = n_f fund + n_A Lambda^2 + n_S Sym^2), sl(N)
adj4 = sp.expand(gen_VVd(4).subs({p[0]: N, p[1]: 0}))
L24 = sp.expand(gen_L2(4).subs({p[0]: N, p[1]: 0, p[3]: 0}))
S24 = sp.expand(gen_S2(4).subs({p[0]: N, p[1]: 0, p[3]: 0}))
S3_sl = sp.expand(-adj4.subs(p[3], 0) + nf * p[4] + nA * L24 + nS * S24)
log("  6d (1,0), sl(N), hypers n_f fund + n_A Lambda2 + n_S Sym2 (p3 terms dropped: absent from")
log("      Tr_adj, Lambda2, Sym2 at k=4 up to p1 p3 = 0):  S_3 =", sp.collect(S3_sl, [p[4], p[2] ** 2]))
log("      irreducible p4 coefficient = 0  <=>", sp.Eq(sp.expand(S3_sl.coeff(p[4])), 0))
so4 = sp.expand(gen_L2(4).subs({p[0]: N, p[1]: 0, p[3]: 0}))
sp4 = sp.expand(gen_S2(4).subs({p[0]: N, p[1]: 0, p[3]: 0}))
log("  6d (1,0), so(N), n_v vector hypers:  S_3 =", sp.collect(sp.expand(-so4 + nv * p[4]), [p[4], p[2] ** 2]))
log("  6d (1,0), sp(N), n_f fund hypers:    S_3 =", sp.collect(sp.expand(-sp4 + nf * p[4]), [p[4], p[2] ** 2]))
# 10d
so6 = sp.expand(gen_L2(6).subs({p[0]: N, p[1]: 0, p[3]: 0, p[5]: 0}))
log("  10d hCS, so(N): S_5 = -Tr_adj x^6 / 2 =", sp.collect(sp.expand(-so6 / 2), [p[6]]),
    " -> N = 32: S_5 =", sp.expand(-so6.subs(N, 32) / 2))
sl6 = sp.expand(gen_VVd(6).subs({p[0]: N, p[1]: 0}))
log("  10d hCS, sl(N): S_5 =", sp.expand(-sl6 / 2))
log("(d) all explicit checks:", ok_all)
log.close()
