"""High-precision evaluation of  c = s * int_0^1 f(xi) dxi,  f = sum_k R_k(xi) log(1-xi)^k,
R_k = N_k(xi)/xi^p_k (exported from Mathematica after Together; no pole at xi = 1).
 [0, 1/2]: exact Taylor series of f at xi = 0 (rational arithmetic; the negative powers must cancel exactly),
           integrated termwise; truncation after M terms, tail bounded by the last terms (geometric, ratio <= 1/2 * growth).
 [1/2, 1]: (i) tanh-sinh (mpmath.quad), (ii) Gauss-Legendre on geometric subintervals [1-2^-j, 1-2^-(j+1)]
           plus the remainder below 2^-J estimated by bounding |f| ~ C log^3."""
import sys, re
from fractions import Fraction as Fr
import mpmath as mp

def load(path):
    lines = open(path).read().strip().split('\n')
    terms = []
    i = 0
    while i < len(lines):
        m = re.match(r'k (\d+) p (\d+)', lines[i])
        k, p = int(m.group(1)), int(m.group(2))
        coeffs = [Fr(c) for c in lines[i + 1].split()] if lines[i + 1].strip() else []
        terms.append((k, p, coeffs)); i += 2
    return terms

def series_part(terms, M):
    """coefficients a_n (n = 0..M-1) of f(xi) = sum_n a_n xi^n; checks cancellation of negative powers."""
    pmax = max(p for k, p, cs in terms)
    L = [Fr(0)] + [Fr(-1, n) for n in range(1, M + pmax + 2)]   # log(1-xi)
    def mul(a, b, n):
        out = [Fr(0)] * n
        for i, x in enumerate(a[:n]):
            if x:
                for j, y in enumerate(b[:n - i]):
                    if y:
                        out[i + j] += x * y
        return out
    n = M + pmax + 1
    Lk = {0: [Fr(1)] + [Fr(0)] * (n - 1)}
    for k in range(1, max(k for k, p, cs in terms) + 1):
        Lk[k] = mul(Lk[k - 1], L, n)
    tot = {}
    for k, p, cs in terms:
        prod = mul(cs + [Fr(0)] * n, Lk[k], n)
        for e, x in enumerate(prod):
            if x:
                tot[e - p] = tot.get(e - p, Fr(0)) + x
    neg = {e: x for e, x in tot.items() if e < 0 and x != 0}
    a = [tot.get(e, Fr(0)) for e in range(M)]
    return a, neg

def run(base, sign=1, M=300, dps=80):
    terms = load(base + '_HLexport.txt')
    a, neg = series_part(terms, M)
    mp.mp.dps = dps + 20
    A = mp.mpf(0)
    for nn, x in enumerate(a):
        A += mp.mpf(x.numerator) / x.denominator * mp.mpf(0.5) ** (nn + 1) / (nn + 1)
    lastterms = [abs(mp.mpf(x.numerator) / x.denominator * mp.mpf(0.5) ** (nn + 1) / (nn + 1)) for nn, x in enumerate(a[-5:], start=M - 5)]
    def f(x):
        Lg = mp.log(1 - x); tot = mp.mpf(0)
        for k, p, cs in terms:
            s = mp.mpf(0)
            for cc in reversed(cs):
                s = s * x + mp.mpf(cc.numerator) / cc.denominator
            tot += s / x ** p * Lg ** k
        return tot
    B1 = mp.quad(f, [mp.mpf(1) / 2, 1], method='tanh-sinh')
    J = 230
    pts = [1 - mp.mpf(2) ** (-j) for j in range(1, J + 1)]
    B2 = mp.mpf(0)
    for lo, hi in zip(pts[:-1], pts[1:]):
        B2 += mp.quad(f, [lo, hi], method='gauss-legendre')
    rem = abs(f(pts[-1])) * (1 - pts[-1]) * 4      # |f| ~ log^3 near 1: crude bound for the last piece
    return sign * (A + B1), sign * (A + B2), neg, max(lastterms), rem

if __name__ == "__main__":
    base = sys.argv[1]; sign = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    v1, v2, neg, tail, rem = run(base, sign)
    mp.mp.dps = 70
    print(base, "negative powers left after cancellation:", neg)
    print("  series[0,1/2] + tanh-sinh[1/2,1]      ", mp.nstr(v1, 60))
    print("  series[0,1/2] + Gauss-Legendre(geom)   ", mp.nstr(v2, 60))
    print("  |difference| %s   last series terms %s   GL remainder bound %s" % (mp.nstr(abs(v1 - v2), 3), mp.nstr(tail, 3), mp.nstr(rem, 3)))
