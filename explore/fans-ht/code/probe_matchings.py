import sys
sys.path.insert(0, ".")
import sympy as sp
from integrate_general import build, pmul, padd, integrate
import fanlib as F
for n in (4, 6):
    S, L, M = build(n)
    lam = sp.symbols("l1:%d" % (n + 2))
    tot = 0
    for m in F.matchings(list(range(1, n + 1))):
        p = None
        for (i, j) in m:
            f = pmul(L[i], M[j])
            p = f if p is None else pmul(p, f)
        res, _ = integrate(n, S, p)
        expr = 0
        for lk, d in res.items():
            for k, q in d.items():
                if q:
                    expr += sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k * sp.Mul(*[x ** pp for x, pp in zip(lam, lk)])
        expr = sp.expand(expr)
        tot += expr
        if n == 4:
            print(m, sp.factor(expr))
    print("n=%d total" % n, sp.factor(sp.expand(tot)) if n == 4 else "")
    # structure: substitute l1 -> x, l2 -> y and check divided-difference form
    K = sp.expand(tot)
    x, y = lam[0], lam[1]
    Kx = sp.Poly(K, x, y)
    print("  degree in (l1,l2):", Kx.total_degree())
