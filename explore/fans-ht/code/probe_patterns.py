import sys, pickle, os
sys.path.insert(0, ".")
import sympy as sp
from integrate_general import kernel_general
out = {}
for n in (2, 4, 6, 8):
    res, minexp, nterms = kernel_general(n)
    lam = sp.symbols("l1:%d" % (n + 2))
    expr = 0
    for lk, d in res.items():
        for k, q in d.items():
            if q:
                expr += sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k * sp.Mul(*[x ** pp for x, pp in zip(lam, lk)])
    expr = sp.expand(expr)
    out[n] = expr
    c_end = sp.Rational(4) ** (n // 2) * sp.pi ** (n // 2) / ((n + 1) * sp.factorial(n) ** 2)
    Kn = sp.expand(expr / c_end)
    print("n=%d: K / [(4 pi)^{n/2}/((n+1) n!^2)] has %d monomials" % (n, len(Kn.as_ordered_terms())))
    print("   endpoint coefficient:", Kn.coeff(lam[0], n // 2).coeff(lam[n], n // 2).subs({x: 0 for x in lam}))
    for w in range(2, n + 1):
        sub = {x: 0 for i, x in enumerate(lam) if i not in (0, w - 1, n)}
        print("   sources 1,%d,%d:" % (w, n + 1), sp.factor(Kn.subs(sub)))
pickle.dump({n: str(e) for n, e in out.items()}, open("../general_kernels.pkl", "wb"))
with open("../general_kernels.txt", "w") as f:
    for n, e in out.items():
        f.write("n = %d\nK = %s\n\n" % (n, sp.factor(e)))
