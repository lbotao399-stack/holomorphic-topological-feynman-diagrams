import sys, itertools
sys.path.insert(0, ".")
import sympy as sp
from integrate_general import build, pmul, padd, integrate
n = int(sys.argv[1])
S, L, M = build(n)
lam = sp.symbols("l1:%d" % (n + 2))
def ballot(O):
    c = 0
    for k in range(1, n + 1):
        c += 1 if k in O else -1
        if c < 0:
            return False
    return c == 0
def N(O):
    avail, cnt = 0, 1
    for k in range(1, n + 1):
        if k in O:
            avail += 1
        else:
            cnt *= avail
            avail -= 1
    return cnt
for O in itertools.combinations(range(1, n + 1), n // 2):
    if not ballot(set(O)):
        continue
    C = [k for k in range(1, n + 1) if k not in O]
    p = None
    for i in O:
        p = L[i] if p is None else pmul(p, L[i])
    for j in C:
        p = pmul(p, M[j])
    res, _ = integrate(n, S, p)
    expr = 0
    for lk, d in res.items():
        for k, q in d.items():
            if q:
                expr += sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k * sp.Mul(*[x ** pp for x, pp in zip(lam, lk)])
    print("O=%s C=%s N=%d:" % (O, C, N(set(O))), sp.factor(sp.expand(expr)))
