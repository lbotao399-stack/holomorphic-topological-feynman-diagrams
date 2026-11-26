# Check: Tr(A dA dA) restricted to A = sum_a A^a T_a with A^a odd (total degree 1), no dz components,
# d = dz_i d_i + (dzbar, dtau parts).  Top (2,*) part in dz: -dz1 dz2 d_abc A^a d1A^b d2A^c, d_abc = tr(T_a{T_b,T_c}).
# Model: Grassmann algebra on generators; A^a, d1A^a, d2A^a odd symbols; dz1, dz2 odd symbols.
import itertools, random
from fractions import Fraction
# minimal Grassmann algebra: monomials = tuples of generator indices in canonical order, sign from sorting
def mul(x, y):
    out = {}
    for mx, cx in x.items():
        for my, cy in y.items():
            if set(mx) & set(my): continue
            seq = list(mx + my)
            # bubble sort sign
            s = 1
            for i in range(len(seq)):
                for j in range(len(seq) - 1 - i):
                    if seq[j] > seq[j + 1]:
                        seq[j], seq[j + 1] = seq[j + 1], seq[j]; s = -s
            k = tuple(seq); out[k] = out.get(k, 0) + s * cx * cy
    return {k: v for k, v in out.items() if v}
def add(x, y, c=1):
    out = dict(x)
    for k, v in y.items(): out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}
n = 3  # dim of algebra (random structure tensor t_abc = Tr(T_a T_b T_c), not symmetric)
random.seed(1)
t = {(a, b, c): Fraction(random.randint(-5, 5)) for a in range(n) for b in range(n) for c in range(n)}
# generator labels: dz1=0, dz2=1, A^a = 10+a, d1A^a = 20+a, d2A^a = 30+a
g = lambda i: {(i,): 1}
A = [g(10 + a) for a in range(n)]
dA = [add(mul(g(0), g(20 + a)), mul(g(1), g(30 + a))) for a in range(n)]  # dz1 d1A + dz2 d2A (holomorphic part of dA)
lhs = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            lhs = add(lhs, mul(mul(A[a], dA[b]), dA[c]), t[(a, b, c)])
rhs = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            dabc = t[(a, b, c)] + t[(a, c, b)]
            rhs = add(rhs, mul(mul(mul(g(0), g(1)), A[a]), mul(g(20 + b), g(30 + c))), -dabc)
print("Tr(A dA dA)|_(dz1dz2) == -dz1dz2 d_abc A^a d1A^b d2A^c :", add(lhs, rhs, -1) == {})
# F_D integrand d_abc d1A^a d2A^b A^c vs d_abc A^a d1A^b d2A^c (d symmetric)
x = {}; y = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            dabc = t[(a, b, c)] + t[(a, c, b)]
            dsym = Fraction(1, 3) * sum(t[p] + t[(p[0], p[2], p[1])] for p in [(a, b, c), (b, c, a), (c, a, b)]) / 2
            x = add(x, mul(mul(g(20 + a), g(30 + b)), g(10 + c)), dsym)
            y = add(y, mul(mul(g(10 + a), g(20 + b)), g(30 + c)), dsym)
print("d_abc d1A^a d2A^b A^c == d_abc A^a d1A^b d2A^c (d symmetric):", add(x, y, -1) == {})
