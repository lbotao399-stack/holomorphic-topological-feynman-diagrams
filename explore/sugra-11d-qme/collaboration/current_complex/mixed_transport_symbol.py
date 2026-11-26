"""Coefficient of degree one in the marked-path reduced Wick tensor at B=tM.

W=H^T is diagonal in a generic basis, with product-one eigenvalues. Exponents
of the first four independent eigenvalues may be negative. A power of 256 per
variable encodes each Laurent monomial: expressions below have fewer than
eight eigenvalue factors and ordinary variable degree below eight.
"""

from itertools import permutations
from collections import Counter


def mono(j):
    return {j: 1}


def add(*terms):
    r = Counter()
    for a in terms:
        r.update(a)
    return {k: n for k, n in r.items() if n}


def neg(a):
    return {k: -v for k, v in a.items()}


def prod2(a, b):
    r = {}
    for k, v in a.items():
        for j, w in b.items():
            r[k + j] = r.get(k + j, 0) + v * w
    return {k: v for k, v in r.items() if v}


def prod(*a):
    v = {0: 1}
    for b in a:
        if not b:
            return {}
        v = prod2(v, b)
    return v


def vadd(*a):
    return [add(*r) for r in zip(*a)]


def vneg(a):
    return list(map(neg, a))


def varvec(start):
    return [mono(256 ** (start + i)) for i in range(5)]


w = [256**i for i in range(4)]
w.append(-sum(w))


def W(x, d=1):
    return [{k + d * a: v for k, v in p.items()} for a, p in zip(w, x)]


C, M, eta = varvec(4), varvec(9), varvec(14)
K = [[mono(256 ** (19 + 5 * i + j)) for j in range(5)] for i in range(5)]
K[4][4] = neg(add(*(K[i][i] for i in range(4))))


def act(x):
    return [add(*(prod2(a, b) for a, b in zip(r, x))) for r in K]


F = W(C)
P = vadd(C, vneg(F))
Q = vneg(P)
p = vadd(P, vneg(act(M)))


def lv(x):
    return vadd(W(x), W(x, -1), vneg(x), vneg(x))


def nv(x):
    return vadd(W(x), vneg(x))


def rv(x):
    return vadd(W(x, -1), vneg(x))


PER = [
    (r, (-1) ** sum(a > b for i, a in enumerate(r) for b in r[i + 1 :]))
    for r in permutations(range(5))
]


def eps(*vs):
    c = Counter()
    for perm, sg in PER:
        y = prod(*(vs[i][j] for i, j in enumerate(perm)))
        c.update(y if sg == 1 else neg(y))
    return {k: v for k, v in c.items() if v}


def dh(a, b, c, op=lambda v: v, left=False):
    x, y = (op(eta), p) if left else (eta, op(p))
    z, t = (op(p), eta) if left else (p, op(eta))
    return add(eps(x, y, a, b, c), neg(eps(z, t, a, b, c)))


def V(x, y):
    return eps(x, y, eta, M, p)


if __name__ == "__main__":
    S = add(
        dh(W(P, -1), M, C),
        dh(W(M, -1), M, C, act),
        dh(P, M, C, lambda x: W(x, -1)),
        dh(W(M), P, F),
        dh(W(M), M, F, act),
        dh(M, P, F, W, left=True),
    )
    trace = add(*(mono(t) for t in w), *(neg(mono(-t)) for t in w))
    B = add(
        V(C, lv(P)),
        V(rv(P), Q),
        prod2(trace, V(C, F)),
        neg(V(C, act(lv(M)))),
        V(lv(M), act(F)),
        neg(V(Q, act(nv(M)))),
        neg(V(rv(M), act(Q))),
    )
    assert len(S) == len(B) == 14280 and not add(S, B)
    print(
        "Source and sink each have",
        len(S),
        "nonzero Laurent monomials. Their sum vanishes.",
    )
