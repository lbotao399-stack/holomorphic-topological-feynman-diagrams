"""Degree-eight bound on five theta Wick tensors by formal logarithms of relative transports."""

from mixed_transport_symbol import mono, add, neg, prod2, prod, vadd, vneg, varvec, eps
from itertools import product as cart
from collections import Counter
from math import comb, factorial, prod as iprod


def verify(N, collinear):
    zero = [{} for _ in range(5)]

    def const(v):
        return [v] + [zero] * N

    def plus(s, t):
        return [vadd(a, b) for a, b in zip(s, t)]

    def minus(s, t):
        return plus(s, [vneg(a) for a in t])

    def mulc(v, c):
        return [{k: x * c for k, x in p.items()} for p in v]

    def scalc(p, c):
        return {k: x * c for k, x in p.items()}

    l = varvec(0)
    l[4] = neg(add(*l[:4]))
    K = [[mono(256 ** (5 + 5 * i + j)) for j in range(5)] for i in range(5)]
    K[4][4] = neg(add(*(K[i][i] for i in range(4))))

    def L(v):
        return [prod2(a, b) for a, b in zip(l, v)]

    def kt(v):
        return [add(*(prod2(a, b) for a, b in zip(r, v))) for r in K]

    def apply(s, op, k, sgn):
        t = s[k]
        for i in range(1, N - k + 1):
            t = op(t)
            yield i + k, mulc(t, sgn**i * comb(i + k, i))

    def exp(s, op, sgn, base=True):
        out = s.copy() if base else [zero] * (N + 1)
        for k, v in enumerate(s):
            if not any(v):
                continue
            for i, w in apply(s, op, k, sgn):
                out[i] = vadd(out[i], w)
        return out

    def W(s):
        return exp(s, L, 1)

    def U(s):
        return exp(s, kt, 1)

    def D(s):
        return exp(s, kt, -1)

    def n(s):
        return exp(s, L, 1, False)

    def r(s):
        return exp(s, L, -1, False)

    def z(s):
        return plus(n(s), r(s))

    M = const([mono(0)] * 5)
    V = const([mono(256**30)] * 5 if collinear else varvec(30))
    eta = const([mono(0), {}, {}, {}, {}])
    F = W(V)
    Q = minus(F, V)
    E = minus(M, Q)
    p = minus(D(E), M)
    up = U(p)
    ue = U(eta)

    def det(args, d):
        acc = Counter()
        slots = [[i for i in range(d + 1) if any(s[i])] for s in args]
        for ix in cart(*slots):
            if sum(ix) != d:
                continue
            mult = factorial(d) // iprod(factorial(i) for i in ix)
            v = eps(*(s[i] for s, i in zip(args, ix)))
            acc.update(scalc(v, mult))
        return {k: x for k, x in acc.items() if x}

    def ds(a, b, c, d, op=lambda x: x):
        return add(det([eta, op(up), a, b, c], d), neg(det([p, op(ue), a, b, c], d)))

    def t(a, b, d):
        return det([a, D(b), eta, M, p], d)

    def trace(i):
        if not i:
            return {}
        vals = [mono(0) for _ in range(5)]
        for _ in range(i):
            vals = [prod2(a, b) for a, b in zip(vals, l)]
        return add(*vals)

    def trace_mul(fn, d, odd):
        return add(
            *(
                scalc(prod2(trace(i), fn(d - i)), comb(d, i) * (2 if odd else 1))
                for i in range(1, d + 1)
                if not odd or i % 2
            )
        )

    for d in [5] if collinear else [3, 4]:
        source = add(
            trace_mul(lambda k: ds(Q, M, V, k), d, False),
            ds(z(E), M, V, d),
            neg(ds(Q, M, V, d, z)),
            ds(n(M), M, Q, d),
        )
        sink = add(
            t(V, z(E), d),
            neg(t(z(M), F, d)),
            t(Q, n(M), d),
            t(r(E), Q, d),
            trace_mul(lambda k: t(V, F, k), d, True),
        )
        assert (
            not add(source, sink)
            and len(source) == len(sink) == {3: 414, 4: 3060, 5: 0}[d]
        )
        print(
            "Collinear" if collinear else "Independent",
            "derivative",
            d,
            "source and sink monomials",
            len(source),
            len(sink),
        )


if __name__ == "__main__":
    verify(4, False)
    verify(5, True)
