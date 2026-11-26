"""Restricted seventh-power formula for the one-form source and cubic sink paths."""

from itertools import permutations
from random import Random
from frobenius import P, add, scale, dot, transverse

PERMS = [
    (p, (-1) ** sum(x > y for i, x in enumerate(p) for y in p[i + 1 :]))
    for p in permutations(range(5))
]


def endpoint(p, eta, incoming, z, kind):
    if kind == "source":
        return add(scale(dot(p, z), eta), scale(-dot(eta, z), p))
    out = [0] * 5
    for (i, j, k, l, m), sgn in PERMS:
        out[i] += sgn * z[j] * eta[k] * incoming[l] * p[m]
    return tuple(x % P for x in out)


class GammaPath:
    def __init__(self, ell, p, u, z, s, kind):
        self.ell, self.p, self.u, self.z, self.s, self.kind = ell, p, u, z, s, kind
        self.sigma = -1 if kind == "source" else 1
        self.zero = (0,) * 5
        self.full = 31
        self.prefix = [self.zero] * 32
        for m in range(1, 32):
            k = m & -m
            self.prefix[m] = add(self.prefix[m ^ k], p[k.bit_length() - 1])
        assert all(dot(q, v) == 0 for q, v in zip(p[1:], u[1:]))
        if kind == "source":
            assert dot(z, self.input(0)) == 0

    def input(self, mask):
        return add(self.ell, self.prefix[31 ^ mask])

    def insert(self, rows, key, v):
        rows[key] = add(rows.get(key, self.zero), v)

    def action(self, q, u, k, z, vector):
        # Unnegated Lie derivative for momentum -q on the K-mode.
        if vector:
            return add(scale(dot(u, k), z), scale(dot(q, z), u))
        return add(scale(dot(u, k), z), scale(-dot(u, z), q))

    def direct(self, rows):
        ans = {}
        for m, v in rows.items():
            k = self.input(m)
            self.insert(ans, m, scale(dot(self.s, self.prefix[31 ^ m]), v))
            for i in range(5):
                bit = 1 << i
                if m & bit:
                    continue
                q = self.p[i]
                u = self.u[i]
                if i == 0:
                    out = endpoint(q, u, k, v, self.kind)
                else:
                    before = bool(m & 1)
                    forward = (self.kind == "source") != before
                    alpha = dot(u, k)
                    if forward:
                        new = add(k, scale(-1, q))
                        out = add(scale(dot(new, v), u), scale(-alpha, v))
                    else:
                        out = add(scale(alpha, v), scale(-dot(u, v), k))
                    if before:
                        out = scale(-1, out)
                self.insert(ans, m | bit, out)
        return ans

    def fused(self):
        ordinary = {
            0: self.s,
            **{1 << i: scale(self.sigma, self.u[i]) for i in range(1, 5)},
        }
        gamma = {1: self.u[0]}
        for _ in range(6):
            n1 = {}
            n2 = {}
            for m, v in ordinary.items():
                q = self.prefix[m]
                self.insert(n1, m, scale(-dot(self.s, q), v))
                for j in range(1, 5):
                    if not m & (1 << j):
                        self.insert(
                            n1, m | 1 << j, scale(-self.sigma * dot(self.u[j], q), v)
                        )
            for m, v in gamma.items():
                k = scale(-1, self.prefix[m])
                self.insert(n2, m, scale(dot(self.s, k), v))
                for j in range(1, 5):
                    if not m & (1 << j):
                        new = self.action(self.p[j], self.u[j], k, v, False)
                        self.insert(n2, m | 1 << j, scale(self.sigma, new))
            ordinary, gamma = n1, n2
        k = self.input(0)
        out = {0: scale(-pow(dot(self.s, self.ell), 7, P), self.z)}
        for m, v in ordinary.items():
            self.insert(
                out, m, self.action(self.prefix[m], v, k, self.z, self.kind == "source")
            )
        for m, v in gamma.items():
            self.insert(out, m, endpoint(self.prefix[m], v, k, self.z, self.kind))
        return out

    def moment(self, d, use_fusion):
        values = self.fused() if use_fusion else {0: self.z}
        for _ in range(5 + d - 7 * int(use_fusion)):
            values = self.direct(values)
        return values.get(31, self.zero)


if __name__ == "__main__":
    rng = Random(78313)

    def vec():
        return tuple(rng.randrange(7) for _ in range(5))

    for kind in ["source", "sink"]:
        for _ in range(32):
            ell, s = vec(), vec()
            p = [vec() for _ in range(5)]
            u = [vec()] + [transverse(q, vec()) for q in p[1:]]
            k = ell
            for q in p:
                k = add(k, q)
            z = transverse(k, vec()) if kind == "source" else vec()
            t = GammaPath(ell, p, u, z, s, kind)
            for d in [2, 3, 4, 5, 6]:
                a = t.moment(d, False)
                b = t.moment(d, True)
                delta = add(a, scale(-1, b))
                if kind == "source":
                    assert all(
                        (x * delta[j] - ell[j] * delta[i]) % 7 == 0
                        for i, x in enumerate(ell)
                        for j in range(5)
                    )
                else:
                    assert not any(delta)
        print(kind, "160 fused path moments")
