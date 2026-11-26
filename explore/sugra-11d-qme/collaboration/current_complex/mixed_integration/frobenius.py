"""Characteristic-seven path moments from one transverse restricted adjoint."""

from random import Random

P = 7


def add(x, y):
    return tuple((a + b) % P for a, b in zip(x, y))


def scale(a, x):
    return tuple(a * b % P for b in x)


def dot(x, y):
    return sum(a * b for a, b in zip(x, y)) % P


class OrdinaryPath:
    def __init__(self, ell, p, u, z, observer):
        self.ell, self.p, self.u, self.s = ell, p, u, observer
        self.n = len(p)
        self.full = (1 << self.n) - 1
        self.zero = (0,) * len(ell)
        self.prefix = [self.zero] * (1 << self.n)
        for t in range(1, 1 << self.n):
            bit = t & -t
            self.prefix[t] = add(self.prefix[t ^ bit], p[bit.bit_length() - 1])
        self.z = z
        assert all(dot(a, b) == 0 for a, b in zip(p, u))
        assert dot(z, add(ell, self.prefix[-1])) == 0

    def total(self, terms, key, vector):
        terms[key] = add(terms.get(key, self.zero), vector)

    def step(self, terms):
        # Fourier input has remaining momentum ell + unconsumed external momenta.
        result = {}
        for subset, v in terms.items():
            rem = self.full ^ subset
            shift = self.prefix[rem]
            self.total(result, subset, scale(dot(self.s, shift), v))
            while rem:
                bit = rem & -rem
                rem ^= bit
                j = bit.bit_length() - 1
                incoming = add(self.ell, self.prefix[self.full ^ subset])
                image = add(
                    scale(-dot(self.u[j], incoming), v),
                    scale(-dot(self.p[j], v), self.u[j]),
                )
                self.total(result, subset | bit, image)
        return result

    def restricted_field(self):
        # Z^7 applied to the coordinate functions, Z=s-sum e_i u_i exp(-p_i.y).
        terms = {0: self.s, **{1 << i: scale(-1, u) for i, u in enumerate(self.u)}}
        for _ in range(P - 1):
            result = {}
            for subset, v in terms.items():
                q = self.prefix[subset]
                self.total(result, subset, scale(-dot(self.s, q), v))
                for j, u in enumerate(self.u):
                    if not subset & (1 << j):
                        self.total(result, subset | 1 << j, scale(dot(u, q), v))
            terms = result
        assert all(dot(self.prefix[subset], v) == 0 for subset, v in terms.items())
        return terms

    def moment(self, d, fused):
        if not fused:
            values = {0: self.z}
            for _ in range(self.n + d):
                values = self.step(values)
            return values[self.full]
        assert self.n + d >= 7
        values = {0: scale(-pow(dot(self.s, self.ell), 7, P), self.z)}
        incoming = add(self.ell, self.prefix[-1])
        for subset, w in self.restricted_field().items():
            v = add(
                scale(dot(w, incoming), self.z),
                scale(dot(self.prefix[subset], self.z), w),
            )
            self.total(values, subset, v)
        for _ in range(self.n + d - 7):
            values = self.step(values)
        return values.get(self.full, self.zero)


def transverse(p, u):
    i = next((i for i, a in enumerate(p) if a % 7), None)
    if i is None:
        return u
    out = list(u)
    out[i] = (out[i] - dot(p, u) * pow(p[i], -1, 7)) % 7
    return tuple(out)


def test():
    rng = Random(7943)

    def vec():
        return tuple(rng.randrange(7) for _ in range(5))

    for _ in range(25):
        ell, s = vec(), vec()
        p = [vec() for i in range(5)]
        u = [transverse(q, vec()) for q in p]
        q = ell
        for k in p:
            q = add(q, k)
        z = transverse(q, vec())
        path = OrdinaryPath(ell, p, u, z, s)
        for d in range(2, 7):
            assert path.moment(d, True) == path.moment(d, False)
    # Separated ordinary inputs but unrestricted observer, terminal vector and initial momentum.
    for _ in range(25):
        h = 2
        ell, s = vec(), vec()
        p = []
        u = []
        for i in range(5):
            q = vec()
            v = vec()
            q = tuple(x if (i < h) == (j < h) else 0 for j, x in enumerate(q))
            v = tuple(x if (i < h) == (j < h) else 0 for j, x in enumerate(v))
            p.append(q)
            u.append(transverse(q, v))
        q = ell
        for k in p:
            q = add(q, k)
        path = OrdinaryPath(ell, p, u, transverse(q, vec()), s)
        for d in [2, 3]:
            assert path.moment(d, False) == path.zero
        for d in [2, 3, 4, 5, 6]:
            assert path.moment(d, False) == path.moment(d, True)
    print("250 projected path identities over F7 and 50 separated zeros")


if __name__ == "__main__":
    test()
