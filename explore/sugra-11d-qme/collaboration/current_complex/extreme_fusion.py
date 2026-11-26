"""Dual bracket coefficients for the extreme Taylor profile, over the rationals."""

from collections import defaultdict
from functools import cache
from fractions import Fraction
from math import comb, lcm


def comp(n, d):
    if d == 1:
        yield (n,)
        return
    for j in range(n + 1):
        for t in comp(n - j, d - 1):
            yield (j,) + t


def shift(a, i, j):
    b = list(a)
    b[i] -= 1
    b[j] += 1
    return tuple(b)


class Jet:
    def __init__(self, d):
        self.basis = []
        self.pos = {}
        self.weight = []
        for a in comp(d, 5):
            for v in range(5):
                if a[v] and v == 4:
                    continue
                p = {(a, v): 1}
                if a[v]:
                    p[shift(a, v, 4), 4] = -1
                self.pos[a, v] = len(self.basis)
                self.basis.append(p)
                self.weight.append(tuple(int(k == v) - x for k, x in enumerate(a)))

    def coords(self, p):
        return {self.pos[m]: x for m, x in p.items() if m in self.pos and x}

    def action(self, t, i, j):
        out = defaultdict(int)
        for (a, v), x in self.basis[t].items():
            if v == j:
                out[a, i] += x
            if a[i]:
                out[shift(a, i, j), v] -= x * (1 + a[j])
        return self.coords(out)

    def bracket(self, a, b):
        c = defaultdict(int)
        for x, y, sg in [(a, b, 1), (b, a, -1)]:
            for (m, v), s in x.items():
                for (n, w), t in y.items():
                    if n[v]:
                        h = tuple(
                            u + z - int(i == v) for i, (u, z) in enumerate(zip(m, n))
                        )
                        coeff = s * t * sg
                        for u, z in zip(m, h):
                            coeff *= comb(z, u)
                        c[h, w] += coeff
        return self.coords(c)


def add(a, b, c):
    for i, x in b.items():
        y = a.get(i, 0) + c * x
        if y:
            a[i] = y
        else:
            a.pop(i, None)


class Dual:
    """Dual of the high vector jet, tensored with the tangent and det^-3."""

    def __init__(self, d):
        self.jets = Jet(d)
        self.groups = defaultdict(list)
        for t, w in enumerate(self.jets.weight):
            for g in range(5):
                self.groups[
                    tuple(-x - 3 + int(i == g) for i, x in enumerate(w))
                ].append((t, g))
        self.acts = {}
        for i in range(5):
            for j in range(5):
                if abs(i - j) != 1:
                    continue
                rows = defaultdict(dict)
                for t in range(len(self.jets.basis)):
                    for h, c in self.jets.action(t, i, j).items():
                        rows[h][t] = -c
                self.acts[i, j] = rows

    def act(self, v, i, j):
        z = {}
        for (t, g), c in v.items():
            for h, d in self.acts[i, j][t].items():
                add(z, {(h, g): d}, c)
            if g == j:
                add(z, {(t, i): 1}, c)
        return z

    def highest(self, wt):
        slots = self.groups[wt]
        pivots = {}
        for k in slots:
            v = {k: 1}
            out = {}
            for i in range(4):
                for s, x in self.act(v, i, i + 1).items():
                    out[i, s] = x
            while out:
                p = min(out)
                c = out[p]
                if p not in pivots:
                    pivots[p] = (v, out)
                    break
                old, test = pivots[p]
                coef = -c * Fraction(1, test[p])
                add(out, test, coef)
                add(v, old, coef)
            else:
                return v
        raise ValueError(wt)


def data():
    # Quadratic vectors acting by the Taylor bracket on seventh-order jets.
    q = Jet(2)
    a = Dual(7)
    b = Dual(8)
    wi = [(5, -3, -3, -3, -4), (4, -2, -3, -3, -4), (4, -3, -3, -3, -3)]
    wo = [(w[0] + 1,) + w[1:] for w in wi]
    tops = {w: a.highest(w) for w in wi}

    @cache
    def descendants(w):
        """Independent lowering words on the three highest vectors at weight w."""
        roots = []
        if w in tops:
            roots.append(((wi.index(w), ()), tops[w]))
        for i in range(4):
            v = shift(w, i + 1, i)
            if not any(
                all(sum(c - d for c, d in zip(h[:k], v[:k])) >= 0 for k in range(1, 5))
                for h in wi
            ):
                continue
            if not a.groups.get(v):
                continue
            for (t, word), old in descendants(v):
                new = a.act(old, i + 1, i)
                if new:
                    roots.append(((t, word + (i,)), new))
        basis = {}
        accepted = []
        for name, col in roots:
            pivcol = col.copy()
            while pivcol:
                p = min(pivcol)
                c = pivcol[p]
                if p not in basis:
                    basis[p] = pivcol
                    accepted.append((name, col))
                    break
                add(pivcol, basis[p], -c * Fraction(1, basis[p][p]))
        assert len(accepted) == len(a.groups[w]), (w, len(accepted), len(a.groups[w]))
        return accepted

    for h, w in enumerate(wo):
        terms = []
        top = b.highest(w)
        for i, weight in enumerate(q.weight):
            v = tuple(x + y for x, y in zip(w, weight))
            rhs = {}
            for t, g in a.groups.get(v, []):
                c = sum(
                    top.get((k, g), 0) * s
                    for k, s in b.jets.bracket(q.basis[i], a.jets.basis[t]).items()
                )
                if c:
                    rhs[t, g] = c
            if not rhs:
                continue
            piv = {}
            for name, y in descendants(v):
                c = {name: 1}
                row = y.copy()
                while row:
                    key = min(row)
                    if key not in piv:
                        piv[key] = row, c
                        break
                    rowp, cp = piv[key]
                    f = -row[key] * Fraction(1, rowp[key])
                    add(row, rowp, f)
                    add(c, cp, f)
            co = {}
            while rhs:
                key = min(rhs)
                rowp, cp = piv[key]
                f = rhs[key] * Fraction(1, rowp[key])
                add(rhs, rowp, -f)
                add(co, cp, f)
            for (t, word), c in co.items():
                terms.append((i, t, Fraction(c), word))
        denominator = lcm(*(c.denominator for _, _, c, _ in terms))
        for i, t, c, word in terms:
            print(h, i, t, int(c * denominator), len(word), *word)


if __name__ == "__main__":
    data()
