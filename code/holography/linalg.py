"""Exact sparse linear algebra over Q.

Vectors are dicts {key: Fraction}; keys are hashable and totally ordered within one call
(tuples of ints).  No floating point anywhere.
"""
from fractions import Fraction

__all__ = ["Fraction", "vadd", "vscale", "clean", "Echelon", "rank", "kernel", "rel_rank"]


def clean(v):
    return {k: Fraction(c) for k, c in v.items() if c}


def vadd(acc, v, c=1):
    """acc += c*v (in place); returns acc."""
    if not c:
        return acc
    for k, x in v.items():
        y = acc.get(k, 0) + c * x
        if y:
            acc[k] = y
        else:
            acc.pop(k, None)
    return acc


def vscale(v, c):
    return {k: c * x for k, x in v.items()} if c else {}


class Echelon:
    """Incremental row echelon form.  piv[col] = row with row[col] = 1 and no earlier pivot columns.

    If track=True every stored row carries the combination (dict over insertion labels) of input
    vectors that produced it; reduce() then also returns the combination of the residue.
    """

    def __init__(self, track=False):
        self.piv = {}
        self.comb = {}
        self.track = track

    def reduce(self, v, label=None):
        row = clean(v)
        comb = {label: Fraction(1)} if (self.track and label is not None) else {}
        while True:
            hit = None
            for c in row:
                if c in self.piv:
                    hit = c
                    break
            if hit is None:
                return row, comb
            a = row[hit]
            vadd(row, self.piv[hit], -a)
            if self.track:
                vadd(comb, self.comb[hit], -a)

    def add(self, v, label=None):
        """Insert v; returns True if v was independent of the stored rows."""
        row, comb = self.reduce(v, label)
        if not row:
            return False, comb
        c0 = min(row)
        inv = 1 / row[c0]
        self.piv[c0] = vscale(row, inv)
        if self.track:
            self.comb[c0] = vscale(comb, inv)
        return True, comb

    def rank(self):
        return len(self.piv)

    def contains(self, v):
        row, _ = self.reduce(v)
        return not row


def rank(vectors):
    E = Echelon()
    for v in vectors:
        E.add(v)
    return E.rank()


def kernel(images):
    """images[i] = image of basis vector i.  Returns a basis of the kernel as dicts {i: coeff}."""
    E = Echelon(track=True)
    ker = []
    for i, v in enumerate(images):
        indep, comb = E.add(v, label=i)
        if not indep:
            ker.append(comb)
    return ker


def rel_rank(U, V):
    """dim (span U + span V) - dim span V."""
    E = Echelon()
    for v in V:
        E.add(v)
    r0 = E.rank()
    for u in U:
        E.add(u)
    return E.rank() - r0
