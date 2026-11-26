"""Graded-commutative letter algebra for superfields C^A(theta_1..theta_k) on R x C^H, any k
(generalization of 6d-holomorphic-twist-note/code/salg.py, which is the case k = 2 with another
encoding of the letters).

Letters (t, A, j1, j2, j3):  t = 0: theta_A (A = 1..k);  t = 1 + mask: component of C^A along
theta_S, S = bits of mask, with holomorphic jet (j1, j2, j3).
C^A(theta) = sum_S theta_S c_S^A,  theta_S = theta_{s1} theta_{s2} ... (increasing), C odd, so
c_S is odd iff |S| is even.
delta C^B(theta') / delta C^A(theta) = delta_A^B prod_{i=1..k} (theta_i - theta'_i):
delta/delta C^A(theta) = sum_S F_S(theta) d_L/d c_S^A, with prod_i (theta_i - theta'_i) = sum_S F_S(theta) theta'_S.
Polynomial: dict {(odd_tuple, even_tuple): Fraction}, both tuples sorted.
"""
from fractions import Fraction as Fr
from itertools import product


class SAlg:
    def __init__(self, k):
        self.k = k
        self.full = (1 << k) - 1
        self._F = {mask: self._Ffactor(mask) for mask in range(1 << k)}

    # ---------------------------------------------------------------- basics
    def is_odd(self, l):
        if l[0] == 0:
            return True
        return bin(l[0] - 1).count('1') % 2 == 0

    @staticmethod
    def const(x=1):
        return {((), ()): Fr(x)}

    def gen(self, l):
        return {((l,), ()): Fr(1)} if self.is_odd(l) else {((), (l,)): Fr(1)}

    @staticmethod
    def add(p, q, s=1):
        r = dict(p)
        for m, c in q.items():
            v = r.get(m, 0) + s * c
            if v:
                r[m] = v
            else:
                r.pop(m, None)
        return r

    @staticmethod
    def scale(p, s):
        return {m: c * s for m, c in p.items()} if s else {}

    @staticmethod
    def _merge_odd(a, b):
        if not a:
            return b, 1
        if not b:
            return a, 1
        out, i, j, inv = [], 0, 0, 0
        la, lb = len(a), len(b)
        while i < la and j < lb:
            if a[i] < b[j]:
                out.append(a[i]); i += 1
            elif a[i] > b[j]:
                out.append(b[j]); j += 1
                inv += la - i
            else:
                return None
        out.extend(a[i:]); out.extend(b[j:])
        return tuple(out), (-1 if inv & 1 else 1)

    def mul(self, p, q):
        r = {}
        for (o1, e1), c1 in p.items():
            for (o2, e2), c2 in q.items():
                mo = self._merge_odd(o1, o2)
                if mo is None:
                    continue
                o, s = mo
                key = (o, tuple(sorted(e1 + e2)))
                v = r.get(key, 0) + s * c1 * c2
                if v:
                    r[key] = v
                else:
                    r.pop(key, None)
        return r

    def lderiv(self, p, l):
        r = {}
        if self.is_odd(l):
            for (o, e), c in p.items():
                if l in o:
                    i = o.index(l)
                    key = (o[:i] + o[i + 1:], e)
                    r[key] = r.get(key, 0) + (-c if i & 1 else c)
        else:
            for (o, e), c in p.items():
                n = e.count(l)
                if n:
                    i = e.index(l)
                    key = (o, e[:i] + e[i + 1:])
                    r[key] = r.get(key, 0) + n * c
        return {m: c for m, c in r.items() if c}

    @staticmethod
    def letters_of(p):
        s = set()
        for (o, e) in p:
            s.update(o); s.update(e)
        return s

    @staticmethod
    def shift(l, a):
        j = list(l[2:]); j[a] += 1
        return (l[0], l[1]) + tuple(j)

    def dpar(self, p, a):
        r = {}
        for l in self.letters_of(p):
            if l[0] == 0:
                continue
            r = self.add(r, self.mul(self.gen(self.shift(l, a)), self.lderiv(p, l)))
        return r

    def theta_word(self, mask):
        out = self.const(1)
        for i in range(self.k):
            if mask >> i & 1:
                out = self.mul(out, self.gen((0, i + 1, 0, 0, 0)))
        return out

    def csuper(self, A, j=(0, 0, 0)):
        r = {}
        for mask in range(1 << self.k):
            r = self.add(r, self.mul(self.theta_word(mask), self.gen((1 + mask, A) + tuple(j))))
        return r

    def theta_coeff(self, p, mask):
        want = tuple((0, i + 1, 0, 0, 0) for i in range(self.k) if mask >> i & 1)
        n = len(want)
        r = {}
        for (o, e), c in p.items():
            nth = sum(1 for l in o if l[0] == 0)
            if nth == n and o[:n] == want:
                r[(o[n:], e)] = c
        return r

    def _Ffactor(self, mask):
        """F_S(theta) from prod_i (theta_i - theta'_i) = sum_S F_S theta'_S (theta'_S on the right)."""
        k = self.k
        # term: theta'_i for i in S (sign -1 each), theta_i otherwise; word in order i = 1..k;
        # move the theta' letters to the right: inversions = #{i < j: i in S, j not in S}
        S = [i for i in range(k) if mask >> i & 1]
        notS = [i for i in range(k) if not mask >> i & 1]
        inv = sum(1 for i in S for j in notS if i < j)
        sign = (-1) ** (len(S) + inv)
        return self.scale(self.theta_word(self.full ^ mask), sign)

    def sderiv(self, p, l):
        d = self.lderiv(p, l)
        return self.mul(self._F[l[0] - 1], d) if d else {}


def structure_constants(name):
    from salg import structure_constants as sc
    return sc(name)
