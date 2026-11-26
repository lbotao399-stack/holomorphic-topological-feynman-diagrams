"""Independent supercommutative algebra for superfields C^A(x, theta_1..theta_k) on R x C^H.

Generators are integers.  theta_i -> i (i = 0..k-1), theta'_i -> k + i (only for the functional
derivative factors), component letters c^A_S with holomorphic jet j (H-tuple) -> ids >= 2k.
Superfield (C odd):  C^A(theta) = sum_S theta_S c^A_S,  theta_S = theta_{s1} theta_{s2} ... (s1 < s2 < ...),
so c^A_S is odd iff |S| is even.  A polynomial is a dict {(odd ids sorted, even ids sorted): Fraction};
the theta's, having the smallest ids, stand on the left of every monomial.

Definitions used (S2 Sections 5-7, explorer's statement of the operator):
  Q0hat C^A = (1/2) f_{ABC} C^B C^C (odd derivation, Q0hat theta = 0, [Q0hat, d_a] = 0)
     =>  Q0hat c^A_S = (-1)^{|S|} [theta_S-coefficient of (1/2) f_{ABC} C^B C^C];
  delta C^B(theta')/delta C^A(theta) = delta^B_A prod_{i=1..k} (theta_i - theta'_i)
     =>  delta/delta C^A(theta) = sum_S F_S(theta) d_L/d c^A_S,  prod_i (theta_i - theta'_i) = sum_S F_S(theta) theta'_S;
  int d^k theta: coefficient of theta_0 theta_1 .. theta_{k-1} standing on the left.
"""
from fractions import Fraction as Fr
from itertools import combinations


class Alg:
    def __init__(self, k, H):
        self.k, self.H = k, H
        self.full = tuple(range(k))
        self.ids = {}
        self.info = {}          # id -> (S tuple, A, jet)
        self.par = {}           # id -> 1 odd, 0 even
        for i in range(2 * k):
            self.par[i] = 1
        self.subsets = [S for r in range(k + 1) for S in combinations(range(k), r)]
        self.F = self._Ffactors()
        self._cs = {}

    # ---------------------------------------------------------------- letters
    def letter(self, S, A, jet):
        S = tuple(S); jet = tuple(jet)
        key = (S, A, jet)
        if key not in self.ids:
            i = 2 * self.k + len(self.ids)
            self.ids[key] = i
            self.info[i] = key
            self.par[i] = 1 if len(S) % 2 == 0 else 0
        return self.ids[key]

    def gen(self, i, c=1):
        if self.par[i]:
            return {((i,), ()): Fr(c)}
        return {((), (i,)): Fr(c)}

    # ---------------------------------------------------------------- arithmetic
    @staticmethod
    def add(p, q, s=1):
        r = dict(p)
        for mn, c in q.items():
            v = r.get(mn, 0) + s * c
            if v:
                r[mn] = v
            else:
                r.pop(mn, None)
        return r

    @staticmethod
    def scale(p, s):
        if not s:
            return {}
        return {mn: c * s for mn, c in p.items()}

    @staticmethod
    def merge(a, b):
        """sorted concatenation of two sorted odd tuples: (tuple, sign) or None."""
        if not a:
            return b, 1
        if not b:
            return a, 1
        out = []
        i = j = 0
        inv = 0
        la, lb = len(a), len(b)
        while i < la and j < lb:
            x, y = a[i], b[j]
            if x < y:
                out.append(x); i += 1
            elif y < x:
                out.append(y); j += 1
                inv += la - i
            else:
                return None
        if i < la:
            out.extend(a[i:])
        if j < lb:
            out.extend(b[j:])
        return tuple(out), (-1 if inv & 1 else 1)

    def mul(self, p, q):
        r = {}
        mg = self.merge
        for (o1, e1), c1 in p.items():
            for (o2, e2), c2 in q.items():
                m = mg(o1, o2)
                if m is None:
                    continue
                o, s = m
                key = (o, tuple(sorted(e1 + e2)) if e1 and e2 else (e1 or e2))
                v = r.get(key, 0) + (c1 * c2 if s > 0 else -c1 * c2)
                if v:
                    r[key] = v
                else:
                    del r[key]
        return r

    def lder(self, p, l):
        r = {}
        if self.par[l]:
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
        return {mn: c for mn, c in r.items() if c}

    @staticmethod
    def letters(p):
        s = set()
        for (o, e) in p:
            s.update(o); s.update(e)
        return s

    def theta_part(self, o):
        n = 0
        while n < len(o) and o[n] < self.k:
            n += 1
        return n

    # ---------------------------------------------------------------- superfields, theta coefficients
    def theta_word(self, S, primed=False):
        p = {((), ()): Fr(1)}
        for i in S:
            p = self.mul(p, self.gen(i + (self.k if primed else 0)))
        return p

    def csuper(self, A, jet):
        key = (A, tuple(jet))
        if key not in self._cs:
            r = {}
            for S in self.subsets:
                r = self.add(r, self.mul(self.theta_word(S), self.gen(self.letter(S, A, jet))))
            self._cs[key] = r
        return self._cs[key]

    def thcoef(self, p, S):
        """coefficient of theta_S (on the left) in p; S a sorted tuple."""
        n = len(S)
        r = {}
        for (o, e), c in p.items():
            if self.theta_part(o) == n and o[:n] == S:
                r[(o[n:], e)] = c
        return r

    def decompose(self, p):
        """{S: p_S} with p = sum_S theta_S p_S."""
        out = {}
        for (o, e), c in p.items():
            n = self.theta_part(o)
            out.setdefault(o[:n], {})[(o[n:], e)] = c
        return out

    def _Ffactors(self):
        k = self.k
        p = {((), ()): Fr(1)}
        for i in range(k):
            p = self.mul(p, self.add(self.gen(i), self.gen(k + i), -1))
        F = {}
        for (o, e), c in p.items():
            th = tuple(x for x in o if x < k)
            pr = tuple(x - k for x in o if x >= k)
            # monomial o = theta-part then theta'-part (ids sorted): F_S theta'_S with S = pr
            F.setdefault(pr, {})[(th, ())] = c
        return F

    def sder(self, p, l):
        """delta/delta C^A(theta) restricted to the component letter l (with its jet)."""
        d = self.lder(p, l)
        if not d:
            return {}
        S = self.info[l][0]
        return self.mul(self.F[S], d)

    def dpar(self, p, a):
        r = {}
        for l in self.letters(p):
            if l < 2 * self.k:
                continue
            S, A, jet = self.info[l]
            j = list(jet); j[a] += 1
            r = self.add(r, self.mul(self.gen(self.letter(S, A, j)), self.lder(p, l)))
        return r


def so_structure_constants(N):
    """so(N), basis L_{ij} (i<j), (L_ij)_{ab} = delta_ia delta_jb - delta_ib delta_ja;
    [L_B, L_C] = sum_A f_{ABC} L_A.  Checks total antisymmetry and the Jacobi identity."""
    basis = [(i, j) for i in range(N) for j in range(i + 1, N)]
    n = len(basis)

    def mat(i, j):
        M = [[0] * N for _ in range(N)]
        M[i][j] = 1; M[j][i] = -1
        return M

    mats = [mat(*b) for b in basis]
    f = {}
    for B in range(n):
        for C in range(n):
            X, Y = mats[B], mats[C]
            Z = [[sum(X[a][c] * Y[c][b] - Y[a][c] * X[c][b] for c in range(N)) for b in range(N)] for a in range(N)]
            for A, (i, j) in enumerate(basis):
                if Z[i][j]:
                    f[(A, B, C)] = Fr(Z[i][j])
    for (A, B, C), v in f.items():
        assert f.get((B, C, A)) == v and f.get((B, A, C)) == -v
    for A in range(n):
        for B in range(n):
            for C in range(n):
                for D in range(n):
                    s = sum(f.get((A, B, E), 0) * f.get((E, C, D), 0) + f.get((A, C, E), 0) * f.get((E, D, B), 0)
                            + f.get((A, D, E), 0) * f.get((E, B, C), 0) for E in range(n))
                    assert s == 0
    return n, f


class Q0:
    def __init__(self, alg, n, f):
        self.S = alg
        self.n, self.f = n, f
        self.fA = {}
        for (A, B, C), v in f.items():
            self.fA.setdefault(A, []).append((B, C, v))
        self._q = {}

    def q0_letter(self, l):
        if l in self._q:
            return self._q[l]
        S_, A, jet = self.S.info[l]
        if any(jet):
            a = next(i for i in range(len(jet)) if jet[i])
            j = list(jet); j[a] -= 1
            r = self.S.dpar(self.q0_letter(self.S.letter(S_, A, j)), a)
        else:
            Rr = {}
            z = (0,) * self.S.H
            for (B, C, v) in self.fA.get(A, []):
                Rr = self.S.add(Rr, self.S.scale(self.S.mul(self.S.csuper(B, z), self.S.csuper(C, z)), v / 2))
            r = self.S.scale(self.S.thcoef(Rr, S_), (-1) ** len(S_))
        self._q[l] = r
        return r

    def __call__(self, p):
        r = {}
        for l in self.S.letters(p):
            if l < 2 * self.S.k:
                continue
            d = self.S.lder(p, l)
            if d:
                r = self.S.add(r, self.S.mul(self.q0_letter(l), d))
        return r
