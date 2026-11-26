"""Graded-commutative jet algebra for the strict model, exact integer arithmetic.

Letters (0-based indices i, k):
  ('X', i, a)  odd      ('n', 0, a)  nu, even     ('g', i, a)  gamma_i, even     ('b', 0, a)  beta, odd
Polynomial: dict {monomial: int}; monomial = tuple of letter ids sorted increasingly
(even ids may repeat, odd ids may not).  The product of the letters in tuple order equals the monomial.

Q (odd derivation, acts from the left, Q(ab) = Q(a) b + (-1)^{|a|} a Q(b)):
  Q X^i     = s_X * X^j X^i_{e_j}                      (s_X = -1 : derived sign; s_X = +1 : the stated sign)
  Q nu      = X^i_{e_i} - nu X^i_{e_i} - X^i nu_{e_i}
  Q gamma_i = beta_{e_i} + X^k (gamma_{k,e_i} - gamma_{i,e_k})
  Q beta    = X^i X^j gamma_{j,e_i}
  Q(phi_a)  = D^a (Q phi)
Quotient model (quot=True): no nu; X^0_a with a_0 >= 1 is replaced by -sum_{i>=1} X^i_{a-e_0+e_i}.
"""
from math import gcd, comb
from bisect import bisect_left
import itertools

TYPES = ('X', 'n', 'g', 'b')
PAR = {'X': 1, 'n': 0, 'g': 0, 'b': 1}


def multi_indices(n, m):
    if n == 1:
        yield (m,)
        return
    for a in range(m, -1, -1):
        for rest in multi_indices(n - 1, m - a):
            yield (a,) + rest


def addv(a, b):
    return tuple(x + y for x, y in zip(a, b))


def unit(n, k):
    return tuple(int(j == k) for j in range(n))


class Alg:
    def __init__(self, n, quot=False, sX=-1):
        self.n = n
        self.quot = quot
        self.sX = sX
        self.ids = {}
        self.info = []
        self.par = []
        self.Gw = []
        self.tw = []
        self._Q = {}
        self._red = {}
        self.zero = (0,) * n

    # ---------- letters ----------
    def L(self, typ, i, a):
        key = (typ, i, tuple(a))
        r = self.ids.get(key)
        if r is None:
            r = len(self.info)
            self.ids[key] = r
            self.info.append(key)
            self.par.append(PAR[typ])
            m = sum(a)
            n = self.n
            if typ == 'X':
                G = m - 1
                t = tuple(int(j == i) - a[j] for j in range(n))
            elif typ == 'n':
                G = m
                t = tuple(-a[j] for j in range(n))
            elif typ == 'g':
                G = m
                t = tuple(-int(j == i) - a[j] for j in range(n))
            else:
                G = m - 1
                t = tuple(-a[j] for j in range(n))
            self.Gw.append(G)
            self.tw.append(t)
        return r

    def is_basis(self, lid):
        typ, i, a = self.info[lid]
        if not self.quot:
            return True
        if typ == 'n':
            return False
        if typ == 'X' and i == 0 and a[0] >= 1:
            return False
        return True

    def name(self, lid):
        typ, i, a = self.info[lid]
        nm = {'X': 'X%d' % (i + 1), 'n': 'nu', 'g': 'g%d' % (i + 1), 'b': 'b'}[typ]
        if sum(a):
            nm += '_' + ''.join(str(x) for x in a)
        return nm

    # ---------- monomial arithmetic ----------
    def mono_mul(self, m1, m2):
        """returns (sign, m) with m1*m2 = sign*m, or (0, None)."""
        if not m1:
            return 1, m2
        if not m2:
            return 1, m1
        par = self.par
        o2 = [x for x in m2 if par[x]]
        sign = 1
        if o2:
            o1 = [x for x in m1 if par[x]]
            if o1:
                s1 = set(o1)
                for x in o2:
                    if x in s1:
                        return 0, None
                # inversions: pairs (a in o1, b in o2) with a > b
                inv = 0
                o2s = o2  # sorted already
                for a in o1:
                    inv += bisect_left(o2s, a)
                if inv & 1:
                    sign = -1
        return sign, tuple(sorted(m1 + m2))

    def pmul(self, p, q):
        r = {}
        for m1, c1 in p.items():
            for m2, c2 in q.items():
                s, m = self.mono_mul(m1, m2)
                if s:
                    v = r.get(m, 0) + s * c1 * c2
                    if v:
                        r[m] = v
                    else:
                        r.pop(m, None)
        return r

    def padd(self, p, q, c=1):
        r = dict(p)
        for m, v in q.items():
            x = r.get(m, 0) + c * v
            if x:
                r[m] = x
            else:
                r.pop(m, None)
        return r

    def letter(self, typ, i, a):
        return {(self.L(typ, i, a),): 1}

    def prod(self, *ps):
        r = {(): 1}
        for p in ps:
            r = self.pmul(r, p)
        return r

    # ---------- total derivative ----------
    def shift(self, lid, k):
        typ, i, a = self.info[lid]
        b = list(a)
        b[k] += 1
        return self.L(typ, i, tuple(b))

    def Dmono(self, m, k):
        r = {}
        par = self.par
        for j in range(len(m)):
            if j > 0 and m[j] == m[j - 1] and not par[m[j]]:
                continue  # repeated even letter: handled with multiplicity below
            lid = m[j]
            mult = 1
            if not par[lid]:
                mult = m.count(lid)
            rest = m[:j] + m[j + 1:] if par[lid] else tuple(x for idx, x in enumerate(m) if idx != j)
            nl = self.shift(lid, k)
            # m = eps * lid * rest, eps from moving lid to front
            eps = 1
            if par[lid]:
                before = sum(1 for x in m[:j] if par[x])
                if before & 1:
                    eps = -1
            s, mm = self.mono_mul((nl,), rest)
            if s:
                v = r.get(mm, 0) + eps * s * mult
                if v:
                    r[mm] = v
                else:
                    r.pop(mm, None)
        return r

    def D(self, p, k):
        r = {}
        for m, c in p.items():
            for mm, v in self.Dmono(m, k).items():
                x = r.get(mm, 0) + c * v
                if x:
                    r[mm] = x
                else:
                    r.pop(mm, None)
        return r

    # ---------- reduction to the quotient basis ----------
    def red_letter(self, lid):
        r = self._red.get(lid)
        if r is not None:
            return r
        typ, i, a = self.info[lid]
        if not self.quot:
            r = {(lid,): 1}
        elif typ == 'n':
            r = {}
        elif typ == 'X' and i == 0 and a[0] >= 1:
            r = {}
            for j in range(1, self.n):
                b = list(a)
                b[0] -= 1
                b[j] += 1
                r[(self.L('X', j, tuple(b)),)] = -1
        else:
            r = {(lid,): 1}
        self._red[lid] = r
        return r

    def reduce(self, p):
        if not self.quot:
            return p
        out = {}
        for m, c in p.items():
            if all(self.is_basis(x) for x in m):
                v = out.get(m, 0) + c
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
                continue
            acc = {(): c}
            for x in m:
                acc = self.pmul(acc, self.red_letter(x))
                if not acc:
                    break
            out = self.padd(out, acc)
        return out

    # ---------- Q ----------
    def Qbase(self, typ, i):
        n = self.n
        z = self.zero
        e = [unit(n, k) for k in range(n)]
        if typ == 'X':
            r = {}
            for j in range(n):
                r = self.padd(r, self.prod(self.letter('X', j, z), self.letter('X', i, e[j])), self.sX)
            return r
        if typ == 'n':
            r = {}
            for j in range(n):
                r = self.padd(r, self.letter('X', j, e[j]))
                r = self.padd(r, self.prod(self.letter('n', 0, z), self.letter('X', j, e[j])), -1)
                r = self.padd(r, self.prod(self.letter('X', j, z), self.letter('n', 0, e[j])), -1)
            return r
        if typ == 'g':
            r = self.letter('b', 0, e[i])
            for k in range(n):
                r = self.padd(r, self.prod(self.letter('X', k, z), self.letter('g', k, e[i])))
                r = self.padd(r, self.prod(self.letter('X', k, z), self.letter('g', i, e[k])), -1)
            return r
        if typ == 'b':
            r = {}
            for a in range(n):
                for b in range(n):
                    r = self.padd(r, self.prod(self.letter('X', a, z), self.letter('X', b, z), self.letter('g', b, e[a])))
            return r
        raise ValueError

    def Qfull_letter(self, lid):
        """Q of a letter in the full algebra (no reduction)."""
        key = ('full', lid)
        r = self._Q.get(key)
        if r is not None:
            return r
        typ, i, a = self.info[lid]
        if sum(a) == 0:
            r = self.Qbase(typ, i)
        else:
            k = next(j for j in range(self.n) if a[j] > 0)
            b = list(a)
            b[k] -= 1
            r = self.D(self.Qfull_letter(self.L(typ, i, tuple(b))), k)
        self._Q[key] = r
        return r

    def Qletter(self, lid):
        if not self.quot:
            return self.Qfull_letter(lid)
        key = ('red', lid)
        r = self._Q.get(key)
        if r is None:
            r = self.reduce(self.Qfull_letter(lid))
            self._Q[key] = r
        return r

    def Qmono(self, m):
        """Q of a monomial (in the basis of the current algebra)."""
        par = self.par
        r = {}
        nodd = 0
        prev = None
        for j, lid in enumerate(m):
            if par[lid]:
                eps = -1 if (nodd & 1) else 1
                nodd += 1
                rest = m[:j] + m[j + 1:]
                mult = 1
            else:
                if lid == prev:
                    continue
                mult = m.count(lid)
                eps = 1
                jj = m.index(lid)
                rest = m[:jj] + m[jj + 1:]
            prev = lid
            ql = self.Qletter(lid)
            for mm, c in ql.items():
                s, mo = self.mono_mul(mm, rest)
                if s:
                    v = r.get(mo, 0) + eps * s * c * mult
                    if v:
                        r[mo] = v
                    else:
                        r.pop(mo, None)
        return r

    def Q(self, p):
        r = {}
        for m, c in p.items():
            for mm, v in self.Qmono(m).items():
                x = r.get(mm, 0) + c * v
                if x:
                    r[mm] = x
                else:
                    r.pop(mm, None)
        return r

    # ---------- derivations given on letters ----------
    def apply_der(self, p, f, odd):
        """derivation determined by f(lid) -> poly (already in the algebra basis); odd: parity of derivation."""
        par = self.par
        r = {}
        for m, c in p.items():
            nodd = 0
            prev = None
            for j, lid in enumerate(m):
                if par[lid]:
                    eps = -1 if (nodd & 1) else 1  # sign (-1)^{|letter| * #odd letters before}, independent of |f|
                    nodd += 1
                    rest = m[:j] + m[j + 1:]
                    mult = 1
                else:
                    if lid == prev:
                        continue
                    mult = m.count(lid)
                    eps = 1
                    jj = m.index(lid)
                    rest = m[:jj] + m[jj + 1:]
                prev = lid
                fl = f(lid)
                for mm, v in fl.items():
                    s, mo = self.mono_mul(mm, rest)
                    if s:
                        x = r.get(mo, 0) + eps * s * v * c * mult
                        if x:
                            r[mo] = x
                        else:
                            r.pop(mo, None)
        return r

    def Dq(self, p, k):
        """total derivative followed by reduction to the quotient basis."""
        return self.reduce(self.D(p, k))

    def weight(self, m):
        G = 0
        t = [0] * self.n
        for x in m:
            G += self.Gw[x]
            tx = self.tw[x]
            for j in range(self.n):
                t[j] += tx[j]
        return G, tuple(t)

    def Ndeg(self, m):
        return sum(1 for x in m if self.info[x][0] in ('X', 'b'))

    def wdeg(self, m):
        return sum(1 for x in m if self.info[x][0] in ('g', 'b'))

    def pstr(self, p, maxterms=40):
        items = sorted(p.items(), key=lambda kv: (len(kv[0]), [self.info[x] for x in kv[0]]))
        out = []
        for m, c in items[:maxterms]:
            nm = '*'.join(self.name(x) for x in m) if m else '1'
            out.append(('%+d ' % c) + nm)
        s = ' '.join(out)
        if len(items) > maxterms:
            s += ' ... (%d terms)' % len(items)
        return s if s else '0'


# ---------- basis enumeration for C_quot (w = 0, 1) ----------

def x_letters_quot(A, maxorder):
    """basis X letters of the quotient up to the given jet order, grouped by G = order-1."""
    n = A.n
    levels = {}
    for m in range(maxorder + 1):
        lst = []
        for a in multi_indices(n, m):
            for i in range(n):
                if i == 0 and a[0] >= 1:
                    continue
                lst.append(A.L('X', i, a))
        levels[m - 1] = sorted(lst)
    return levels


def xsets(A, g, levels):
    """all sets of basis X letters with total G = g. returns list of sorted tuples."""
    n = A.n
    out = []
    maxlev = max(levels)
    levs = sorted(levels, reverse=True)  # high G first

    def rec(li, budget, chosen):
        lev = levs[li]
        if lev == -1:
            k = -budget
            if 0 <= k <= len(levels[-1]):
                for c in itertools.combinations(levels[-1], k):
                    out.append(tuple(sorted(chosen + list(c))))
            return
        if lev == 0:
            for k in range(len(levels[0]) + 1):
                for c in itertools.combinations(levels[0], k):
                    rec(li + 1, budget, chosen + list(c))
            return
        kmax = (budget + n) // lev if budget + n >= 0 else -1
        for k in range(0, min(kmax, len(levels[lev])) + 1):
            for c in itertools.combinations(levels[lev], k):
                rec(li + 1, budget - k * lev, chosen + list(c))

    rec(0, g, [])
    return out


def wletters_quot(A, maxorder):
    """gamma and beta letters with G <= maxorder."""
    n = A.n
    out = []
    for m in range(maxorder + 2):
        for a in multi_indices(n, m):
            out.append(A.L('b', 0, a))
            if m <= maxorder:
                for i in range(n):
                    out.append(A.L('g', i, a))
    return out


def quot_space(A, w, G):
    """all basis monomials of C_quot with given w (0 or 1) and G, grouped: dict (t, N) -> sorted list."""
    n = A.n
    maxorder = G + n + 2 + 1
    levels = x_letters_quot(A, maxorder)
    groups = {}
    if w == 0:
        for xs in xsets(A, G, levels):
            Gm, t = A.weight(xs)
            assert Gm == G
            groups.setdefault((t, len(xs)), []).append(xs)
    elif w == 1:
        cache = {}
        for wl in wletters_quot(A, maxorder):
            gw = A.Gw[wl]
            g = G - gw
            if g < -n:
                continue
            if g not in cache:
                cache[g] = xsets(A, g, levels)
            for xs in cache[g]:
                m = tuple(sorted(xs + (wl,)))
                Gm, t = A.weight(m)
                assert Gm == G
                groups.setdefault((t, A.Ndeg(m)), []).append(m)
    else:
        raise ValueError
    for k in groups:
        groups[k].sort()
    return groups


# ---------- exact sparse echelon over Z ----------

class Echelon:
    """incremental row-echelon form of integer vectors (dict col -> int); exact."""

    def __init__(self):
        self.piv = {}

    def insert(self, v):
        v = {k: x for k, x in v.items() if x}
        piv = self.piv
        while v:
            c = max(v)
            p = piv.get(c)
            if p is None:
                g = 0
                for x in v.values():
                    g = gcd(g, x)
                    if g == 1:
                        break
                if g > 1:
                    v = {k: x // g for k, x in v.items()}
                if v[c] < 0:
                    v = {k: -x for k, x in v.items()}
                piv[c] = v
                return True
            a = p[c]
            b = v[c]
            g = gcd(a, b)
            a //= g
            b //= g
            if a != 1:
                nv = {k: a * x for k, x in v.items()}
            else:
                nv = dict(v)
            for k, x in p.items():
                y = nv.get(k, 0) - b * x
                if y:
                    nv[k] = y
                else:
                    nv.pop(k, None)
            v = nv
        return False

    def reduce(self, v):
        """reduced remainder of v (exact; integer multiple of the true remainder)."""
        v = {k: x for k, x in v.items() if x}
        piv = self.piv
        while v:
            c = max(v)
            p = piv.get(c)
            if p is None:
                return v
            a = p[c]
            b = v[c]
            g = gcd(a, b)
            a //= g
            b //= g
            nv = {k: a * x for k, x in v.items()}
            for k, x in p.items():
                y = nv.get(k, 0) - b * x
                if y:
                    nv[k] = y
                else:
                    nv.pop(k, None)
            v = nv
        return v

    @property
    def rank(self):
        return len(self.piv)


class EchelonModP:
    def __init__(self, p):
        self.p = p
        self.piv = {}

    def insert(self, v):
        p = self.p
        v = {k: x % p for k, x in v.items() if x % p}
        piv = self.piv
        while v:
            c = max(v)
            row = piv.get(c)
            if row is None:
                inv = pow(v[c], p - 2, p)
                piv[c] = {k: (x * inv) % p for k, x in v.items()}
                return True
            b = v[c]
            for k, x in row.items():
                y = (v.get(k, 0) - b * x) % p
                if y:
                    v[k] = y
                else:
                    v.pop(k, None)
        return False

    @property
    def rank(self):
        return len(self.piv)


def vec(p, index):
    out = {}
    for m, c in p.items():
        out[index[m]] = c
    return out
