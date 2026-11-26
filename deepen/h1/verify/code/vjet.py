"""Verifier's own jet model (independent of deepen/h1/code).

Directions: u = 0 tau, 1 zb1, 2 zb2, 3 z1, 4 z2.  Form directions v in {0,1,2}; T is a bitmask (bit v <-> dx_v),
dx_T ordered 0 < 1 < 2.  Generator (a, T, al): jet coordinate c^a_{T,al} of the superfield
A^a = sum_{T,al} c^a_{T,al} x^al/al! dx_T, parity (1-|T|) mod 2, ghost 1-|T|.
Superalgebra R = C[c] (x) Lambda[dx] (x) C[[x]], (c1 y1)(c2 y2) = (-1)^{|y1||c2|} c1 c2 y1 y2.
Q c^a_{T,b} := coefficient of x^b/b! dx_T in  2 dA^a + f^a_{bc} A^b A^c   (= 2 x [dA + 1/2 [A,A]]);
  the factor 2 does not change kernels or images.
D_u c_{T,al} = c_{T,al+e_u} (even derivation).
eta_v c_{T,al} := coefficient of iota_v A  (odd derivation, iota_v = 1 (x) iota_v on R).
Arithmetic: integer coefficients; ranks mod primes and exactly over Q (fractions.Fraction).
"""
import sys
from fractions import Fraction
from math import comb

P1 = 1000000007
P2 = 998244353


def popc(T):
    return bin(T).count("1")


def eps_vT(v, T):
    """dx_v dx_T = eps * dx_{T|v} (v not in T)."""
    return -1 if popc(T & ((1 << v) - 1)) % 2 else 1


def shuffle_sign(T1, T2):
    """dx_T1 dx_T2 = sign * dx_{T1|T2} (disjoint)."""
    s = 0
    for i in range(3):
        if (T1 >> i) & 1:
            s += popc(T2 & ((1 << i) - 1))
    return -1 if s % 2 else 1


def submasks(T):
    S = T
    while True:
        yield S
        if S == 0:
            return
        S = (S - 1) & T


def sub_multi(be):
    """all al <= be componentwise."""
    out = [()]
    for b in be:
        out = [o + (k,) for o in out for k in range(b + 1)]
    return out


def lie(name):
    """returns dim, f (dict (a,b)-> {c: coef}, [t_a,t_b] = sum_c f^c_ab t_c), charge (eigenvalue of ad h on t_a or None)."""
    f = {}

    def setb(a, b, d):
        f[(a, b)] = dict(d)
        f[(b, a)] = {c: -v for c, v in d.items()}
    if name == "ab1":
        return 1, f, None
    if name == "ab2":
        return 2, f, None
    if name in ("sl2", "gl2"):
        e, ff, h = 0, 1, 2
        setb(e, ff, {h: 1})
        setb(h, e, {e: 2})
        setb(h, ff, {ff: -2})
        if name == "sl2":
            return 3, f, [2, -2, 0]
        return 4, f, [2, -2, 0, 0]
    if name == "aff1":       # [x, y] = y  (non-unimodular, no invariant form)
        setb(0, 1, {1: 1})
        return 2, f, [0, 1]  # ad x eigenvalues (diagonalizable)
    if name == "heis":       # [x, y] = z
        setb(0, 1, {2: 1})
        return 3, f, None
    raise ValueError(name)


class Jets:
    def __init__(self, name):
        self.name = name
        self.dim, self.f, self.charge = lie(name)
        self.gid = {}
        self.ginfo = []
        self.par = []
        self.cost = []
        self.gh = []
        self.chg = []
        self._Q = {}

    def g(self, a, T, al):
        key = (a, T, al)
        i = self.gid.get(key)
        if i is None:
            i = len(self.ginfo)
            self.gid[key] = i
            self.ginfo.append(key)
            self.par.append((1 - popc(T)) % 2)
            self.gh.append(1 - popc(T))
            self.cost.append((al[0] + (T & 1), al[1] + ((T >> 1) & 1), al[2] + ((T >> 2) & 1), al[3], al[4]))
            self.chg.append(-self.charge[a] if self.charge else 0)
        return i

    # ---------------------------------------------------------------- monomials
    def norm(self, seq):
        arr = list(seq)
        par = self.par
        sign = 1
        for i in range(1, len(arr)):
            x = arr[i]
            px = par[x]
            j = i
            while j > 0 and arr[j - 1] > x:
                if px and par[arr[j - 1]]:
                    sign = -sign
                arr[j] = arr[j - 1]
                j -= 1
            arr[j] = x
        for i in range(len(arr) - 1):
            if arr[i] == arr[i + 1] and par[arr[i]]:
                return 0, None
        return sign, tuple(arr)

    def mono_par(self, m):
        return sum(self.par[x] for x in m) % 2

    def ghost(self, m):
        return sum(self.gh[x] for x in m)

    def pmul(self, P, Qp):
        R = {}
        for m1, c1 in P.items():
            for m2, c2 in Qp.items():
                s, m = self.norm(m1 + m2)
                if s:
                    add(R, m, s * c1 * c2)
        return R

    # ---------------------------------------------------------------- Q on generators
    def Qgen(self, i):
        r = self._Q.get(i)
        if r is not None:
            return r
        a, T, be = self.ginfo[i]
        res = {}
        for v in range(3):
            if (T >> v) & 1:
                Tp = T & ~(1 << v)
                s = (-1 if (1 - popc(Tp)) % 2 else 1) * eps_vT(v, Tp)
                al = list(be)
                al[v] += 1
                add(res, (self.g(a, Tp, tuple(al)),), 2 * s)
        for (b, c), fd in self.f.items():
            coef = fd.get(a, 0)
            if not coef:
                continue
            for T1 in submasks(T):
                T2 = T ^ T1
                s = shuffle_sign(T1, T2) * (-1 if (popc(T1) * ((1 - popc(T2)) % 2)) % 2 else 1)
                for al1 in sub_multi(be):
                    al2 = tuple(x - y for x, y in zip(be, al1))
                    m = 1
                    for x, y in zip(be, al1):
                        m *= comb(x, y)
                    g1 = self.g(b, T1, al1)
                    g2 = self.g(c, T2, al2)
                    sg, mono = self.norm([g1, g2])
                    if sg:
                        add(res, mono, coef * s * m * sg)
        self._Q[i] = res
        return res

    def Qmono(self, m):
        res = {}
        pre = 1
        for i, x in enumerate(m):
            qx = self.Qgen(x)
            left = m[:i]
            right = m[i + 1:]
            for mono, c in qx.items():
                sg, nm = self.norm(left + mono + right)
                if sg:
                    add(res, nm, pre * sg * c)
            if self.par[x]:
                pre = -pre
        return res

    def Qpoly(self, P):
        R = {}
        for m, c in P.items():
            for mm, cc in self.Qmono(m).items():
                add(R, mm, c * cc)
        return R

    def Dmono(self, m, u):
        res = {}
        for i, x in enumerate(m):
            a, T, al = self.ginfo[x]
            al2 = list(al)
            al2[u] += 1
            y = self.g(a, T, tuple(al2))
            sg, nm = self.norm(m[:i] + (y,) + m[i + 1:])
            if sg:
                add(res, nm, sg)
        return res

    def Dpoly(self, P, u):
        R = {}
        for m, c in P.items():
            for mm, cc in self.Dmono(m, u).items():
                add(R, mm, c * cc)
        return R

    def etagen(self, x, v):
        a, T, al = self.ginfo[x]
        if (T >> v) & 1:
            return None
        Tn = T | (1 << v)
        s = (-1 if (1 - popc(Tn)) % 2 else 1) * eps_vT(v, T)
        return s, self.g(a, Tn, al)

    def etamono(self, m, v):
        res = {}
        pre = 1
        for i, x in enumerate(m):
            r = self.etagen(x, v)
            if r is not None:
                s, y = r
                sg, nm = self.norm(m[:i] + (y,) + m[i + 1:])
                if sg:
                    add(res, nm, pre * s * sg)
            if self.par[x]:
                pre = -pre
        return res

    def etapoly(self, P, v):
        R = {}
        for m, c in P.items():
            for mm, cc in self.etamono(m, v).items():
                add(R, mm, c * cc)
        return R

    # ---------------------------------------------------------------- slices
    def gens_upto(self, t):
        out = []
        for a in range(self.dim):
            for T in range(8):
                lim = [t[v] - ((T >> v) & 1) for v in range(3)] + [t[3], t[4]]
                if min(lim) < 0:
                    continue
                for al in sub_multi(tuple(lim)):
                    i = self.g(a, T, al)
                    if any(self.cost[i]):
                        out.append(i)
        return sorted(out)

    def slice(self, t, charge0=False):
        """dict ghost -> sorted list of monomials of total cost t (reduced)."""
        t = tuple(t)
        pos = self.gens_upto(t)
        zero = [self.g(a, 0, (0, 0, 0, 0, 0)) for a in range(self.dim)]
        cost = self.cost
        res = []

        def rec(j, rem, cur):
            if not any(rem):
                res.append(tuple(cur))
                return
            for k in range(j, len(pos)):
                x = pos[k]
                cx = cost[x]
                if all(cx[u] <= rem[u] for u in range(5)):
                    nrem = tuple(rem[u] - cx[u] for u in range(5))
                    if self.par[x]:
                        rec(k + 1, nrem, cur + [x])
                    else:
                        rec(k, nrem, cur + [x])
        rec(0, t, [])
        out = {}
        nz = len(zero)
        for base in res:
            for s in range(1 << nz):
                z = [zero[i] for i in range(nz) if (s >> i) & 1]
                seq = list(base) + z
                if not seq:
                    continue
                if charge0 and sum(self.chg[x] for x in seq) != 0:
                    continue
                sg, m = self.norm(seq)
                if not sg:
                    continue
                out.setdefault(self.ghost(m), []).append(m)
        for k in out:
            out[k] = sorted(set(out[k]))
        return out


def add(D, k, v):
    if v:
        w = D.get(k, 0) + v
        if w:
            D[k] = w
        else:
            del D[k]


# ---------------------------------------------------------------- linear algebra
class Echelon:
    """row echelon basis over F_p (p int) or Q (p = None); rows dict col->val, leading = min col, leading coeff 1."""

    def __init__(self, p=None):
        self.p = p
        self.piv = {}

    def reduce(self, v):
        p = self.p
        if p:
            v = {k: x % p for k, x in v.items() if x % p}
        else:
            v = {k: Fraction(x) for k, x in v.items() if x}
        piv = self.piv
        while v:
            c = min(v)
            r = piv.get(c)
            if r is None:
                return v, c
            a = v[c]
            for k, x in r.items():
                y = v.get(k, 0) - a * x
                if p:
                    y %= p
                if y:
                    v[k] = y
                else:
                    v.pop(k, None)
        return v, None

    def nf(self, v):
        """normal form modulo the span (all pivot columns eliminated)."""
        p = self.p
        if p:
            v = {k: x % p for k, x in v.items() if x % p}
        else:
            v = {k: Fraction(x) for k, x in v.items() if x}
        piv = self.piv
        done = {}
        while v:
            c = min(v)
            r = piv.get(c)
            a = v.pop(c)
            if r is None:
                done[c] = a
                continue
            for k, x in r.items():
                if k == c:
                    continue
                y = v.get(k, 0) - a * x
                if p:
                    y %= p
                if y:
                    v[k] = y
                else:
                    v.pop(k, None)
        return done

    def add(self, v):
        w, c = self.reduce(v)
        if c is None:
            return False
        p = self.p
        a = w[c]
        if p:
            inv = pow(a, p - 2, p)
            w = {k: (x * inv) % p for k, x in w.items()}
        else:
            w = {k: x / a for k, x in w.items()}
        self.piv[c] = w
        return True

    def rank(self):
        return len(self.piv)
