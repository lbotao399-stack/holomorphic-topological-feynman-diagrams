"""Independent cochain library for the o3 verification (does not import any explorer / deform code).

Cochains on g[[z1,z2]], g = gl_N or sl_N, kappa = tr, formal N (coefficients Laurent polynomials in N).
k-cochain = dict {(perm, mono): LP}
  perm: tuple of length k, perm[i] = successor of slot i in its trace cycle; perm[i] = i <-> tr(X_i).
        tau_perm(X_0..X_{k-1}) = prod_cycles tr(X_{i} X_{perm(i)} X_{perm(perm(i))} ...).
  mono: exponents (a_0..a_{k-2}, b_0..b_{k-2}) of p_0^1..p_{k-2}^1, p_0^2..p_{k-2}^2;  p_{k-1} = -sum_{i<k-1} p_i.
sl mode: perms with a fixed point are dropped (traceless arguments).
Conventions (definitions, as in the task statement):
  (sigma.f)(x_0..x_{k-1}) = f(x_sigma(0)..x_sigma(k-1)),  alt f = sum sgn(sigma) sigma.f
  d w(x_0..x_k) = sum_{i<j} (-1)^{i+j} w([x_i,x_j], x_0..^i..^j..x_k)
  {a,b} = alt(gamma_0)/((m-1)!(n-1)!),  gamma_0 = a(x_0..x_{m-2}, X_e) kappa^{ee'} b(X_e', x_{m-1}..x_{m+n-3})
  gl Casimir: tr(A X_e) tr(B X^e) = tr(AB), tr(A X_e B X^e) = trA trB;  sl: minus (1/N) trA trB, (1/N) tr(AB).
"""
from fractions import Fraction as Fr
import itertools, math


# ------------------------------------------------------------------ Laurent polynomials in N
class LP:
    __slots__ = ("c",)

    def __init__(self, c=None):
        self.c = {} if c is None else {e: v for e, v in c.items() if v != 0}

    @staticmethod
    def const(x):
        return LP({0: Fr(x)})

    @staticmethod
    def N(k=1):
        return LP({k: Fr(1)})

    def __add__(self, o):
        if not isinstance(o, LP):
            o = LP.const(o)
        r = dict(self.c)
        for e, v in o.c.items():
            r[e] = r.get(e, 0) + v
        return LP(r)
    __radd__ = __add__

    def __neg__(self):
        return LP({e: -v for e, v in self.c.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, LP) else LP.const(-Fr(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, LP):
            o = Fr(o)
            if o == 0:
                return LP()
            return LP({e: v * o for e, v in self.c.items()})
        r = {}
        for e1, v1 in self.c.items():
            for e2, v2 in o.c.items():
                r[e1 + e2] = r.get(e1 + e2, 0) + v1 * v2
        return LP(r)
    __rmul__ = __mul__

    def __truediv__(self, o):
        return self * (Fr(1) / Fr(o))

    def __bool__(self):
        return bool(self.c)

    def __eq__(self, o):
        if not isinstance(o, LP):
            o = LP.const(o)
        return self.c == o.c

    def __hash__(self):
        return hash(tuple(sorted(self.c.items())))

    def at(self, N):
        N = Fr(N)
        return sum((v * N ** e for e, v in self.c.items()), Fr(0))

    def __repr__(self):
        if not self.c:
            return "0"
        return " + ".join("(%s)N^%d" % (v, e) for e, v in sorted(self.c.items()))


ONE = LP.const(1)


# ------------------------------------------------------------------ permutations / traces
def perm_sign(p):
    p = list(p)
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def seq_sign(seq):
    """sign of the permutation taking sorted(seq) to seq (entries distinct)."""
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def cycles_of(perm):
    k = len(perm)
    seen, out = [False] * k, []
    for i in range(k):
        if not seen[i]:
            c, j = [], i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = perm[j]
            out.append(c)
    return out


def perm_from_words(words, K):
    """words: list of cyclic words over target slots 0..K-1 (each slot exactly once) -> perm."""
    nxt = [None] * K
    for w in words:
        for a, b in zip(w, w[1:] + w[:1]):
            nxt[a] = b
    assert None not in nxt, (words, K)
    return tuple(nxt)


def has_fixed(perm):
    return any(perm[i] == i for i in range(len(perm)))


def perm_str(perm):
    return "".join("tr(" + "".join(str(i + 1) for i in c) + ")" for c in cycles_of(perm))


# ------------------------------------------------------------------ momentum polynomials
_powcache = {}


def _lin_pow(r, e):
    """(sum_j r_j u_j)^e as dict {exp tuple: Fr}."""
    key = (r, e)
    if key in _powcache:
        return _powcache[key]
    n = len(r)
    if e == 0:
        res = {tuple([0] * n): Fr(1)}
    else:
        prev = _lin_pow(r, e - 1)
        res = {}
        for m, c in prev.items():
            for j in range(n):
                if r[j]:
                    mm = list(m)
                    mm[j] += 1
                    mm = tuple(mm)
                    res[mm] = res.get(mm, 0) + c * r[j]
        res = {m: c for m, c in res.items() if c}
    _powcache[key] = res
    return res


def _pmul(A, B):
    r = {}
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            r[m] = r.get(m, 0) + c1 * c2
    return {m: c for m, c in r.items() if c}


_subcache = {}


def subst_mono(mono, k, rows, K):
    """mono of a k-slot cochain; rows[s] (s < k-1): linear form of slot s's momentum over the FULL target
    momenta (length K).  Returns dict {target mono (length 2(K-1)): Fr}."""
    key = (mono, k, rows, K)
    if key in _subcache:
        return _subcache[key]
    red = [tuple(row[j] - row[K - 1] for j in range(K - 1)) for row in rows]
    a, b = mono[:k - 1], mono[k - 1:]
    P1 = {tuple([0] * (K - 1)): Fr(1)}
    P2 = {tuple([0] * (K - 1)): Fr(1)}
    for s in range(k - 1):
        if a[s]:
            P1 = _pmul(P1, _lin_pow(red[s], a[s]))
        if b[s]:
            P2 = _pmul(P2, _lin_pow(red[s], b[s]))
    res = {}
    for m1, c1 in P1.items():
        for m2, c2 in P2.items():
            res[m1 + m2] = res.get(m1 + m2, 0) + c1 * c2
    res = {m: c for m, c in res.items() if c}
    _subcache[key] = res
    return res


def unit(i, K):
    return tuple(1 if j == i else 0 for j in range(K))


# ------------------------------------------------------------------ cochain basics
def add_to(out, key, c):
    v = out.get(key)
    v = c if v is None else v + c
    if v:
        out[key] = v
    elif key in out:
        del out[key]


def clean(w):
    return {k: v for k, v in w.items() if v}


def lin(*pairs):
    """linear combination sum c_i w_i."""
    out = {}
    for c, w in pairs:
        for key, v in w.items():
            add_to(out, key, v * c)
    return out


def scal(c, w):
    return {key: v * c for key, v in w.items() if v * c}


def is_zero(w):
    return not clean(w)


def arity(w):
    for (p, m) in w:
        return len(p)
    return None


def sl_filter(w):
    return {key: v for key, v in w.items() if not has_fixed(key[0])}


def from_raw(k, terms, sl):
    """terms: list of (coeff, words over slots 0..k-1, poly over full momenta of k slots as
    {full mono (a_0..a_{k-1}, b_0..b_{k-1}): Fr}) -> cochain (reduced)."""
    out = {}
    for c, words, poly in terms:
        perm = perm_from_words(words, k)
        if sl and has_fixed(perm):
            continue
        for fm, pc in poly.items():
            # substitute p_{k-1} = -sum
            a, b = fm[:k], fm[k:]
            rows = tuple(unit(s, k) for s in range(k))
            # express full monomial as product over all k slots; handle slot k-1 via linear form
            red_last = tuple([-1] * (k - 1))
            P1 = {tuple([0] * (k - 1)): Fr(1)}
            P2 = {tuple([0] * (k - 1)): Fr(1)}
            for s in range(k):
                r = unit(s, k - 1) if s < k - 1 else red_last
                if a[s]:
                    P1 = _pmul(P1, _lin_pow(r, a[s]))
                if b[s]:
                    P2 = _pmul(P2, _lin_pow(r, b[s]))
            for m1, c1 in P1.items():
                for m2, c2 in P2.items():
                    add_to(out, (perm, m1 + m2), c * (c1 * c2 * pc))
    return out


def substitute(w, k, slot_words, slot_rows, K, sl):
    """General substitution: slot s of w gets the word slot_words[s] (tuple of target slots) and the momentum
    linear form slot_rows[s] (full target, length K).  Returns K-cochain dict."""
    out = {}
    rows = tuple(slot_rows[s] for s in range(k - 1))
    for (perm, mono), c in w.items():
        words = []
        for cyc in cycles_of(perm):
            word = []
            for s in cyc:
                word.extend(slot_words[s])
            words.append(word)
        nperm = perm_from_words(words, K)
        if sl and has_fixed(nperm):
            continue
        for m2, c2 in subst_mono(mono, k, rows, K).items():
            add_to(out, (nperm, m2), c * c2)
    return out


def act(w, sigma, sl):
    """(sigma.w)(x_0..) = w(x_sigma(0), ..)."""
    k = len(sigma)
    return substitute(w, k, [(sigma[s],) for s in range(k)], [unit(sigma[s], k) for s in range(k)], k, sl)


def alt(w, k, sl):
    out = {}
    for sigma in itertools.permutations(range(k)):
        sg = perm_sign(sigma)
        for key, v in act(w, sigma, sl).items():
            add_to(out, key, v * sg)
    return out


def d(w, k, sl):
    """CE differential of a k-cochain (k >= 1)."""
    K = k + 1
    out = {}
    for i in range(K):
        for j in range(i + 1, K):
            sg = (-1) ** (i + j)
            rest = [l for l in range(K) if l not in (i, j)]
            rows = [tuple(1 if (t == i or t == j) else 0 for t in range(K))] + [unit(r, K) for r in rest]
            for word0, s0 in (((i, j), 1), ((j, i), -1)):
                words = [word0] + [(r,) for r in rest]
                part = substitute(w, k, words, rows, K, sl)
                for key, v in part.items():
                    add_to(out, key, v * (sg * s0))
    return out


def bracket(a, m, b, n, sl):
    """{a,b}, a of arity m, b of arity n (both alternating)."""
    K = m + n - 2
    g0 = {}
    a_rows = tuple(unit(t, K) for t in range(m - 1))
    qe_p = tuple(1 if t < m - 1 else 0 for t in range(K))           # momentum of X_e' = sum of a's args
    b_rows = (qe_p,) + tuple(unit(m - 1 + t, K) for t in range(n - 2))
    amap = {s: s for s in range(m - 1)}
    bmap = {s: m - 2 + s for s in range(1, n)}
    Ninv = LP.N(-1)
    NN = LP.N(1)
    for (pa, ma), ca in a.items():
        cyc_a = cycles_of(pa)
        A, other_a = None, []
        for c in cyc_a:
            if (m - 1) in c:
                i = c.index(m - 1)
                rot = c[i + 1:] + c[:i]
                A = [amap[s] for s in rot]
            else:
                other_a.append([amap[s] for s in c])
        Pa = subst_mono(ma, m, a_rows, K)
        for (pb, mb), cb in b.items():
            B, other_b = None, []
            for c in cycles_of(pb):
                if 0 in c:
                    i = c.index(0)
                    rot = c[i + 1:] + c[:i]
                    B = [bmap[s] for s in rot]
                else:
                    other_b.append([bmap[s] for s in c])
            Pb = subst_mono(mb, n, b_rows, K)
            P = _pmul(Pa, Pb)
            base = other_a + other_b
            terms = []
            # gl part
            if A or B:
                terms.append((ONE, base + [A + B]))
            else:
                terms.append((NN, base))
            if sl:
                fac = ONE
                ws = list(base)
                for part in (A, B):
                    if part:
                        ws.append(part)
                    else:
                        fac = fac * NN
                terms.append((-(Ninv * fac), ws))
            for cf, ws in terms:
                perm = perm_from_words(ws, K)
                if sl and has_fixed(perm):
                    continue
                cc = ca * cb * cf
                for mono, pc in P.items():
                    add_to(g0, (perm, mono), cc * pc)
    res = alt(g0, K, sl)
    return scal(Fr(1, math.factorial(m - 1) * math.factorial(n - 1)), res)


def Cop(w, k):
    """C: X -> -X^T:  tau_perm -> (-1)^k tau_{perm^{-1}}."""
    out = {}
    s = (-1) ** k
    for (perm, mono), v in w.items():
        inv = [0] * k
        for i in range(k):
            inv[perm[i]] = i
        add_to(out, (tuple(inv), mono), v * s)
    return out


def at_N(w, N):
    return {key: (v.at(N) if isinstance(v, LP) else Fr(v)) for key, v in w.items() if (v.at(N) if isinstance(v, LP) else v) != 0}


def sl2_E(w, k):
    """E = sum_j p_j^1 d/dp_j^2 applied to the momentum polynomials (reduced variables)."""
    out = {}
    n = k - 1
    for (perm, mono), v in w.items():
        a, b = list(mono[:n]), list(mono[n:])
        for j in range(n):
            if b[j]:
                aa, bb = list(a), list(b)
                c = bb[j]
                bb[j] -= 1
                aa[j] += 1
                add_to(out, (perm, tuple(aa) + tuple(bb)), v * c)
    return out


def sl2_F(w, k):
    out = {}
    n = k - 1
    for (perm, mono), v in w.items():
        a, b = list(mono[:n]), list(mono[n:])
        for j in range(n):
            if a[j]:
                aa, bb = list(a), list(b)
                c = aa[j]
                aa[j] -= 1
                bb[j] += 1
                add_to(out, (perm, tuple(aa) + tuple(bb)), v * c)
    return out


# ------------------------------------------------------------------ standard cochains
def bracket_poly(i, j, k):
    """[p_i p_j] as full-momentum polynomial on k slots."""
    def m(s1, s2):
        e = [0] * (2 * k)
        e[s1] += 1
        e[k + s2] += 1
        return tuple(e)
    out = {}
    add_to(out, m(i, j), Fr(1))
    add_to(out, m(j, i), Fr(-1))
    return out


def fpmul(A, B):
    return _pmul(A, B)


def fpow(A, e, k):
    r = {tuple([0] * (2 * k)): Fr(1)}
    for _ in range(e):
        r = _pmul(r, A)
    return r


TENS3 = {
    'f': [(1, [[0, 1, 2]]), (-1, [[0, 2, 1]])],
    'd': [(1, [[0, 1, 2]]), (1, [[0, 2, 1]])],
    's1': [(1, [[0], [1, 2]]), (1, [[1], [0, 2]]), (1, [[2], [0, 1]])],
    's0': [(1, [[0], [1], [2]])],
}


def cubic(tname, w, sl):
    """[p_0 p_1]^w * T(X_0,X_1,X_2) as a reduced 3-cochain (alternating on sum p = 0)."""
    poly = fpow(bracket_poly(0, 1, 3), w, 3)
    terms = [(LP.const(c), words, poly) for c, words in TENS3[tname]]
    return from_raw(3, terms, sl)


# ------------------------------------------------------------------ spanning sets
def momentum_monos(k, deg):
    """all reduced monomials of bidegree deg in k-1 variables."""
    n = k - 1
    def comps(total, n):
        if n == 0:
            if total == 0:
                yield ()
            return
        for x in range(total + 1):
            for rest in comps(total - x, n - 1):
                yield (x,) + rest
    return [a + b for a in comps(deg[0], n) for b in comps(deg[1], n)]


def cycle_type_reps(k, sl):
    """one perm per cycle type (partition of k); sl: no parts equal to 1."""
    def parts(n, mx):
        if n == 0:
            yield ()
            return
        for p in range(min(n, mx), 0, -1):
            for rest in parts(n - p, p):
                yield (p,) + rest
    out = []
    for pt in parts(k, k):
        if sl and 1 in pt:
            continue
        words, s = [], 0
        for p in pt:
            words.append(list(range(s, s + p)))
            s += p
        out.append(perm_from_words(words, k))
    return out


def all_perms(k, sl):
    return [p for p in itertools.permutations(range(k)) if not (sl and has_fixed(p))]


def span_set(k, deg, sl, sl2inv=False, parity=None):
    """spanning set of C^k_deg: alt(tau_rep * m).  sl2inv: m over products of brackets (deg must be (w,w)).
    parity: None, +1, -1 (C-projection (w +- Cw)/2)."""
    out = []
    reps = cycle_type_reps(k, sl)
    if sl2inv:
        assert deg[0] == deg[1]
        pairs = list(itertools.combinations(range(k), 2))
        polys = []
        for tup in itertools.combinations_with_replacement(pairs, deg[0]):
            P = {tuple([0] * (2 * k)): Fr(1)}
            for (i, j) in tup:
                P = _pmul(P, bracket_poly(i, j, k))
            polys.append(P)
        reps = all_perms(k, sl)     # brackets are not stable under relabelling: use all taus
        for tau in reps:
            for P in polys:
                w = from_raw(k, [(ONE, [c for c in cycles_of(tau)], P)], sl)
                a = alt(w, k, sl)
                if parity is not None:
                    a = scal(Fr(1, 2), lin((1, a), (parity, Cop(a, k))))
                if a:
                    out.append(a)
        return out
    for tau in reps:
        for mono in momentum_monos(k, deg):
            w = {(tau, mono): ONE}
            a = alt(w, k, sl)
            if parity is not None:
                a = scal(Fr(1, 2), lin((1, a), (parity, Cop(a, k))))
            if a:
                out.append(a)
    return out
