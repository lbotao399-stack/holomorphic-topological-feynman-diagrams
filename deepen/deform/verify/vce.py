"""Independent CE engine (verifier).  Not importing explorer code.

Cochain of arity k: dict  key=(sigma, mono, npow) -> Fraction
  sigma: permutation tuple of {0..k-1};  tau_sigma(X_0..X_{k-1}) = prod over cycles (c, s(c), s^2(c), ..)
         of tr(X_c X_{s(c)} X_{s^2(c)} ...)
  mono : exponent tuple of length 2k, index 2*slot + comp (comp 0 = first component, 1 = second);
         reduced form: exponents of slot k-1 are zero (p_{k-1} = -(p_0+..+p_{k-2}))
  npow : integer power of N
Value: w(X_0 e^{p_0 z},..) = sum coeff * N^npow * mono(p) * tau_sigma(X)  on sum p = 0.

Conventions:
  (pi.w)(x_0..x_{k-1}) = w(x_{pi(0)}, .., x_{pi(k-1)});  alt w = sum_pi sgn(pi) pi.w
  d w(x_0..x_k) = sum_{i<j} (-1)^{i+j} w([x_i,x_j], x_0..^i..^j..x_k)
  {a,b} = alt(gamma_0)/((m-1)!(n-1)!),
      gamma_0(x_1..x_{m-1},y_1..y_{n-1}) = a(x.., X_e) kappa^{ee'} b(X_e', y..)   (kappa = tr)
"""
import itertools
from fractions import Fraction as Fr
from math import factorial, comb
from functools import lru_cache


# ---------------------------------------------------------------- permutations
def perm_sign(p):
    s = 1
    seen = [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def perm_cycles(p):
    seen = [False] * len(p)
    out = []
    for i in range(len(p)):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = p[j]
            out.append(tuple(c))
    return out


def perm_from_cycles(cyc, k):
    p = [None] * k
    for c in cyc:
        for a in range(len(c)):
            p[c[a]] = c[(a + 1) % len(c)]
    assert None not in p, (cyc, k)
    return tuple(p)


def has_fixed(p):
    return any(p[i] == i for i in range(len(p)))


ALLPERMS = {k: list(itertools.permutations(range(k))) for k in range(0, 7)}


# ---------------------------------------------------------------- polynomials in momenta
def _pow_linear(coefs, n, nvars):
    """(sum_s coefs[s] x_s)^n as dict exponent-tuple(len nvars)->int."""
    items = [(s, c) for s, c in coefs.items() if c != 0]
    out = {}
    if n == 0:
        return {tuple([0] * nvars): 1}
    if not items:
        return {}
    idx = [s for s, _ in items]
    for alpha in _compositions(n, len(items)):
        c = factorial(n)
        e = [0] * nvars
        for (s, a), al in zip(items, alpha):
            c //= factorial(al)
            c *= a ** al
            e[s] += al
        out[tuple(e)] = out.get(tuple(e), 0) + c
    return out


def _compositions(n, m):
    if m == 1:
        yield (n,)
        return
    for a in range(n + 1):
        for rest in _compositions(n - a, m - 1):
            yield (a,) + rest


def _polymul(A, B):
    out = {}
    for ea, ca in A.items():
        for eb, cb in B.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c != 0}


_subst_cache = {}


def subst(mono, M, knew):
    """mono: exponents over old slots (len 2*kold).  M: tuple over old slots of tuple((newslot, coef),..).
    New momentum of old slot b = sum coef * p_newslot.  The last new slot is eliminated.
    Returns dict mono(len 2*knew) -> int."""
    key = (mono, M, knew)
    r = _subst_cache.get(key)
    if r is not None:
        return r
    last = knew - 1
    res = {tuple([0] * (2 * knew)): 1}
    for b, row in enumerate(M):
        lin = {}
        for s, c in row:
            if s == last:
                for t in range(last):
                    lin[t] = lin.get(t, 0) - c
            else:
                lin[s] = lin.get(s, 0) + c
        for comp in (0, 1):
            e = mono[2 * b + comp]
            if e == 0:
                continue
            # power of the linear form in component comp
            pw = _pow_linear(lin, e, knew)  # exponents per slot
            P = {}
            for es, c in pw.items():
                ee = [0] * (2 * knew)
                for s in range(knew):
                    ee[2 * s + comp] = es[s]
                P[tuple(ee)] = c
            res = _polymul(res, P)
    _subst_cache[key] = res
    return res


def monos(k, a, b):
    """reduced monomials of bidegree (a,b) for arity k (slots 0..k-2)."""
    nv = k - 1
    out = []
    for ea in _compositions(a, nv) if nv > 0 else ([()] if a == 0 else []):
        for eb in _compositions(b, nv) if nv > 0 else ([()] if b == 0 else []):
            e = [0] * (2 * k)
            for s in range(nv):
                e[2 * s] = ea[s]
                e[2 * s + 1] = eb[s]
            out.append(tuple(e))
    return out


def bideg(mono):
    return (sum(mono[0::2]), sum(mono[1::2]))


# ---------------------------------------------------------------- cochains
def add(out, key, c):
    v = out.get(key, 0) + c
    if v == 0:
        out.pop(key, None)
    else:
        out[key] = v


def lin(*pairs):
    out = {}
    for s, w in pairs:
        for k, c in w.items():
            add(out, k, s * c)
    return out


def scale(w, s):
    return {k: s * c for k, c in w.items() if s * c != 0}


def arity(w):
    for (sig, m, n) in w:
        return len(sig)
    return None


def act(pi, w):
    """(pi.w)(x_0..) = w(x_pi(0),..)"""
    k = len(pi)
    inv = [0] * k
    for i, x in enumerate(pi):
        inv[x] = i
    M = tuple(((pi[b], 1),) for b in range(k))
    out = {}
    for (sig, m, n), c in w.items():
        sig2 = [0] * k
        for cc in range(k):
            sig2[pi[cc]] = pi[sig[cc]]
        sig2 = tuple(sig2)
        for m2, c2 in subst(m, M, k).items():
            add(out, (sig2, m2, n), c * c2)
    return out


def alt(w, k=None):
    if k is None:
        k = arity(w)
    out = {}
    for pi in ALLPERMS[k]:
        s = perm_sign(pi)
        for key, c in act(pi, w).items():
            add(out, key, s * c)
    return out


def d(w, sl=False):
    k = arity(w)
    if k is None:
        return {}
    out = {}
    K = k + 1
    for i in range(K):
        for j in range(i + 1, K):
            sgn = (-1) ** (i + j)
            rest = [t for t in range(K) if t != i and t != j]
            mp = [None] + rest  # w-slot b>=1 -> rest[b-1]
            M = (((i, 1), (j, 1)),) + tuple(((rest[b - 1], 1),) for b in range(1, k))
            for (sig, m, n), c in w.items():
                pol = subst(m, M, K)
                for (I, J, s2) in ((i, j, 1), (j, i, -1)):
                    sp = [None] * K
                    for b in range(1, k):
                        if sig[b] != 0:
                            sp[mp[b]] = mp[sig[b]]
                    u = sig.index(0)
                    if u != 0:
                        sp[mp[u]] = I
                        sp[I] = J
                        sp[J] = mp[sig[0]]
                    else:
                        sp[I] = J
                        sp[J] = I
                    sp = tuple(sp)
                    if sl and has_fixed(sp):
                        continue
                    for m2, c2 in pol.items():
                        add(out, (sp, m2, n), sgn * s2 * c * c2)
    return out


def _cycle_with(sig, x):
    """word after x in its cycle (excluding x), and the other cycles."""
    cyc = perm_cycles(sig)
    word = None
    others = []
    for c in cyc:
        if x in c:
            a = c.index(x)
            word = list(c[a + 1:] + c[:a])
        else:
            others.append(list(c))
    return word, others


def gamma0(a, b, sl=False):
    m = arity(a)
    n = arity(b)
    r = m + n - 2
    out = {}
    # b slot 0 -> contracted; b slot c (>=1) -> result slot m-2+c
    Mb = (tuple((s, 1) for s in range(m - 1)),) + tuple(((m - 2 + c, 1),) for c in range(1, n))
    Ma = tuple(((s, 1),) for s in range(m - 1)) + (tuple((s, -1) for s in range(m - 1)),)
    # a polynomial: slot m-1 exponent zero (reduced); map identically
    for (sa, ma, na), ca in a.items():
        Pw, Oa = _cycle_with(sa, m - 1)
        pa = subst(ma, Ma, r)
        for (sb, mb, nb), cb in b.items():
            Qw, Ob = _cycle_with(sb, 0)
            Qw = [m - 2 + c for c in Qw]
            Ob = [[m - 2 + c for c in cyc] for cyc in Ob]
            pb = subst(mb, Mb, r)
            pol = _polymul(pa, pb)
            base = Oa + Ob
            terms = []
            merged = Pw + Qw
            if merged:
                terms.append((base + [merged], 0, Fr(1)))
            else:
                terms.append((base, 1, Fr(1)))
            if sl:
                assert Pw and Qw
                terms.append((base + [Pw, Qw], -1, Fr(-1)))
            for cyc, dn, cf in terms:
                sp = perm_from_cycles(cyc, r)
                if sl and has_fixed(sp):
                    continue
                for mm, cp in pol.items():
                    add(out, (sp, mm, na + nb + dn), ca * cb * cf * cp)
    return out


def bracket(a, b, sl=False):
    m = arity(a)
    n = arity(b)
    g = gamma0(a, b, sl)
    return scale(alt(g, m + n - 2), Fr(1, factorial(m - 1) * factorial(n - 1)))


def evalN(w, N):
    """substitute numeric N -> dict (sigma, mono) -> Fraction"""
    out = {}
    for (s, m, n), c in w.items():
        add(out, (s, m), c * Fr(N) ** n)
    return out


def npieces(w):
    out = {}
    for (s, m, n), c in w.items():
        out.setdefault(n, {})[(s, m, 0)] = c
    return out


def restrict_sl(w):
    return {k: c for k, c in w.items() if not has_fixed(k[0])}


# ---------------------------------------------------------------- constructors
def bracket_pq(k, a, b):
    """[p_a p_b] = p_a^1 p_b^2 - p_a^2 p_b^1 as reduced poly for arity k."""
    M = tuple(((s, 1),) for s in range(k))
    e1 = [0] * (2 * k); e1[2 * a] += 1; e1[2 * b + 1] += 1
    e2 = [0] * (2 * k); e2[2 * a + 1] += 1; e2[2 * b] += 1
    out = {}
    for mm, c in subst(tuple(e1), M, k).items():
        out[mm] = out.get(mm, 0) + c
    for mm, c in subst(tuple(e2), M, k).items():
        out[mm] = out.get(mm, 0) - c
    return {mm: c for mm, c in out.items() if c != 0}


def polyprod(*ps):
    r = ps[0]
    for p in ps[1:]:
        r = _polymul(r, p)
    return r


def from_poly_tau(poly, taus):
    """poly: dict mono->coef; taus: list of (coef, sigma)."""
    out = {}
    for mm, c in poly.items():
        for ct, sig in taus:
            add(out, (tuple(sig), mm, 0), Fr(c) * ct)
    return out


def cyc_perm(k, *cycles):
    return perm_from_cycles([tuple(c) for c in cycles], k)
