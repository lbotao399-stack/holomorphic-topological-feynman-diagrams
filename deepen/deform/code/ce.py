"""Chevalley-Eilenberg cochains of g[[z1,z2]] modulo translations, invariant, in the trace basis.

Cochain of arity k (alternating, Taylor convention):  dict  (tau, mono) -> Npoly
  tau  : tuple of cycles over the slot labels 0..k-1 (cycle = tuple, rotated to start at its minimum),
         tau(X_0..X_{k-1}) = prod_c tr(prod_{i in c} X_i)
  mono : canonical monomial, exponent tuple over (slot, component) for slots 0..k-2, components 1,2
         (slot-major); the momentum of slot k-1 is eliminated by p_{k-1} = -(p_0 + .. + p_{k-2}).
  Npoly: dict {power of N (int, may be negative): Fraction}.
Value:  omega(X_0 e^{p_0 z}, .., X_{k-1} e^{p_{k-1} z}) = sum tau(X) mono(p) Npoly(N)  on sum p = 0.

Conventions (alternating cochains of even arguments x_i = X_i e^{p_i z}):
  (sigma.f)(x_1..x_k) = f(x_sigma(1), .., x_sigma(k)),  alt f = sum_sigma sgn(sigma) sigma.f
  d omega(x_0..x_k)   = sum_{i<j} (-1)^{i+j} omega([x_i,x_j], x_0, ..^i..^j.., x_k)
  {a,b}(x_1..x_{m+n-2}) = sum_{sigma in Sh(m-1,n-1)} sgn(sigma) a(x_s(1)..x_s(m-1), X_e) kappa^{ee'} b(X_e', x_s(m)..)
                        = alt(gamma_0) / ((m-1)!(n-1)!),  gamma_0 = a(x_1..x_{m-1}, X_e) kappa^{ee'} b(X_e', y_1..y_{n-1}),
  the contracted pair carries momenta q, -q fixed by conservation on each side.
  Casimir: gl_N: sum_e (X_e)_ij (X^e)_kl = d_il d_jk  ->  tr(A X_e) tr(B X^e) = tr(AB), tr(1) = N.
           sl_N: minus (1/N) d_ij d_kl            ->  extra term -(1/N) tr(A) tr(B).
  sl_N cochains are represented by their restriction to sl_N (taus with a 1-cycle dropped).
"""
import itertools
from fractions import Fraction as Fr
from math import factorial

H = 2


# ------------------------------------------------------------------ N polynomials
def np_add(a, b, s=1):
    out = dict(a)
    for p, v in b.items():
        out[p] = out.get(p, 0) + s * v
    return {p: v for p, v in out.items() if v != 0}


def np_mul(a, b):
    out = {}
    for p, v in a.items():
        for q, w in b.items():
            out[p + q] = out.get(p + q, 0) + v * w
    return {p: v for p, v in out.items() if v != 0}


def np_scale(a, c):
    return {p: v * c for p, v in a.items() if v * c != 0}


def np_eval(a, N):
    return sum(Fr(v) * Fr(N) ** p for p, v in a.items())


def add_term(out, key, npoly, c=1):
    e = out.get(key)
    if e is None:
        e = {}
        out[key] = e
    for p, v in npoly.items():
        e[p] = e.get(p, 0) + c * v
        if e[p] == 0:
            del e[p]
    if not e:
        del out[key]


def clean(c):
    return {k: {p: v for p, v in d.items() if v != 0} for k, d in c.items() if any(v != 0 for v in d.values())}


def scale(c, s):
    if s == 0:
        return {}
    return {k: np_scale(d, s) for k, d in c.items()}


def lin(*pairs):
    """sum of s_i * c_i for pairs (s_i, c_i); s_i Fraction or Npoly-dict."""
    out = {}
    for s, c in pairs:
        for k, d in c.items():
            if isinstance(s, dict):
                add_term(out, k, np_mul(s, d))
            else:
                add_term(out, k, d, s)
    return out


# ------------------------------------------------------------------ taus
def canon_cycle(c):
    c = tuple(c)
    i = c.index(min(c))
    return c[i:] + c[:i]


def canon_tau(cycles):
    return tuple(sorted(canon_cycle(c) for c in cycles if len(c) > 0))


def relabel(tau, f):
    """each label i replaced by the sequence f[i]."""
    out = []
    for c in tau:
        seq = []
        for i in c:
            seq += list(f[i])
        out.append(seq)
    return canon_tau(out)


def all_taus(k):
    out = set()
    for s in itertools.permutations(range(k)):
        seen, cyc = set(), []
        for i in range(k):
            if i in seen:
                continue
            c, j = [], i
            while j not in seen:
                seen.add(j)
                c.append(j)
                j = s[j]
            cyc.append(c)
        out.add(canon_tau(cyc))
    return sorted(out, key=lambda t: (len(t), t))


def has_singleton(tau):
    return any(len(c) == 1 for c in tau)


def tau_str(tau):
    return " ".join("tr(%s)" % "".join(str(i + 1) for i in c) for c in tau)


# ------------------------------------------------------------------ monomials
_MULTI = {}


def compositions(total, parts):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if parts == 1:
        yield (total,)
        return
    for a in range(total + 1):
        for rest in compositions(total - a, parts - 1):
            yield (a,) + rest


def multinom(total, alpha):
    r = factorial(total)
    for a in alpha:
        r //= factorial(a)
    return r


def _elim(a, nslots):
    """(-sum_{i<nslots} q_i)^a = list of (coef, alpha)."""
    key = (a, nslots)
    if key not in _MULTI:
        _MULTI[key] = [((-1) ** a * multinom(a, al), al) for al in compositions(a, nslots)]
    return _MULTI[key]


_CANON = {}


def canon_full(fm, k):
    """full exponent tuple over k slots -> dict canonical mono -> int."""
    key = (fm, k)
    if key in _CANON:
        return _CANON[key]
    a, b = fm[2 * (k - 1)], fm[2 * (k - 1) + 1]
    base = list(fm[:2 * (k - 1)])
    out = {}
    if k == 1:
        out = {(): 1} if (a, b) == (0, 0) else {}
    else:
        for c1, al in _elim(a, k - 1):
            for c2, be in _elim(b, k - 1):
                m = list(base)
                for i in range(k - 1):
                    m[2 * i] += al[i]
                    m[2 * i + 1] += be[i]
                m = tuple(m)
                out[m] = out.get(m, 0) + c1 * c2
        out = {m: v for m, v in out.items() if v != 0}
    _CANON[key] = out
    return out


def mono_degree(m):
    return (sum(m[0::2]), sum(m[1::2]))


def monos_canonical(k, deg):
    """all canonical monomials of arity k with bidegree deg."""
    nv = k - 1
    out = []
    for e1 in compositions(deg[0], nv) if nv else ([()] if deg[0] == 0 else []):
        for e2 in compositions(deg[1], nv) if nv else ([()] if deg[1] == 0 else []):
            m = []
            for i in range(nv):
                m += [e1[i], e2[i]]
            out.append(tuple(m))
    return sorted(out)


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


_PERMS = {}


def perms(k):
    if k not in _PERMS:
        _PERMS[k] = [(s, perm_sign(s)) for s in itertools.permutations(range(k))]
    return _PERMS[k]


# ------------------------------------------------------------------ alternation
def alt_raw(raw, k, sl=False):
    """raw: dict (tau, fullmono over k slots) -> Npoly;  returns canonical alt(raw)."""
    out = {}
    for (tau, fm), npl in raw.items():
        for s, sg in perms(k):
            t2 = relabel(tau, {i: (s[i],) for i in range(k)})
            if sl and has_singleton(t2):
                continue
            m2 = [0] * (2 * k)
            for i in range(k):
                m2[2 * s[i]] = fm[2 * i]
                m2[2 * s[i] + 1] = fm[2 * i + 1]
            for cm, c in canon_full(tuple(m2), k).items():
                add_term(out, (t2, cm), npl, sg * c)
    return out


def canon_raw(raw, k, sl=False):
    """canonical form of a raw function (no alternation)."""
    out = {}
    for (tau, fm), npl in raw.items():
        if sl and has_singleton(tau):
            continue
        for cm, c in canon_full(fm, k).items():
            add_term(out, (tau, cm), npl, c)
    return out


def full_of_canon(m):
    return tuple(m) + (0, 0)


def as_raw(c):
    """canonical cochain viewed as a raw function (last slot exponent 0)."""
    return {(t, full_of_canon(m)): d for (t, m), d in c.items()}


def is_alternating(c, k, sl=False):
    """checks sigma.c = sgn(sigma) c for the generating transpositions."""
    raw = as_raw(c)
    for i in range(k - 1):
        s = list(range(k))
        s[i], s[i + 1] = s[i + 1], s[i]
        tr = {}
        for (tau, fm), npl in raw.items():
            t2 = relabel(tau, {j: (s[j],) for j in range(k)})
            m2 = [0] * (2 * k)
            for j in range(k):
                m2[2 * s[j]] = fm[2 * j]
                m2[2 * s[j] + 1] = fm[2 * j + 1]
            for cm, cc in canon_full(tuple(m2), k).items():
                add_term(tr, (t2, cm), npl, cc)
        if clean(lin((1, tr), (1, c))) != {}:
            return False
    return True


# ------------------------------------------------------------------ differential
_DEXP = {}


def _binom_split(a):
    """(q_i + q_j)^a = sum_r C(a,r) q_i^r q_j^(a-r)."""
    if a not in _DEXP:
        _DEXP[a] = [(factorial(a) // (factorial(r) * factorial(a - r)), r, a - r) for r in range(a + 1)]
    return _DEXP[a]


def d_ce(omega, k, sl=False):
    """CE differential of a canonical k-cochain -> canonical (k+1)-cochain."""
    out = {}
    K = k + 1
    for i, j in itertools.combinations(range(K), 2):
        rest = [r for r in range(K) if r not in (i, j)]
        sgn = (-1) ** (i + j)
        f1 = {0: (i, j)}
        f2 = {0: (j, i)}
        for s_ in range(1, k):
            f1[s_] = (rest[s_ - 1],)
            f2[s_] = (rest[s_ - 1],)
        for (tau, m), npl in omega.items():
            t1 = relabel(tau, f1)
            t2 = relabel(tau, f2)
            # monomial: slot 0 -> p_i + p_j, slot s (1..k-2) -> rest[s-1]
            base = [0] * (2 * K)
            for s_ in range(1, k - 1):
                base[2 * rest[s_ - 1]] += m[2 * s_]
                base[2 * rest[s_ - 1] + 1] += m[2 * s_ + 1]
            a0, b0 = (m[0], m[1]) if k >= 2 else (0, 0)
            for c1, r1, s1 in _binom_split(a0):
                for c2, r2, s2 in _binom_split(b0):
                    fm = list(base)
                    fm[2 * i] += r1
                    fm[2 * j] += s1
                    fm[2 * i + 1] += r2
                    fm[2 * j + 1] += s2
                    for cm, c in canon_full(tuple(fm), K).items():
                        cc = sgn * c1 * c2 * c
                        if not (sl and has_singleton(t1)):
                            add_term(out, (t1, cm), npl, cc)
                        if not (sl and has_singleton(t2)):
                            add_term(out, (t2, cm), npl, -cc)
    return out


# ------------------------------------------------------------------ bracket
def _contract(ta, m, tb, n, sl):
    """contraction of slot m-1 of ta with slot n-1 of tb; tb labels j < n-1 -> m-1+j.
    returns list of (tau, Npoly) on labels 0..m+n-3."""
    e, ep = m - 1, n - 1
    ca = next(c for c in ta if e in c)
    cb = next(c for c in tb if ep in c)
    ia = ca.index(e)
    A = list(ca[ia + 1:] + ca[:ia])
    ib = cb.index(ep)
    B = [m - 1 + x for x in (cb[ib + 1:] + cb[:ib])]
    resta = [list(c) for c in ta if c is not ca]
    restb = [[m - 1 + x for x in c] for c in tb if c is not cb]
    others = resta + restb
    res = []
    merged = A + B
    if merged:
        res.append((canon_tau(others + [merged]), {0: Fr(1)}))
    else:
        res.append((canon_tau(others), {1: Fr(1)}))
    if sl:
        npw = -1
        cyc = list(others)
        for part in (A, B):
            if part:
                cyc.append(part)
            else:
                npw += 1
        res.append((canon_tau(cyc), {npw: Fr(-1)}))
    return res


def bracket(a, m, b, n, sl=False):
    """{a, b} for canonical cochains a (arity m), b (arity n); result arity m+n-2."""
    K = m + n - 2
    raw = {}
    sgn = (-1) ** (n - 1)
    for (ta, ma), da in a.items():
        for (tb, mb), db in b.items():
            fm = tuple(ma) + tuple(mb)
            assert len(fm) == 2 * K
            dd = np_mul(da, db)
            for tau, npl in _contract(ta, m, tb, n, sl):
                if sl and has_singleton(tau):
                    continue
                add_term(raw, (tau, fm), np_mul(dd, npl), sgn)
    out = alt_raw(raw, K, sl=sl)
    return scale(out, Fr(1, factorial(m - 1) * factorial(n - 1)))


# ------------------------------------------------------------------ basic cochains
def bracket_poly(i, j, k):
    """[p_i p_j] = p_i^1 p_j^2 - p_i^2 p_j^1 as raw dict fullmono -> int over k slots."""
    out = {}
    m = [0] * (2 * k)
    m[2 * i] += 1
    m[2 * j + 1] += 1
    out[tuple(m)] = out.get(tuple(m), 0) + 1
    m = [0] * (2 * k)
    m[2 * i + 1] += 1
    m[2 * j] += 1
    out[tuple(m)] = out.get(tuple(m), 0) - 1
    return {k_: v for k_, v in out.items() if v}


def poly_mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            m = tuple(u + v for u, v in zip(a, b))
            out[m] = out.get(m, 0) + x * y
    return {k_: v for k_, v in out.items() if v}


def raw_from(poly, W):
    """poly: dict fullmono -> number; W: dict tau -> Npoly.  raw product."""
    raw = {}
    for fm, c in poly.items():
        for tau, npl in W.items():
            add_term(raw, (tau, fm), npl, c)
    return raw


def monomial_full(k, exps):
    """exps: dict (slot, comp in {0,1}) -> power."""
    m = [0] * (2 * k)
    for (s, c), p in exps.items():
        m[2 * s + c] += p
    return tuple(m)


def s3_cochain():
    """kappa(x1,[x2,x3]) = tr(123) - tr(132) (trace form), arity 3, weight (0,0)."""
    return {(((0, 1, 2),), (0, 0, 0, 0)): {0: Fr(1)}, (((0, 2, 1),), (0, 0, 0, 0)): {0: Fr(-1)}}


def restrict_sl(c):
    return {k: d for k, d in c.items() if not has_singleton(k[0])}


def at_N(c, N):
    return {k: np_eval(d, N) for k, d in c.items() if np_eval(d, N) != 0}


def weight_of(c):
    ws = {mono_degree(m) for (t, m) in c}
    return ws
