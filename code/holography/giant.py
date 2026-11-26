"""1/16-BPS index of twisted N=4 U(N)/SU(N) and Murthy's giant graviton expansion (exact integers).

Level grading: x^{2(E+J_L)}.  Letters (relative complex, constant ghost removed):
    f = 1 - (1-t1)(1-t2)(1-t3)/((1-p)(1-q)),   unrefined t_I = x^2, p = q = x^3.
Series are dicts {exponent tuple: int}; exponent[0] is the level, truncated at level <= L.

Z_N    = Molien-Weyl integral of PE[f chi_adj]                  (Budzik-Gaiotto (129), v = vbar = 0)
Z_inf  = prod_{n>=1} 1/(1 - f(x^n))                             ((130))
G_N^(k) = Z^+_{k|k}[f^!; N],  f^! = -f/(1-f)                     ((132), (135)-(139))
Murthy: Z_N = Z_inf sum_k G_N^(k)                               ((131))
"""
from math import comb
from fractions import Fraction


# ---------------------------------------------------------------- univariate helpers (level only)

def ser_mul1(a, b, L):
    c = [0] * (L + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, L + 1 - i):
                y = b[j]
                if y:
                    c[i + j] += x * y
    return c


def ser_inv1(a, L):
    assert a[0] in (1, -1)
    inv = [0] * (L + 1)
    inv[0] = a[0]
    for n in range(1, L + 1):
        s = 0
        for k in range(1, n + 1):
            if a[k]:
                s += a[k] * inv[n - k]
        inv[n] = -s * a[0]
    return inv


def letter_index_unrefined(L):
    """f_j, j = 0..L, for f = 1 - (1-x^2)^3/(1-x^3)^2."""
    num = [0] * (L + 1)
    for k, c in enumerate([1, -3, 3, -1]):
        if 2 * k <= L:
            num[2 * k] += c
    den = [0] * (L + 1)
    for k, c in enumerate([1, -2, 1]):
        if 3 * k <= L:
            den[3 * k] += c
    r = ser_mul1(num, ser_inv1(den, L), L)
    f = [-v for v in r]
    f[0] += 1
    return f


def dual_letter_index_unrefined(L):
    """f^! = -f/(1-f) = 1 - (1-x^3)^2/(1-x^2)^3."""
    num = [0] * (L + 1)
    for k, c in enumerate([1, -2, 1]):
        if 3 * k <= L:
            num[3 * k] += c
    den = [0] * (L + 1)
    for k, c in enumerate([1, -3, 3, -1]):
        if 2 * k <= L:
            den[2 * k] += c
    r = ser_mul1(num, ser_inv1(den, L), L)
    g = [-v for v in r]
    g[0] += 1
    return g


# ---------------------------------------------------------------- sparse multivariate series

def smul(A, B, L):
    C = {}
    for ea, ca in A.items():
        la = ea[0]
        for eb, cb in B.items():
            if la + eb[0] > L:
                continue
            e = tuple(x + y for x, y in zip(ea, eb))
            v = C.get(e, 0) + ca * cb
            if v:
                C[e] = v
            else:
                C.pop(e, None)
    return C


def pe_factor(m, s, L):
    """Series of (1 - X^m)^(-s) truncated at level L; m an exponent tuple with m[0] > 0."""
    out = {}
    lev = m[0]
    nmax = L // lev
    for n in range(nmax + 1):
        if s > 0:
            c = comb(n + s - 1, s - 1) if s > 0 else 0
        else:
            c = (-1) ** n * comb(-s, n)
        if c:
            out[tuple(n * x for x in m)] = c
    return out


def pe(terms, L, prune=None):
    """PE of sum_m s_m X^m (dict m -> s); every m has m[0] > 0.  prune(e) -> False drops a monomial."""
    dim = len(next(iter(terms)))
    R = {tuple([0] * dim): 1}
    for m in sorted(terms, key=lambda e: e[0]):
        s = terms[m]
        if s == 0:
            continue
        F = pe_factor(m, s, L)
        R = smul(R, F, L)
        if prune is not None:
            R = {e: c for e, c in R.items() if prune(e)}
    return R


def plethystic_terms_unrefined(f, chars, L):
    """Terms of f(x) * sum_g chars[g] Y^g as {(level, g...): coeff}.  chars: dict gauge-exponent -> int."""
    T = {}
    for j in range(1, L + 1):
        if f[j]:
            for g, cg in chars.items():
                e = (j,) + tuple(g)
                T[e] = T.get(e, 0) + f[j] * cg
    return {e: c for e, c in T.items() if c}


def to_dense(S, L):
    out = [0] * (L + 1)
    for e, c in S.items():
        out[e[0]] += c
    return out


def constant_term(S, L):
    out = [0] * (L + 1)
    for e, c in S.items():
        if all(x == 0 for x in e[1:]):
            out[e[0]] += c
    return out


# ---------------------------------------------------------------- Molien-Weyl

def pe_univariate(f, L):
    T = {(j,): f[j] for j in range(1, L + 1) if f[j]}
    return to_dense(pe(T, L), L)


def Z_U(N, L, f=None):
    """Unrefined U(N) index, N = 1, 2, 3."""
    if f is None:
        f = letter_index_unrefined(L)
    if N == 1:
        return pe_univariate(f, L)
    if N == 2:
        # adjoint = 2 + z + 1/z (z = mu1/mu2), Haar (1/2)(1-z)(1-1/z)
        T = plethystic_terms_unrefined(f, {(1,): 1, (-1,): 1}, L)
        P = pe(T, L)
        haar = {(0, 0): 2, (0, 1): -1, (0, -1): -1}
        I = constant_term(smul(P, haar, L), L)
        I = [Fraction(v, 2) for v in I]
        assert all(v.denominator == 1 for v in I)
        I = [int(v) for v in I]
        return ser_mul1(pe_univariate([2 * v for v in f], L), I, L)
    if N == 3:
        # y1 = mu1/mu2, y2 = mu2/mu3; roots +-(1,0), +-(0,1), +-(1,1)
        roots = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
        T = plethystic_terms_unrefined(f, {r: 1 for r in roots}, L)
        P = pe(T, L)
        haar = {(0, 0, 0): 1}
        for r in roots:
            haar = smul(haar, {(0, 0, 0): 1, (0,) + r: -1}, L)
        I = constant_term(smul(P, haar, L), L)
        I = [Fraction(v, 6) for v in I]
        assert all(v.denominator == 1 for v in I)
        I = [int(v) for v in I]
        return ser_mul1(pe_univariate([3 * v for v in f], L), I, L)
    raise ValueError(N)


def Z_SU(N, L):
    f = letter_index_unrefined(L)
    ZU = Z_U(N, L, f)
    return ser_mul1(ZU, ser_inv1(pe_univariate(f, L), L), L)


def Z_inf(L, f=None):
    if f is None:
        f = letter_index_unrefined(L)
    R = [1] + [0] * L
    for n in range(1, L + 1):
        fn = [0] * (L + 1)
        for j in range(1, L // n + 1):
            fn[j * n] = f[j]
        one_minus = [-v for v in fn]
        one_minus[0] += 1
        R = ser_mul1(R, ser_inv1(one_minus, L), L)
    return R


# ---------------------------------------------------------------- Murthy's G_N^(1), G_N^(2)

def G1(N, L, fd=None):
    """G_N^(1) = CT_w[ -sum_n (n+1) w^{N+1+n} PE[f^!(2 - w - 1/w)] ],  |w| < 1."""
    if fd is None:
        fd = dual_letter_index_unrefined(L)
    T = {}
    for j in range(1, L + 1):
        if fd[j]:
            T[(j, 0)] = T.get((j, 0), 0) + 2 * fd[j]
            T[(j, 1)] = T.get((j, 1), 0) - fd[j]
            T[(j, -1)] = T.get((j, -1), 0) - fd[j]
    P = pe(T, L)
    out = [0] * (L + 1)
    for (lev, m), c in P.items():
        # need m = -(N+1+n), n >= 0
        n = -m - (N + 1)
        if n >= 0:
            out[lev] += -(n + 1) * c
    return out


def _Fcoef(A, B, C):
    """Coefficient of a^A b^B c^{-C} in 1/((1-a)^2(1-a/c)^2(1-b)^2(1-b/c)^2)."""
    if A < 0 or B < 0 or C < 0:
        return 0
    s = 0
    for k in range(0, min(A, C) + 1):
        l = C - k
        if l < 0 or l > B:
            continue
        s += (A - k + 1) * (k + 1) * (B - l + 1) * (l + 1)
    return s


def G2(N, L, fd=None):
    """G_N^(2) via mu=(a,b), nu=(1,c):
       (1/4) CT[(ab/c)^{N+2} (1-a/b)(1-b/a)(1-1/c)(1-c) PE[f^! S Sbar] F(a,b,c)]."""
    if fd is None:
        fd = dual_letter_index_unrefined(L)
    # S Sbar = 4 + a/b + b/a + c + 1/c - a - b - 1/a - 1/b - a/c - b/c - c/a - c/b
    SS = {(0, 0, 0): 4, (1, -1, 0): 1, (-1, 1, 0): 1, (0, 0, 1): 1, (0, 0, -1): 1,
          (1, 0, 0): -1, (0, 1, 0): -1, (-1, 0, 0): -1, (0, -1, 0): -1,
          (1, 0, -1): -1, (0, 1, -1): -1, (-1, 0, 1): -1, (0, -1, 1): -1}
    # pe() includes all windings x^n, a^n, b^n, c^n of each monomial term
    T = {}
    for j in range(1, L + 1):
        if fd[j]:
            for g, cg in SS.items():
                e = (j,) + g
                T[e] = T.get(e, 0) + fd[j] * cg
    T = {e: c for e, c in T.items() if c}
    K = N + 2

    def prune(e):
        lev, al, be, ga = e
        R = (L - lev) // 2 + 1          # max further change of each exponent (+1 for the Haar factor)
        if al + K - R > 0:
            return False
        if be + K - R > 0:
            return False
        if ga - K + R < 0:
            return False
        if (ga - K) + (al + K) + (be + K) - R > 0:
            return False
        return True

    P = pe(T, L, prune)
    pref = {(0, K, K, -K): 1}
    haar = {(0, 0, 0, 0): 1}
    for m in [(1, -1, 0), (-1, 1, 0), (0, 0, -1), (0, 0, 1)]:
        haar = smul(haar, {(0, 0, 0, 0): 1, (0,) + m: -1}, L)
    P = smul(smul(P, pref, L), haar, L)
    out = [Fraction(0)] * (L + 1)
    for (lev, al, be, ga), c in P.items():
        Fc = _Fcoef(-al, -be, ga)
        if Fc:
            out[lev] += Fraction(c * Fc, 4)
    assert all(v.denominator == 1 for v in out)
    return [int(v) for v in out]


def murthy_partial(N, L, kmax=2):
    """Z_inf * (1 + G^(1) + ... + G^(kmax)) truncated at L."""
    Zi = Z_inf(L)
    S = [1] + [0] * L
    if kmax >= 1:
        S = [a + b for a, b in zip(S, G1(N, L))]
    if kmax >= 2:
        S = [a + b for a, b in zip(S, G2(N, L))]
    return ser_mul1(Zi, S, L), S


# ---------------------------------------------------------------- refined series in a charge box
# exponent tuple (level, th1, th2, th3, n+, n-, gauge...) truncated to th_I <= box_I, n <= box.

def _level(th, n):
    return 2 * sum(th) + 3 * sum(n)


def letter_terms_refined(box, dual=False):
    """Monomial terms of f (or f^! if dual) inside box = (th1,th2,th3,n+,n-).
    f   = 1 - (1-t1)(1-t2)(1-t3)/((1-p)(1-q))     (letters: c-jets, gamma, beta, b)
    f^! = 1 - (1-p)(1-q)/((1-t1)(1-t2)(1-t3))."""
    out = {}
    if not dual:
        # numerator (1-t1)(1-t2)(1-t3), times sum p^a q^b
        for S in range(8):
            th = tuple((S >> i) & 1 for i in range(3))
            sgn = (-1) ** sum(th)
            for a in range(box[3] + 1):
                for b in range(box[4] + 1):
                    if any(th[i] > box[i] for i in range(3)):
                        continue
                    e = th + (a, b)
                    out[e] = out.get(e, 0) - sgn
        out[(0, 0, 0, 0, 0)] = out.get((0, 0, 0, 0, 0), 0) + 1
    else:
        for a in range(2):
            for b in range(2):
                sgn = (-1) ** (a + b)
                if a > box[3] or b > box[4]:
                    continue
                for m1 in range(box[0] + 1):
                    for m2 in range(box[1] + 1):
                        for m3 in range(box[2] + 1):
                            e = (m1, m2, m3, a, b)
                            out[e] = out.get(e, 0) - sgn
        out[(0, 0, 0, 0, 0)] = out.get((0, 0, 0, 0, 0), 0) + 1
    return {(_level(e[:3], e[3:]),) + e: c for e, c in out.items() if c}


def _inbox(e, box):
    return all(e[1 + i] <= box[i] for i in range(5))


def pe_box(terms, box, gauge_dim=0):
    L = _level(box[:3], box[3:])
    return pe(terms, L, prune=lambda e: _inbox(e, box))


def Z_U2_box(box):
    """Refined U(2) index restricted to the box: dict charge -> coeff."""
    f = letter_terms_refined(box)
    L = _level(box[:3], box[3:])
    T = {}
    for e, c in f.items():
        for g, cg in ((1, 1), (-1, 1), (0, 2)):
            ee = e + (g,)
            T[ee] = T.get(ee, 0) + c * cg
    P = pe({e: c for e, c in T.items() if c}, L, prune=lambda e: _inbox(e, box))
    haar = {tuple([0] * 6) + (0,): 2, tuple([0] * 6) + (1,): -1, tuple([0] * 6) + (-1,): -1}
    out = {}
    for e, c in P.items():
        for h, ch in haar.items():
            if e[-1] + h[-1] == 0:
                out[e[1:6]] = out.get(e[1:6], 0) + c * ch
    assert all(v % 2 == 0 for v in out.values())
    return {k: v // 2 for k, v in out.items() if v}


def Zinf_box(box):
    f = letter_terms_refined(box)
    L = _level(box[:3], box[3:])
    R = {(0,) * 6: 1}
    for n in range(1, L + 1):
        # 1/(1 - f(x^n)) = sum_m f(x^n)^m
        fn = {}
        for e, c in f.items():
            ee = tuple(n * x for x in e)
            if _inbox(ee, box):
                fn[ee] = fn.get(ee, 0) + c
        if not fn:
            continue
        geo = {(0,) * 6: 1}
        powm = {(0,) * 6: 1}
        while True:
            powm = {e: c for e, c in smul(powm, fn, L).items() if _inbox(e, box)}
            if not powm:
                break
            for e, c in powm.items():
                geo[e] = geo.get(e, 0) + c
        R = {e: c for e, c in smul(R, geo, L).items() if _inbox(e, box)}
    return {e[1:]: c for e, c in R.items() if c}


def G1_box(N, box):
    fd = letter_terms_refined(box, dual=True)
    L = _level(box[:3], box[3:])
    T = {}
    for e, c in fd.items():
        for g, cg in ((0, 2), (1, -1), (-1, -1)):
            ee = e + (g,)
            T[ee] = T.get(ee, 0) + c * cg
    P = pe({e: c for e, c in T.items() if c}, L, prune=lambda e: _inbox(e, box))
    out = {}
    for e, c in P.items():
        n = -e[-1] - (N + 1)
        if n >= 0:
            out[e[1:6]] = out.get(e[1:6], 0) - (n + 1) * c
    return {k: v for k, v in out.items() if v}


# ---------------------------------------------------------------- G_N^(3) (pruned constant term)

def _kernel3(A, B, memo={}):
    """sum over 3x3 nonnegative integer matrices n with row sums A, column sums B of prod (n_ij+1)."""
    key = (A, B)
    if key in memo:
        return memo[key]
    if sum(A) != sum(B) or min(A) < 0 or min(B) < 0:
        return 0
    if len(A) == 1:
        r = 1
        for b in B:
            r *= b + 1
        memo[key] = r
        return r
    a0 = A[0]
    tot = 0
    # distribute a0 over columns
    for n0 in range(min(a0, B[0]) + 1):
        for n1 in range(min(a0 - n0, B[1]) + 1):
            n2 = a0 - n0 - n1
            if n2 > B[2]:
                continue
            w = (n0 + 1) * (n1 + 1) * (n2 + 1)
            tot += w * _kernel3(A[1:], (B[0] - n0, B[1] - n1, B[2] - n2))
    memo[key] = tot
    return tot


def G3(N, L, fd=None):
    """G_N^(3): mu = (a1,a2,a3), nu = (1,c2,c3);
       -(1/36) CT[(a1a2a3/(c2c3))^{N+3} Delta(mu) Delta(nu) PE[f^! S Sbar] prod_{ij}(1-mu_i/nu_j)^{-2}]."""
    if fd is None:
        fd = dual_letter_index_unrefined(L)
    K = N + 3
    # exponents (a1,a2,a3,c2,c3); mu_i -> unit vector, nu_1 -> 0, nu_j -> unit c_j
    mu = [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0)]
    nu = [(0, 0, 0, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1)]
    s = [(v, 1) for v in mu] + [(v, -1) for v in nu]
    SS = {}
    for (u, su) in s:
        for (v, sv) in s:
            e = tuple(x - y for x, y in zip(u, v))
            SS[e] = SS.get(e, 0) + su * sv
    T = {}
    for j in range(1, L + 1):
        if fd[j]:
            for e, c in SS.items():
                if c:
                    ee = (j,) + e
                    T[ee] = T.get(ee, 0) + fd[j] * c
    T = {e: c for e, c in T.items() if c}

    def prune(e):
        lev = e[0]
        Rr = (L - lev) // 2 + 2
        a = e[1:4]
        c = e[4:6]
        for x in a:
            if x + K - Rr > 0:
                return False
        for x in c:
            if x - K + Rr < 0:
                return False
        # B_1 = -sum(alpha_tot) - sum(gamma_tot) >= 0 (degree-zero constraint with nu_1 = 1)
        if (sum(a) + 3 * K) + (sum(c) - 2 * K) - Rr > 0:
            return False
        return True

    P = pe(T, L, prune)
    pref = {(0, K, K, K, -K, -K): 1}
    haar = {(0, 0, 0, 0, 0, 0): 1}
    for i in range(3):
        for j in range(3):
            if i != j:
                e1 = tuple(x - y for x, y in zip(mu[i], mu[j]))
                haar = smul(haar, {(0,) * 6: 1, (0,) + e1: -1}, L)
                e2 = tuple(x - y for x, y in zip(nu[i], nu[j]))
                haar = smul(haar, {(0,) * 6: 1, (0,) + e2: -1}, L)
    P = smul(smul(P, pref, L), haar, L)
    out = [Fraction(0)] * (L + 1)
    for e, c in P.items():
        lev = e[0]
        A = tuple(-x for x in e[1:4])
        B23 = e[4:6]
        if min(A) < 0 or min(B23) < 0:
            continue
        B1 = sum(A) - sum(B23)
        if B1 < 0:
            continue
        kc = _kernel3(A, (B1,) + tuple(B23))
        if kc:
            out[lev] += Fraction(-c * kc, 36)
    assert all(v.denominator == 1 for v in out)
    return [int(v) for v in out]

