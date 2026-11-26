"""(Copied unchanged from 6d-holomorphic-twist-note/code/salg.py.)
Graded-commutative polynomial algebra for the letters of holomorphic CS on C^{3|2}.

Letters: (t, A, k1, k2, k3); t: 0 = theta (A = 1, 2), 1 = c, 2 = b (odd); 3 = gamma^1, 4 = gamma^2 (even).
Polynomial: dict {(odd_tuple, even_tuple): Fraction}, both tuples sorted.
Superfield: C^A(theta) = c^A + theta_1 gamma^{1A} + theta_2 gamma^{2A} + theta_1 theta_2 b^A.
"""
from fractions import Fraction as Fr
from itertools import product

TH1 = (0, 1, 0, 0, 0)
TH2 = (0, 2, 0, 0, 0)
ODD_T = (0, 1, 2)


def is_odd(l):
    return l[0] in ODD_T


def const(x=1):
    return {((), ()): Fr(x)}


def gen(l):
    return {((l,), ()): Fr(1)} if is_odd(l) else {((), (l,)): Fr(1)}


def add(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        v = r.get(m, 0) + s * c
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def scale(p, s):
    return {m: c * s for m, c in p.items()} if s else {}


def _merge_odd(a, b):
    """Sorted merge of two sorted odd tuples; returns (tuple, sign) or None if a generator repeats."""
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


def mul(p, q):
    r = {}
    for (o1, e1), c1 in p.items():
        for (o2, e2), c2 in q.items():
            mo = _merge_odd(o1, o2)
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


def lderiv(p, l):
    """Left derivative with respect to the letter l."""
    r = {}
    if is_odd(l):
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


def letters_of(p):
    s = set()
    for (o, e) in p:
        s.update(o); s.update(e)
    return s


def shift(l, a):
    k = list(l[2:]); k[a] += 1
    return (l[0], l[1]) + tuple(k)


def dpar(p, a):
    """Holomorphic derivative d/dx^a (even derivation raising the jet index)."""
    r = {}
    for l in letters_of(p):
        if l[0] == 0:
            continue
        r = add(r, mul(gen(shift(l, a)), lderiv(p, l)))
    return r


def csuper(A, k=(0, 0, 0)):
    """Superfield jet d^k C^A(theta)."""
    th1, th2 = gen(TH1), gen(TH2)
    r = gen((1, A) + tuple(k))
    r = add(r, mul(th1, gen((3, A) + tuple(k))))
    r = add(r, mul(th2, gen((4, A) + tuple(k))))
    r = add(r, mul(mul(th1, th2), gen((2, A) + tuple(k))))
    return r


def theta_coeff(p, which):
    """Coefficient of 1 ('0'), theta_1 ('1'), theta_2 ('2'), theta_1 theta_2 ('12'), thetas standing on the left."""
    want = {'0': (), '1': (TH1,), '2': (TH2,), '12': (TH1, TH2)}[which]
    n = len(want)
    r = {}
    for (o, e), c in p.items():
        nth = sum(1 for l in o if l[0] == 0)
        if nth == n and o[:n] == want:
            r[(o[n:], e)] = c
    return r


def THETA(t):
    """theta factor of delta/delta C(theta) on the component of type t:
    delta/delta C(theta) = th1 th2 d/dc + th2 d/dgamma^1 - th1 d/dgamma^2 + d/db  (odd derivation),
    so that delta C(theta')/delta C(theta) = (th1 - th1')(th2 - th2')."""
    if t == 1:
        return mul(gen(TH1), gen(TH2))
    if t == 3:
        return gen(TH2)
    if t == 4:
        return scale(gen(TH1), -1)
    return const(1)


def sderiv(p, l):
    """Component l of delta/delta C^A(z, theta): THETA(type) times the left derivative."""
    d = lderiv(p, l)
    return mul(THETA(l[0]), d) if d else {}


def eps3():
    out = {}
    for a, b, c in product(range(3), repeat=3):
        s = (a - b) * (b - c) * (c - a)
        if s:
            out[(a, b, c)] = 1 if s in (2, -4) or (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
    for key in out:
        out[key] = 1 if key in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
    return out


EPS = eps3()


def structure_constants(name):
    """Totally antisymmetric f_{ABC} as dict; indices from 0."""
    if name == 'su2':
        return 3, {k: Fr(v) for k, v in EPS.items()}
    if name == 'su3':
        # f_{ABC} of su(3) in the Gell-Mann basis, rescaled by 2/sqrt(3)-free rational trick:
        # use the equivalent real form with rational structure constants obtained from gl(3) commutators.
        raise NotImplementedError
    if name == 'so4':
        n, f = 6, {}
        for off in (0, 3):
            for k, v in EPS.items():
                f[tuple(i + off for i in k)] = Fr(v)
        return n, f
    if name.startswith('so') and name[2:].isdigit():
        N = int(name[2:])
        idx = {}
        for i in range(N):
            for j in range(i + 1, N):
                idx[(i, j)] = len(idx)
        def Lmat(i, j):
            M = [[0] * N for _ in range(N)]
            M[i][j] = 1; M[j][i] = -1
            return M
        def comm(X, Y):
            return [[sum(X[a][k] * Y[k][b] - Y[a][k] * X[k][b] for k in range(N)) for b in range(N)] for a in range(N)]
        f = {}
        gens = {v: Lmat(*k) for k, v in idx.items()}
        for B in gens:
            for C in gens:
                Z = comm(gens[B], gens[C])
                for (i, j), A in idx.items():
                    if Z[i][j]:
                        f[(A, B, C)] = Fr(Z[i][j])
        # total antisymmetry check
        for (A, B, C), v in f.items():
            assert f.get((B, A, C), 0) == -v and f.get((A, C, B), 0) == -v
        return len(idx), f
    raise ValueError(name)
