"""(L) parametric annihilator relations.  A syzygy (a, b) of the graph polynomial,
   sum_i a_i(t) dU/dt_i = b(t) U,   a_i homogeneous of degree k, b of degree k-1,
gives for every base monomial m = t^(beta-1) with sum beta = n q - k + 1 the total derivative
   sum_i d/dt_i ( a_i m U^(-q) ) = [ sum_i d_i(a_i m) - q b m ] U^(-q),
   int_{RP_{>0}} sum_i d_i(F_i) Omega = - sum_i int_{t_i = 0} F_i     (as for (B), cores/msys.py)
so  sum_gamma c_gamma I_q(G; gamma) + sum_i sum_delta d_delta I_q(G/i; delta) = 0,
the face terms being periods of the contraction G/i at the same level q (F_i|_{t_i=0} = a_i|_{t_i=0} m U_{G/i}^(-q)).
Syzygy bases: exact rational nullspace of the linear map (a, b) -> sum a_i dU_i - b U."""
import itertools
from fractions import Fraction as F

import gutil


def upoly(E):
    N = len(E)
    out = {}
    for C in gutil.cotrees(tuple(E)):
        m = [0] * N
        for e in C:
            m[e] = 1
        out[tuple(m)] = 1
    return out


def monos(N, d):
    res = []
    for c in itertools.combinations_with_replacement(range(N), d):
        m = [0] * N
        for i in c:
            m[i] += 1
        res.append(tuple(m))
    return res


def nullspace(cols, ncol):
    """cols: list over unknowns of dict row -> coef; exact nullspace basis (list of dict unknown -> Fraction)"""
    # build row-major
    rows = {}
    for j, col in enumerate(cols):
        for r, v in col.items():
            rows.setdefault(r, {})[j] = F(v)
    pivots = {}          # pivot column -> row dict (normalized, reduced)
    for r, row in rows.items():
        row = dict(row)
        for c in list(row):
            pass
        changed = True
        while changed:
            changed = False
            for c in [c for c in row if c in pivots]:
                if c not in row:
                    continue
                f = row[c]
                for k, v in pivots[c].items():
                    nv = row.get(k, 0) - f * v
                    if nv:
                        row[k] = nv
                    else:
                        row.pop(k, None)
                changed = True
        if not row:
            continue
        c = max(row)
        f = row[c]
        row = {k: v / f for k, v in row.items()}
        for c2, r2 in pivots.items():
            if c in r2:
                g = r2[c]
                for k, v in row.items():
                    nv = r2.get(k, 0) - g * v
                    if nv:
                        r2[k] = nv
                    else:
                        r2.pop(k, None)
        pivots[c] = row
    free = [j for j in range(ncol) if j not in pivots]
    basis = []
    for fcol in free:
        vec = {fcol: F(1)}
        for c, row in pivots.items():
            if fcol in row:
                vec[c] = -row[fcol]
        basis.append(vec)
    return basis


def syzygies(E, k):
    """basis of syzygies of degree k: list of (a: list over i of dict mono->F, b: dict mono->F)"""
    N = len(E)
    U = upoly(E)
    dU = []
    for i in range(N):
        d = {}
        for m, c in U.items():
            if m[i]:
                mm = list(m)
                mm[i] -= 1
                d[tuple(mm)] = d.get(tuple(mm), 0) + c * m[i]
        dU.append(d)
    am = monos(N, k)
    bm = monos(N, k - 1)
    cols, keys = [], []
    for i in range(N):
        for m in am:
            col = {}
            for u, c in dU[i].items():
                t = tuple(x + y for x, y in zip(m, u))
                col[t] = col.get(t, 0) + c
            cols.append(col)
            keys.append(('a', i, m))
    for m in bm:
        col = {}
        for u, c in U.items():
            t = tuple(x + y for x, y in zip(m, u))
            col[t] = col.get(t, 0) - c
        cols.append(col)
        keys.append(('b', None, m))
    out = []
    for vec in nullspace(cols, len(cols)):
        a = [dict() for _ in range(N)]
        b = {}
        for j, v in vec.items():
            kind, i, m = keys[j]
            if kind == 'a':
                a[i][m] = v
            else:
                b[m] = v
        out.append((a, b))
    return out


def contract_face(E, i):
    """edges of G/i and the map old edge index -> new index (i removed)"""
    u, v = E[i]
    ne, mp = [], {}
    for j, (a, b) in enumerate(E):
        if j == i:
            continue
        mp[j] = len(ne)
        ne.append((u if a == v else a, u if b == v else b))
    return ne, mp


def relation(E, syz, beta, q):
    """items [(edges, gamma, coef)] of the (L) relation, gamma = exponent + 1 (period powers);
    None if some term F_i = a_i m U^-q of the total derivative is not integrable, i.e. some
    beta + mu - e_i (mu a monomial of a_i, entry i >= 1) is a divergent power vector (cancellations
    among main terms can hide such terms: check_syz5.log)."""
    a, b = syz
    N = len(E)
    for i in range(N):
        for m in a[i]:
            g = [x + y for x, y in zip(m, beta)]
            g[i] -= 1
            if g[i] >= 1 and not gutil.convergent(tuple(E), g, q):
                return None
    main = {}
    for i in range(N):
        for m, c in a[i].items():
            # d_i ( c t^m t^(beta-1) ) = c (m_i + beta_i - 1) t^(m + beta - 1 - e_i)
            ex = [x + y - 1 for x, y in zip(m, beta)]
            f = ex[i]
            if f:
                ex[i] -= 1
                g = tuple(x + 1 for x in ex)
                main[g] = main.get(g, 0) + c * f
    for m, c in b.items():
        g = tuple(x + y for x, y in zip(m, beta))
        main[g] = main.get(g, 0) - q * c
    items = [(list(E), list(g), v) for g, v in main.items() if v]
    for i in range(N):
        if beta[i] != 1:
            continue
        # face t_i = 0: F_i = a_i m U^-q, only monomials of a_i without t_i survive
        face = {}
        for m, c in a[i].items():
            if m[i]:
                continue
            face[m] = face.get(m, 0) + c
        if not face:
            continue
        ne, mp = contract_face(E, i)
        for m, c in face.items():
            g = [0] * (N - 1)
            for j in range(N):
                if j != i:
                    g[mp[j]] = m[j] + beta[j]
            items.append((ne, g, F(c)))           # + int_face  (moved to the left: sum ... + face = 0)
    return items


def _flat(s, N, k):
    """syzygy -> dict key -> coef, keys ('a', i, mono) / ('b', mono)"""
    a, b = s
    out = {}
    for i in range(N):
        for m, c in a[i].items():
            out[('a', i, m)] = c
    for m, c in b.items():
        out[('b', m)] = c
    return out


def _shift(s, mu):
    a, b = s
    return ([{tuple(x + y for x, y in zip(m, mu)): c for m, c in ai.items()} for ai in a],
            {tuple(x + y for x, y in zip(m, mu)): c for m, c in b.items()})


def minimal_generators(E, kmax):
    """{k: list of syzygies of degree k not in the span of monomial multiples of lower-degree generators}"""
    N = len(E)
    gens = {}
    for k in range(1, kmax + 1):
        full = syzygies(E, k)
        span = []
        for j, gl in gens.items():
            for mu in monos(N, k - j):
                for s in gl:
                    span.append(_flat(_shift(s, mu), N, k))
        # eliminate: reduced echelon of span, then add full vectors not in span
        piv = {}

        def reduce(v):
            v = dict(v)
            while True:
                hit = [c for c in v if c in piv]
                if not hit:
                    return v
                for c in hit:
                    if c not in v:
                        continue
                    f = v[c]
                    for kk, x in piv[c].items():
                        nv = v.get(kk, 0) - f * x
                        if nv:
                            v[kk] = nv
                        else:
                            v.pop(kk, None)

        def insert(v):
            v = reduce(v)
            if not v:
                return False
            c = max(v)
            f = v[c]
            piv[c] = {kk: x / f for kk, x in v.items()}
            return True
        for v in span:
            insert(v)
        new = []
        for s in full:
            if insert(_flat(s, N, k)):
                new.append(s)
        gens[k] = new
    return gens
