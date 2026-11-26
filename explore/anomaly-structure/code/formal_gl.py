"""gl_N for general N: invariant cochains in the trace basis.

tau = product of traces, encoded as a sorted tuple of cycles (each cycle a tuple of slot labels,
rotated to start at its minimum):  tau(X_0, .., X_{k-1}) = prod_c tr(prod_{i in c} X_i).
For N >= k the tau's of S_k are linearly independent on gl_N^{(x)k}.

Colour weights by the double-line rules (trace form, dual basis X^e):
    tr(A X_e) tr(B X^e) = tr(A B),   tr(A X_e B X^e) = tr(A) tr(B),   tr(1) = N,
    f(h1,h2,h3) = tr(h1 h2 h3) - tr(h1 h3 h2),
with the orientation sign of cochain.colour_weight.

Cochains: dict (tau, monomial) -> {Npow: int}, canonical modulo momentum conservation.
"""
import itertools
import numpy as np
import sympy as sp
import cochain as cc


def canon_cycle(c):
    c = tuple(c)
    i = c.index(min(c))
    return c[i:] + c[:i]


def canon_tau(cycles):
    return tuple(sorted(canon_cycle(c) for c in cycles))


# ------------------------------------------------------------------ colour weight (double line)
def colour_weight_formal(nv, edges, ext):
    """returns dict tau (labels 0..k-1 = positions in ext) -> {Npow: int}."""
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    local = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v],
                     key=lambda h: (h[1], h[2]))
        local.append(hs)
    Lp = [h for hs in local for h in hs]
    pos = {h: i for i, h in enumerate(L)}
    sgn = cc.perm_sign([pos[h] for h in Lp])
    # terms: (coef, Npow, list of words)
    terms = [(sgn, 0, [])]
    for v in range(nv):
        h1, h2, h3 = local[v]
        new = []
        for (c, npw, ws) in terms:
            new.append((c, npw, ws + [[h1, h2, h3]]))
            new.append((-c, npw, ws + [[h1, h3, h2]]))
        terms = new
    for e in range(len(edges)):
        h, hp = ('e', e, 0), ('e', e, 1)
        new = []
        for (c, npw, ws) in terms:
            i = next(k for k, w in enumerate(ws) if h in w)
            j = next(k for k, w in enumerate(ws) if hp in w)
            ws = [list(w) for w in ws]
            if i != j:
                w1, w2 = ws[i], ws[j]
                a = w1.index(h); w1 = w1[a + 1:] + w1[:a]          # A with h removed (A h cyclic)
                b = w2.index(hp); w2 = w2[b + 1:] + w2[:b]         # B
                merged = w1 + w2
                rest = [w for k, w in enumerate(ws) if k not in (i, j)]
                if merged:
                    rest.append(merged)
                    new.append((c, npw, rest))
                else:
                    new.append((c, npw + 1, rest))
            else:
                w = ws[i]
                a = w.index(h)
                w = w[a:] + w[:a]                                 # h B hp A
                b = w.index(hp)
                B = w[1:b]
                A = w[b + 1:]
                rest = [x for k, x in enumerate(ws) if k != i]
                add = 0
                for part in (A, B):
                    if part:
                        rest.append(part)
                    else:
                        add += 1
                new.append((c, npw + add, rest))
        terms = new
    out = {}
    for (c, npw, ws) in terms:
        cyc = [tuple(ext.index(h[1]) for h in w) for w in ws]
        tau = canon_tau(cyc)
        d = out.setdefault(tau, {})
        d[npw] = d.get(npw, 0) + c
    return {t: {p: v for p, v in d.items() if v} for t, d in out.items() if any(d.values())}


def evaluate_tau(tau, mats):
    val = 1
    for c in tau:
        m = np.eye(mats[0].shape[0], dtype=np.int64)
        for i in c:
            m = m @ mats[i]
        val *= int(np.trace(m))
    return val


def tau_tensor(tau, g, k):
    """explicit tensor tau[a_0..a_{k-1}] = prod_c tr(prod_{i in c} X_{a_i}) on the basis g.mats."""
    M = np.array(g.mats, dtype=np.int64)          # (dim, N, N)
    letters = 'abcdefghij'
    mid = 'klmnopqrstuvwxyzABCDEFGH'
    terms, ops, used = [], [], 0
    for c in tau:
        idx = [mid[used + j] for j in range(len(c))]
        used += len(c)
        for j, i in enumerate(c):
            terms.append(letters[i] + idx[j] + idx[(j + 1) % len(c)])
            ops.append(M)
    return np.einsum(','.join(terms) + '->' + letters[:k], *ops)


def tensor_from_formal(W, g, k):
    """explicit tensor W[a_1..a_k] on the basis g.mats, N = matrix size."""
    N = g.mats[0].shape[0]
    T = np.zeros((g.dim,) * k, dtype=np.int64)
    for tau, d in W.items():
        T += sum(v * N ** p for p, v in d.items()) * tau_tensor(tau, g, k)
    return T


# ------------------------------------------------------------------ cochains
def relabel(tau, f):
    """tau with each label i replaced by the sequence f[i] (tuple)."""
    cyc = []
    for c in tau:
        seq = []
        for i in c:
            seq += list(f[i])
        cyc.append(seq)
    return canon_tau(cyc)


def add(out, key, d, s=1):
    e = out.setdefault(key, {})
    for p, v in d.items():
        e[p] = e.get(p, 0) + s * v


def clean(c):
    return {k: {p: v for p, v in d.items() if v} for k, d in c.items() if any(d.values())}


_ALT = {}
_DCE = {}


def alternate_formal(mom, P, W, k):
    """canonical form of sum_sigma sgn(sigma) P(p_sigma) W(X_sigma); W: dict tau -> {Npow: int}."""
    out = {}
    for sigma in itertools.permutations(range(k)):
        s = cc.perm_sign(sigma)
        key = (id(mom), sp.srepr(P), sigma, k)
        if key not in _ALT:
            sub = {mom.p[i][c]: mom.p[sigma[i]][c] for i in range(k) for c in range(mom.H)}
            Ps = sp.sympify(P).subs(sub, simultaneous=True)
            _ALT[key] = mom.canon(Ps, k)
        can = _ALT[key]
        for tau, d in W.items():
            ts = relabel(tau, {i: (sigma[i],) for i in range(k)})
            for m, coef in can.items():
                q = sp.Rational(coef)
                assert q.q == 1
                add(out, (ts, m), d, s * int(q.p))
    return clean(out)


def d_formal(mom, omega, k):
    out = {}
    Pm = mom.p
    for i, j in itertools.combinations(range(k + 1), 2):
        rest = [r for r in range(k + 1) if r not in (i, j)]
        args = [[Pm[i][c] + Pm[j][c] for c in range(mom.H)]] + [[Pm[r][c] for c in range(mom.H)] for r in rest]
        sgn = (-1) ** (i + j)
        f1 = {0: (i, j)}
        f2 = {0: (j, i)}
        for s_, r in enumerate(rest):
            f1[s_ + 1] = (r,)
            f2[s_ + 1] = (r,)
        for (tau, m), d in omega.items():
            key = (id(mom), k, i, j, m)
            if key not in _DCE:
                _DCE[key] = mom.canon(mom.mono_expr(m, k, args[:k - 1]), k + 1)
            t1 = relabel(tau, f1)
            t2 = relabel(tau, f2)
            for M, coef in _DCE[key].items():
                q = sp.Rational(coef)
                assert q.q == 1
                add(out, (t1, M), d, sgn * int(q.p))
                add(out, (t2, M), d, -sgn * int(q.p))
    return clean(out)


def all_taus(k):
    out = set()
    for s in itertools.permutations(range(k)):
        seen, cyc = set(), []
        for i in range(k):
            if i in seen:
                continue
            c, j = [], i
            while j not in seen:
                seen.add(j); c.append(j); j = s[j]
            cyc.append(c)
        out.add(canon_tau(cyc))
    return sorted(out)


def at_N(c, N):
    """numeric coefficients at a given N: dict key -> int."""
    return {k: sum(v * N ** p for p, v in d.items()) for k, d in c.items()}
