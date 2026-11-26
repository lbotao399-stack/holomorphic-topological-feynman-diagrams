"""Completed K5 = conformal K6 family in d = 2q (position space, beta = q - alpha).

J_q(beta) := Pi(K6; beta) = Per_x(K6 - c; beta|_{K6-c})  (any c; symmetry/REPORT.md 1.2)
          = Gamma(q) I_q(K6 - c; q - beta) / prod_{e in K6-c} Gamma(beta_e)   (if beta in [1, q-1] on K6 - c)
balance: sum_{e at v} beta_e = 2q for all 6 vertices; convergence: beta_e < q (15 edges), triangle sums < 2q.
Vertices 0..5, edges (i<j) in lexicographic order.
"""
import itertools
from fractions import Fraction as F

import numpy as np

V = 6
EDGES = [(i, j) for i in range(V) for j in range(i + 1, V)]
EIDX = {}
for _k, (_i, _j) in enumerate(EDGES):
    EIDX[(_i, _j)] = _k
    EIDX[(_j, _i)] = _k
TRI = [(a, b, c) for a, b, c in itertools.combinations(range(V), 3)]
TRI_E = [(EIDX[a, b], EIDX[a, c], EIDX[b, c]) for a, b, c in TRI]
STAR = [[EIDX[v, w] for w in range(V) if w != v] for v in range(V)]

PERMS = list(itertools.permutations(range(V)))
# beta'[k] = beta[PM[p][k]]  for the relabelling v -> sigma(v):  beta'_{sigma i, sigma j} = beta_{ij}
PM = np.zeros((len(PERMS), 15), dtype=np.int64)
for _p, _s in enumerate(PERMS):
    for _k, (_i, _j) in enumerate(EDGES):
        PM[_p, EIDX[_s[_i], _s[_j]]] = _k


def balanced(b, q):
    return all(sum(b[k] for k in STAR[v]) == 2 * q for v in range(V))


def convergent(b, q):
    if max(b) >= q:
        return False
    return all(b[x] + b[y] + b[z] < 2 * q for x, y, z in TRI_E)


def is_T1(b):
    return 0 in b


def is_T2(b, q):
    return any(b[x] + b[y] + b[z] == q for x, y, z in TRI_E)


def canon(b):
    b = np.asarray(b, dtype=np.int64)
    W = b[PM]                     # (720, 15)
    order = np.lexsort(W.T[::-1])
    return tuple(int(x) for x in W[order[0]])


def canon_batch(B, q):
    """canonical forms (lexicographic minimum over S6) of the rows of B (N x 15)"""
    B = np.asarray(B, dtype=np.int64)
    off = 2 * q
    base = 3 * q + 1
    X = B[:, PM] + off            # (N, 720, 15), digits in [0, base)
    k1 = np.zeros(X.shape[:2], dtype=np.int64)
    k2 = np.zeros(X.shape[:2], dtype=np.int64)
    for k in range(8):
        k1 = k1 * base + X[:, :, k]
    for k in range(8, 15):
        k2 = k2 * base + X[:, :, k]
    m1 = k1.min(axis=1)
    k2m = np.where(k1 == m1[:, None], k2, np.iinfo(np.int64).max)
    am = k2m.argmin(axis=1)
    out = X[np.arange(len(B)), am] - off
    return [tuple(int(x) for x in r) for r in out]


def enumerate_convergent(q, lo=None):
    """all balanced convergent integer beta (including beta <= 0); generator of tuples"""
    if lo is None:
        lo = 4 - 2 * q
    hi = q - 1
    b = [None] * 15
    # fill order: vertex 0 edges, vertex 1 remaining, ...
    order = []
    for v in range(V):
        order.append([EIDX[v, w] for w in range(v + 1, V)])
    tri_by_last = {}
    pos = {}
    cnt = 0
    for v in range(V):
        for k in order[v]:
            pos[k] = cnt
            cnt += 1
    for t in TRI_E:
        last = max(t, key=lambda k: pos[k])
        tri_by_last.setdefault(last, []).append(t)

    def tri_ok(k):
        for x, y, z in tri_by_last.get(k, []):
            if b[x] + b[y] + b[z] >= 2 * q:
                return False
        return True

    def rec(v):
        if v == V - 1:
            # vertex 5: all its edges fixed; check balance
            if sum(b[k] for k in STAR[5]) == 2 * q:
                yield tuple(b)
            return
        ks = order[v]
        fixed = sum(b[EIDX[v, w]] for w in range(v))
        rem = 2 * q - fixed
        n = len(ks)

        def fill(i, rem):
            if i == n - 1:
                if lo <= rem <= hi:
                    b[ks[i]] = rem
                    if tri_ok(ks[i]):
                        yield from rec(v + 1)
                    b[ks[i]] = None
                return
            left = n - i - 1
            for x in range(max(lo, rem - hi * left), min(hi, rem - lo * left) + 1):
                b[ks[i]] = x
                if tri_ok(ks[i]):
                    yield from fill(i + 1, rem - x)
                b[ks[i]] = None
        yield from fill(0, rem)

    yield from rec(0)


def classes(q, batch=200000):
    """dict canonical beta -> orbit size (number of labelled vectors)"""
    out = {}
    buf = []
    for b in enumerate_convergent(q):
        buf.append(b)
        if len(buf) >= batch:
            for c in canon_batch(buf, q):
                out[c] = out.get(c, 0) + 1
            buf = []
    if buf:
        for c in canon_batch(buf, q):
            out[c] = out.get(c, 0) + 1
    return out


def k5_decompletions(b, q):
    """vertices c such that K6 - c has all beta in [1, q-1] (a K5 period with alpha in [1, q-1])"""
    out = []
    for c in range(V):
        if all(1 <= b[k] <= q - 1 for k, (i, j) in enumerate(EDGES) if c not in (i, j)):
            out.append(c)
    return out


def kind(b, q):
    if is_T1(b):
        return 'T1'
    if is_T2(b, q):
        return 'T2'
    return 'G'


# ---------------- (X) position-space vertex IBP in K6 form ----------------
def x_relations(b, q):
    """(beta_{xc} - beta_{x x1}) J(b) = sum_i beta_{x xi} [J(b - e_{x x1} + e_{x xi} + e_{x1 c} - e_{xi c})
                                                         - J(b + e_{x xi} - e_{x1 xi} + e_{x1 c} - e_{xc})]
    for ordered distinct (x, x1, c); i over the remaining 3 vertices.  Returns list of [(coef, beta)], only
    relations whose terms are all convergent."""
    out = []
    for x, x1, c in itertools.permutations(range(V), 3):
        others = [i for i in range(V) if i not in (x, x1, c)]
        items = [(F(b[EIDX[x, c]] - b[EIDX[x, x1]]), tuple(b))]
        ok = True
        for i in others:
            w = b[EIDX[x, i]]
            if w == 0:
                continue
            t1 = list(b)
            t1[EIDX[x, x1]] -= 1
            t1[EIDX[x, i]] += 1
            t1[EIDX[x1, c]] += 1
            t1[EIDX[i, c]] -= 1
            t2 = list(b)
            t2[EIDX[x, i]] += 1
            t2[EIDX[x1, i]] -= 1
            t2[EIDX[x1, c]] += 1
            t2[EIDX[x, c]] -= 1
            for t, s in ((t1, -1), (t2, 1)):
                if not convergent(t, q):
                    ok = False
                    break
                items.append((F(s * w), tuple(t)))
            if not ok:
                break
        if ok:
            out.append(items)
    return out


# ---------------- (DS) dimension shift in K6 form ----------------
def spanning_trees_K5(verts):
    vs = list(verts)
    es = [(a, b) for a, b in itertools.combinations(vs, 2)]
    out = []
    for T in itertools.combinations(es, len(vs) - 1):
        # check acyclic & spanning via union-find
        par = {v: v for v in vs}

        def f(v):
            while par[v] != v:
                par[v] = par[par[v]]
                v = par[v]
            return v
        ok = True
        for a, bb in T:
            ra, rb = f(a), f(bb)
            if ra == rb:
                ok = False
                break
            par[ra] = rb
        if ok:
            out.append(T)
    return out


_TREES = {c: spanning_trees_K5([v for v in range(V) if v != c]) for c in range(V)}


def ds_relation(b, q, c):
    """J_{q-1}(b) = 1/(q-1) sum_T prod_{e in T} b_e J_q(b + 1_T, completed),  T spanning trees of K6 - c.
    b at level q-1 (balanced with 2(q-1)).  Returns (items_q, ok) with items_q = [(coef, beta_q)];
    None if some term with nonzero coefficient is divergent."""
    items = {}
    for T in _TREES[c]:
        coef = 1
        for e in T:
            coef *= b[EIDX[e]]
        if coef == 0:
            continue
        nb = list(b)
        for e in T:
            nb[EIDX[e]] += 1
        # completion edges at c: vertex sums must be 2q
        for v in range(V):
            if v == c:
                continue
            s = sum(nb[k] for k in STAR[v])
            nb[EIDX[v, c]] += 2 * q - s
        assert balanced(nb, q)
        if not convergent(nb, q):
            return None
        t = tuple(nb)
        items[t] = items.get(t, 0) + F(coef, q - 1)
    return list((v, k) for k, v in items.items())
