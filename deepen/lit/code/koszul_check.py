# Exact check (sympy DomainMatrix over QQ): for M_k = sign-twisted S_k-invariants of
# C[p_1..p_k] (p_j in C^2) (x) W^{(x)k}, W = C^w, the Koszul complex of (P^1,P^2),
# P^i = sum_j p_j^i, is exact in degrees 1,2 in bidegree (a,b), and the cokernel
# M/(P^1,P^2)M has the dimension of the restriction of M to the hyperplane sum_j p_j = 0.
import itertools, sys
from sympy import symbols, Poly, Rational, expand
from sympy.polys.matrices import DomainMatrix
from sympy import QQ

def sign(perm):
    s, p = 1, list(perm)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def run(k, w, a, b):
    X = symbols(f'x1:{k+1}'); Y = symbols(f'y1:{k+1}')   # p_j = (x_j, y_j)
    def monos(deg, V):
        return [m for m in itertools.product(range(deg+1), repeat=k) if sum(m) == deg]
    def basis(a, b):
        # spanning set: antisymmetrize x^m y^n (x) e_{t1..tk}
        vecs = {}
        for mx in monos(a, X):
            for my in monos(b, Y):
                for t in itertools.product(range(w), repeat=k):
                    v = {}
                    for g in itertools.permutations(range(k)):
                        s = sign(g)
                        # g acts: argument j -> g(j) simultaneously on momenta and tensor slots
                        key = (tuple(mx[g.index(j)] for j in range(k)),
                               tuple(my[g.index(j)] for j in range(k)),
                               tuple(t[g.index(j)] for j in range(k)))
                        v[key] = v.get(key, 0) + s
                    v = {kk: c for kk, c in v.items() if c != 0}
                    if v:
                        vecs[frozenset(v.items())] = v
        return list(vecs.values())
    def keys(a, b):
        return [(mx, my, t) for mx in monos(a, X) for my in monos(b, Y)
                for t in itertools.product(range(w), repeat=k)]
    def mat(vs, a, b):
        ks = keys(a, b); idx = {kk: i for i, kk in enumerate(ks)}
        rows = [[QQ(0)]*len(ks) for _ in vs]
        for r, v in enumerate(vs):
            for kk, c in v.items():
                rows[r][idx[kk]] = QQ(c)
        return DomainMatrix(rows, (len(vs), len(ks)), QQ) if vs else None
    def indep(vs, a, b):
        if not vs: return []
        M = mat(vs, a, b); rref, piv = M.transpose().rref()
        return [vs[i] for i in piv]
    def mulP(v, comp):
        out = {}
        for (mx, my, t), c in v.items():
            for j in range(k):
                if comp == 0:
                    m2 = list(mx); m2[j] += 1; kk = (tuple(m2), my, t)
                else:
                    m2 = list(my); m2[j] += 1; kk = (mx, tuple(m2), t)
                out[kk] = out.get(kk, 0) + c
        return out
    B = {}
    for (aa, bb) in [(a, b), (a-1, b), (a, b-1), (a-1, b-1)]:
        B[(aa, bb)] = indep(basis(aa, bb), aa, bb) if aa >= 0 and bb >= 0 else []
    dab = len(B[(a, b)]); d10 = len(B[(a-1, b)]); d01 = len(B[(a, b-1)]); d11 = len(B[(a-1, b-1)])
    # d1 : M(a-1,b) + M(a,b-1) -> M(a,b), (m1,m2) -> P^1 m1 + P^2 m2
    imgs = [mulP(v, 0) for v in B[(a-1, b)]] + [mulP(v, 1) for v in B[(a, b-1)]]
    r1 = mat(imgs, a, b).rank() if imgs else 0
    # d2 : M(a-1,b-1) -> M(a-1,b) + M(a,b-1), m -> (-P^2 m, P^1 m); rank via stacked coordinates
    if d11:
        ks1 = keys(a-1, b); ks2 = keys(a, b-1)
        i1 = {kk: i for i, kk in enumerate(ks1)}; i2 = {kk: i for i, kk in enumerate(ks2)}
        rows = []
        for v in B[(a-1, b-1)]:
            row = [QQ(0)]*(len(ks1)+len(ks2))
            for kk, c in mulP(v, 1).items(): row[i1[kk]] -= QQ(c)
            for kk, c in mulP(v, 0).items(): row[len(ks1)+i2[kk]] += QQ(c)
            rows.append(row)
        r2 = DomainMatrix(rows, (len(rows), len(ks1)+len(ks2)), QQ).rank()
    else:
        r2 = 0
    # restriction to hyperplane: substitute p_k = -(p_1+...+p_{k-1}); rank of the image
    def restrict(v):
        out = {}
        for (mx, my, t), c in v.items():
            e = c
            for j in range(k-1):
                e *= X[j]**mx[j] * Y[j]**my[j]
            e *= (-sum(X[:k-1]))**mx[k-1] * (-sum(Y[:k-1]))**my[k-1]
            P = Poly(expand(e), *X[:k-1], *Y[:k-1])
            for mon, cc in P.terms():
                kk = (mon, t); out[kk] = out.get(kk, 0) + cc
        return out
    res = [restrict(v) for v in B[(a, b)]]
    allk = sorted({kk for v in res for kk in v})
    ix = {kk: i for i, kk in enumerate(allk)}
    rows = [[QQ(0)]*len(allk) for _ in res]
    for r, v in enumerate(res):
        for kk, c in v.items(): rows[r][ix[kk]] = QQ(int(c))
    rres = DomainMatrix(rows, (len(res), len(allk)), QQ).rank() if res and allk else 0
    exact_mid = (r1 == d10 + d01 - r2)
    inj = (r2 == d11)
    print(f"k={k} w={w} (a,b)=({a},{b}): dim M_ab={dab}, M_(a-1,b)={d10}, M_(a,b-1)={d01}, M_(a-1,b-1)={d11}; "
          f"rank d2={r2} (injective {inj}), rank d1={r1} (exact in middle {exact_mid}); "
          f"dim coker={dab-r1}, rank of restriction to sum p=0: {rres}, equal {dab-r1==rres}")
    return inj and exact_mid and (dab - r1 == rres)

ok = True
for (k, w, a, b) in [(3, 2, 1, 1), (3, 2, 2, 2), (4, 1, 2, 2), (4, 2, 1, 1), (4, 2, 2, 1), (4, 2, 2, 2)]:
    ok &= run(k, w, a, b)
print("all checks:", ok)
