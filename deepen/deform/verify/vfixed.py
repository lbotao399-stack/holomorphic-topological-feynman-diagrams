"""Faithful evaluation at fixed N (gl1..gl5, sl2..sl5): H^3 at total weight 2, H^4_(2,2), classes,
quadratic map.  Trace relations at N handled by evaluation at random integer points, certified exact:
rank of the evaluation matrix (exact, Fractions) == dimension of span of tau's as functions
(gl: rank of the Gram matrix G[s,r] = #{a in [N]^k : a o s = a o r}; sl: binomial inversion)."""
import sys, time, random
from math import gcd, comb
from fractions import Fraction as Fr
import numpy as np
from vce import *
from vce import _polymul
from vla import Echelon, kernel, rank_of
from vw import W133, W223

sys.setrecursionlimit(10000)
rng = random.Random(7)
t0 = time.time()

# --------------------------------------------------------------- faithful reduction
_gl_dim = {}


def gl_dim(N, k):
    key = (N, k)
    if key in _gl_dim:
        return _gl_dim[key]
    if k == 0:
        _gl_dim[key] = 1
        return 1
    P = ALLPERMS[k]
    import itertools
    G = [[0] * len(P) for _ in P]
    for a in itertools.product(range(N), repeat=k):
        img = {}
        for i, s in enumerate(P):
            v = tuple(a[s[t]] for t in range(k))
            img.setdefault(v, []).append(i)
        for grp in img.values():
            for i in grp:
                for j in grp:
                    G[i][j] += 1
    r = rank_of([{j: Fr(x) for j, x in enumerate(row) if x} for row in G])
    _gl_dim[key] = r
    return r


def sl_dim(N, k):
    return sum((-1) ** (k - j) * comb(k, j) * gl_dim(N, j) for j in range(k + 1))


def tau_np(sig, X):
    v = 1
    for c in perm_cycles(sig):
        M = X[c[0]]
        for t in c[1:]:
            M = M @ X[t]
        v *= int(np.trace(M))
    return v


_L = {}


def Lmat(N, k, sl):
    key = (N, k, sl)
    if key in _L:
        return _L[key]
    P = ALLPERMS[k]
    target = sl_dim(N, k) if sl else gl_dim(N, k)
    rows = []
    E = Echelon()
    tries = 0
    while E.rank() < target:
        tries += 1
        X = []
        for _ in range(k):
            M = np.array([[rng.randint(-2, 2) for _ in range(N)] for _ in range(N)], dtype=np.int64)
            if sl:
                M[N - 1, N - 1] = -sum(M[i, i] for i in range(N - 1))
            X.append(M)
        row = [tau_np(s, X) for s in P]
        ok, _ = E.add({i: Fr(x) for i, x in enumerate(row) if x})
        if ok:
            rows.append(row)
        assert tries < 50 * target + 100, (N, k, sl, E.rank(), target)
    L = np.array(rows, dtype=np.int64)  # r x |S_k|
    _L[key] = (L, target)
    return _L[key]


def red(w, N, sl):
    """formal Laurent cochain -> faithful coordinates dict (mono, j) -> Fraction"""
    if not w:
        return {}
    k = arity(w)
    L, r = Lmat(N, k, sl)
    pidx = {s: i for i, s in enumerate(ALLPERMS[k])}
    vals = {}
    for (s, m, n), c in w.items():
        key = (pidx[s], m)
        vals[key] = vals.get(key, 0) + c * Fr(N) ** n
    den = 1
    for v in vals.values():
        den = den * v.denominator // gcd(den, v.denominator)
    ms = sorted({m for (_, m) in vals})
    midx = {m: i for i, m in enumerate(ms)}
    C = np.zeros((len(ALLPERMS[k]), len(ms)), dtype=object)
    for (i, m), v in vals.items():
        C[i, midx[m]] += int(v * den)
    R = L.astype(object) @ C
    out = {}
    for j in range(R.shape[0]):
        for mi, m in enumerate(ms):
            x = R[j, mi]
            if x != 0:
                out[(m, j)] = Fr(int(x), den)
    return out


# --------------------------------------------------------------- formal ingredients
def cbasis(k, a, b, sl):
    E = Echelon()
    basis = []
    ps = [p for p in ALLPERMS[k] if not (sl and has_fixed(p))]
    for m in monos(k, a, b):
        for s in ps:
            w = alt({(s, m, 0): Fr(1)}, k)
            if sl:
                w = restrict_sl(w)
            if not w:
                continue
            ok, _ = E.add(w)
            if ok:
                basis.append(w)
    return basis


def J_of(Tl, sl):
    w = from_poly_tau(bracket_pq(3, 0, 1), Tl)
    return restrict_sl(w) if sl else w


T_D = [(1, (1, 2, 0)), (1, (2, 0, 1))]
T_S1 = [(1, (0, 2, 1)), (1, (2, 1, 0)), (1, (1, 0, 2))]
T_S0 = [(1, (0, 1, 2))]


def F_of(W, sl):
    P = polyprod(bracket_pq(4, 0, 1), bracket_pq(4, 2, 3))
    out = {}
    for mm, c in P.items():
        for (sig, n), v in W.items():
            add(out, (sig, mm, n), Fr(c) * v)
    out = alt(out, 4)
    return restrict_sl(out) if sl else out


FORMAL = {}
for sl in (False, True):
    FORMAL[sl] = dict(
        C3={(a, b): cbasis(3, a, b, sl) for (a, b) in ((1, 1), (2, 0), (0, 2), (2, 2))},
        C4=cbasis(4, 2, 2, sl),
        F133=F_of(W133(sl), sl), F223=F_of(W223(sl), sl),
        J={"D": J_of(T_D, sl), "S1": J_of(T_S1, sl), "S0": J_of(T_S0, sl)} if not sl else {"D": J_of(T_D, sl)},
    )
    FORMAL[sl]["dC3"] = {ab: [d(w, sl) for w in FORMAL[sl]["C3"][ab]] for ab in FORMAL[sl]["C3"]}
    FORMAL[sl]["dC4"] = [d(w, sl) for w in FORMAL[sl]["C4"]]
    names = list(FORMAL[sl]["J"])
    BR = {}
    for i, x in enumerate(names):
        for y in names[i:]:
            BR[(x, y)] = bracket(FORMAL[sl]["J"][x], FORMAL[sl]["J"][y], sl)
    FORMAL[sl]["BR"] = BR
print(f"formal ingredients built ({time.time() - t0:.1f}s)")
sys.stdout.flush()

import sympy
x0, x1, x2 = sympy.symbols('x0 x1 x2')
XS = [x0, x1, x2]

ALGS = [(1, False), (2, False), (3, False), (4, False), (5, False), (2, True), (3, True), (4, True), (5, True)]
if len(sys.argv) > 1:
    ALGS = [a for a in ALGS if f"{'sl' if a[1] else 'gl'}{a[0]}" in sys.argv[1:]]

for (N, sl) in ALGS:
    nm = f"{'sl' if sl else 'gl'}{N}"
    F = FORMAL[sl]
    print(f"==================== {nm} ====================")
    print(f"  tau-span dims (exact): " + ", ".join(f"k={k}: {sl_dim(N, k) if sl else gl_dim(N, k)}" for k in range(2, 6)))
    # H^3 at weights (1,1),(2,0),(0,2)
    H3basis = None
    for ab in ((1, 1), (2, 0), (0, 2)):
        C = F["C3"][ab]
        rC = [red(w, N, sl) for w in C]
        rdC = [red(w, N, sl) for w in F["dC3"][ab]]
        dimC = rank_of(rC)
        # Z = kernel of c -> red(dc) modulo kernel of c -> red(c)
        dimZ = dimC - rank_of(rdC)
        print(f"  ({ab[0]},{ab[1]}): dim C^3 = {dimC}, dim Z^3 = dim H^3 = {dimZ}")
        if ab == (1, 1):
            # which J's are closed / span
            Jn = [n for n in F["J"]]
            rJ = {n: red(F["J"][n], N, sl) for n in Jn}
            closed = {n: red(d(F["J"][n], sl), N, sl) == {} for n in Jn}
            # choose independent subset of J's
            E = Echelon(); basisJ = []
            for n in Jn:
                ok, _ = E.add(rJ[n])
                if ok:
                    basisJ.append(n)
            print(f"      J closed: {closed};  independent J's at N: {basisJ}; span dim {len(basisJ)} == dim Z^3: {len(basisJ) == dimZ}")
            H3basis = basisJ
    # H^4_(2,2)
    C4 = F["C4"]
    rC4 = [red(w, N, sl) for w in C4]
    rdC4 = [red(w, N, sl) for w in F["dC4"]]
    rB = [red(w, N, sl) for w in F["dC3"][(2, 2)]]
    dimC4 = rank_of(rC4)
    dimZ4 = dimC4 - rank_of(rdC4)
    dimB4 = rank_of(rB)
    print(f"  (2,2): dim C^4 = {dimC4}, dim Z^4 = {dimZ4}, dim B^4 = {dimB4}, dim H^4 = {dimZ4 - dimB4}")
    # explicit Z^4 vectors
    ker = kernel(rdC4)
    Zvec = []
    for rel in ker:
        v = {}
        for i, c in rel.items():
            for kk, x in rC4[i].items():
                add(v, kk, c * x)
        if v:
            Zvec.append(v)
    EB = Echelon()
    for i, b in enumerate(rB):
        EB.add(b, ('B', i))
    reps = []
    rF223 = red(F["F223"], N, sl)
    if EB.add(rF223, ('h', 'F223'))[0]:
        reps.append('F223')
    if not sl:
        F1 = {(s_, m_, 0): c_ for (s_, m_, n_), c_ in F["F223"].items() if n_ == 1}
        if EB.add(red(F1, N, sl), ('h', 'F1'))[0]:
            reps.append('F1')
    for i, z in enumerate(Zvec):
        if EB.add(z, ('h', f'z{i}'))[0]:
            reps.append(f'z{i}')
    print(f"  H^4 representatives: {reps} (count {len(reps)} == dim H^4: {len(reps) == dimZ4 - dimB4})")

    def cls(v):
        co = EB.coords(v)
        if co is None:
            return None
        return {nm2: c for (t, nm2), c in co.items() if t == 'h' and c != 0}

    print(f"  [F_223] != 0: {'F223' in reps};  [F_133 + F_223] = {cls(red(lin((1, F['F133']), (1, F['F223'])), N, sl))}")
    # quadratic map on H^3 basis
    if not H3basis:
        print("  H^3_(1,1) = 0: quadratic map trivial")
        continue
    q = {}
    for i, a in enumerate(H3basis):
        for j, b in enumerate(H3basis):
            if j < i:
                continue
            key = (a, b) if (a, b) in F["BR"] else (b, a)
            br = F["BR"][key]
            c = cls(red(br, N, sl))
            assert c is not None, "bracket not closed at N"
            fac = sympy.Rational(1, 2) if i == j else 1
            for rep, val in c.items():
                q[rep] = q.get(rep, 0) + fac * sympy.Rational(val.numerator, val.denominator) * XS[i] * XS[j]
    q = {k: sympy.factor(v) for k, v in q.items() if sympy.expand(v) != 0}
    print(f"  J = " + " + ".join(f"{XS[i]}*{n}" for i, n in enumerate(H3basis)))
    print(f"  (1/2){{J,J}} class: {q}")
    if not sl and N >= 3:
        lam = sympy.expand(q.get('F223', 0))
        ref = sympy.expand((2 * x0 * x2 + x1 ** 2 + N * x1 * x2) / 48)
        others = [k for k in q if k != 'F223']
        print(f"    F_223-coefficient == (2x0x2 + x1^2 + N x1x2)/48: {sympy.expand(lam - ref) == 0}")
        for k in others:
            print(f"    {k}-coefficient / ((2x0 + N x1)(2x1 + N x2)) = {sympy.simplify(q[k] / ((2 * x0 + N * x1) * (2 * x1 + N * x2)))}")
        print(f"    claimed mu = -(2x0 + N x1)(2x1 + N x2)/48 on [F_1]: {sympy.expand(q.get('F1', 0) + (2 * x0 + N * x1) * (2 * x1 + N * x2) / 48) == 0}; other components zero: {set(others) <= {'F1'}}")
        # (alpha, beta, gamma) basis
        al, be, ga = sympy.symbols('alpha beta gamma')
        sub = {x0: al, x1: -2 * al / N + be, x2: 4 * al / N ** 2 - 3 * be / N + ga}
        lam2 = sympy.expand(q.get('F223', 0).subs(sub, simultaneous=True))
        mu2 = sympy.expand(q.get('F1', 0).subs(sub, simultaneous=True))
        print(f"    in (alpha,beta,gamma): lam == (4a^2 - 2N^2 b^2 + N^3 b g)/(48N^2): {sympy.expand(lam2 - (4*al**2 - 2*N**2*be**2 + N**3*be*ga)/(48*N**2)) == 0};  mu == -N b (N g - b)/48: {sympy.expand(mu2 + N*be*(N*ga - be)/48) == 0}")
        # pi^* D_sl point
        sub = {x0: 1, x1: sympy.Rational(-2, N), x2: sympy.Rational(4, N * N)}
        print(f"    at J = pi^*D_sl = D - (2/N)S1 + (4/N^2)S0: " + str({k: sympy.nsimplify(v.subs(sub)) for k, v in q.items()}))
        sub = {x0: 0, x1: 1, x2: sympy.Rational(-2, N)}
        print(f"    at J = S1 - (2/N)S0: " + str({k: sympy.nsimplify(v.subs(sub)) for k, v in q.items()}))
        sub = {x0: sympy.Rational(1, 2), x1: 0, x2: 0}
        print(f"    at J = j_M = D/2: " + str({k: sympy.nsimplify(v.subs(sub)) for k, v in q.items()}))
    if sl and N >= 3:
        print(f"    F_223-coefficient == x0^2/(12N^2) = x0^2/{12 * N * N}: {sympy.expand(q.get('F223', 0) - x0 ** 2 / (12 * N * N)) == 0}")
    print(f"  time {time.time() - t0:.1f}s")
    sys.stdout.flush()
