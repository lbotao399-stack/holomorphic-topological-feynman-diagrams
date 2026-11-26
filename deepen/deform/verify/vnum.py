"""Explicit evaluation of cochains on matrices (exact integers/Fractions), explicit d and bracket
by direct structure-constant / dual-basis sums.  Independent of the formal contraction rules."""
import itertools, random
from fractions import Fraction as Fr
from vce import perm_cycles, perm_sign


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def comm(A, B):
    P = matmul(A, B); Q = matmul(B, A)
    n = len(A)
    return [[P[i][j] - Q[i][j] for j in range(n)] for i in range(n)]


def tau_eval(sig, X):
    v = 1
    for c in perm_cycles(sig):
        M = X[c[0]]
        for t in c[1:]:
            M = matmul(M, X[t])
        v *= tr(M)
    return v


def mono_eval(m, p):
    v = 1
    for s in range(len(p)):
        for comp in (0, 1):
            e = m[2 * s + comp]
            if e:
                v *= p[s][comp] ** e
    return v


def ev(w, X, p, N):
    tot = Fr(0)
    cache = {}
    for (sig, m, n), c in w.items():
        if sig not in cache:
            cache[sig] = tau_eval(sig, X)
        t = cache[sig]
        if t == 0:
            continue
        tot += c * Fr(N) ** n * mono_eval(m, p) * t
    return tot


def basis_gl(N):
    B = []
    for a in range(N):
        for b in range(N):
            E = [[0] * N for _ in range(N)]
            E[a][b] = 1
            B.append(E)
    return B


def basis_sl(N):
    B = []
    for a in range(N):
        for b in range(N):
            if a != b:
                E = [[0] * N for _ in range(N)]
                E[a][b] = 1
                B.append(E)
    for a in range(N - 1):
        E = [[0] * N for _ in range(N)]
        E[a][a] = 1
        E[a + 1][a + 1] = -1
        B.append(E)
    return B


def inv_gram(B):
    import sympy
    G = sympy.Matrix(len(B), len(B), lambda i, j: tr(matmul(B[i], B[j])))
    Gi = G.inv()
    return [[Fr(int(sympy.fraction(Gi[i, j])[0]), int(sympy.fraction(Gi[i, j])[1])) for j in range(len(B))]
            for i in range(len(B))]


def d_explicit(w, X, p, N):
    K = len(X)
    tot = Fr(0)
    for i in range(K):
        for j in range(i + 1, K):
            rest = [t for t in range(K) if t != i and t != j]
            X2 = [comm(X[i], X[j])] + [X[t] for t in rest]
            p2 = [(p[i][0] + p[j][0], p[i][1] + p[j][1])] + [p[t] for t in rest]
            tot += (-1) ** (i + j) * ev(w, X2, p2, N)
    return tot


def bracket_explicit(a, m, b, n, X, p, N, B, Gi):
    r = m + n - 2
    tot = Fr(0)
    for S in itertools.combinations(range(r), m - 1):
        Sc = [t for t in range(r) if t not in S]
        perm = list(S) + Sc
        sg = perm_sign(tuple(perm))
        q = (-sum(p[t][0] for t in S), -sum(p[t][1] for t in S))
        mq = (-q[0], -q[1])
        av = [ev(a, [X[t] for t in S] + [E], [p[t] for t in S] + [q], N) for E in B]
        bv = [ev(b, [E] + [X[t] for t in Sc], [mq] + [p[t] for t in Sc], N) for E in B]
        s = Fr(0)
        for e in range(len(B)):
            if av[e] == 0:
                continue
            for e2 in range(len(B)):
                if Gi[e][e2] != 0 and bv[e2] != 0:
                    s += av[e] * Gi[e][e2] * bv[e2]
        tot += sg * s
    return tot


def rand_point(N, k, sl, rng, R=4):
    X = []
    for _ in range(k):
        M = [[rng.randint(-R, R) for _ in range(N)] for _ in range(N)]
        if sl:
            t = sum(M[i][i] for i in range(N - 1))
            M[N - 1][N - 1] = -t
        X.append(M)
    p = [(rng.randint(-R, R), rng.randint(-R, R)) for _ in range(k - 1)]
    p.append((-sum(x[0] for x in p), -sum(x[1] for x in p)))
    return X, p
