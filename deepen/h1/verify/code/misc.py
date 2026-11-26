"""(1) abelian C_g^0 dimensions by subset counting + Koszul and by Molien; (2) heat-kernel operator identities (sympy);
(3) character enumeration for claim 'chars'."""
import sys, itertools
from fractions import Fraction
from math import factorial
import sympy as sp

# ---------------------------------------------------------------- (1) abelian counts
def subset_count(d, a, b, n):
    """number of n-subsets of {(color, al1, al2)} with sum al = (a,b)  (dim of C_g^n at bidegree (a,b), d colours)."""
    gens = [(c, i, j) for c in range(d) for i in range(a + 1) for j in range(b + 1)]
    # dp over generators: poly in (deg, A, B)
    dp = {(0, 0, 0): 1}
    for (_, i, j) in gens:
        nd = dict(dp)
        for (k, A, B), v in dp.items():
            if k < n and A + i <= a and B + j <= b:
                key = (k + 1, A + i, B + j)
                nd[key] = nd.get(key, 0) + v
        dp = nd
    return dp.get((n, a, b), 0)


def koszul_dim(d, a, b, n):
    f = lambda x, y: subset_count(d, x, y, n) if x >= 0 and y >= 0 else 0
    return f(a, b) - f(a - 1, b) - f(a, b - 1) + f(a - 1, b - 1)


def molien_dim(d, a, b, n):
    """dim of (C[p_1..p_n]_{(a,b)} (x) (g^v)^{(x)n} (x) sgn)^{S_n}, restricted to sum p = 0 (times (1-s)(1-t))."""
    s, t = sp.symbols("s t")
    tot = 0
    for perm in itertools.permutations(range(n)):
        # cycle type
        seen = [False] * n
        cyc = []
        for i in range(n):
            if not seen[i]:
                l = 0
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = perm[j]
                    l += 1
                cyc.append(l)
        sgn = (-1) ** sum(l - 1 for l in cyc)
        term = sp.Integer(sgn) * d ** len(cyc)
        for l in cyc:
            term = term / ((1 - s ** l) * (1 - t ** l))
        tot += term
    gf = sp.series(sp.series(tot * (1 - s) * (1 - t) / factorial(n), s, 0, a + 1).removeO(), t, 0, b + 1).removeO()
    return sp.Poly(sp.expand(gf), s, t).coeff_monomial(s ** a * t ** b)


print("(1) abelian C_g^0 dims: subset-count+Koszul vs Molien")
for d in (1, 2):
    for (a, b) in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (2, 2), (3, 1)]:
        row = []
        for n in range(1, 7):
            k1 = koszul_dim(d, a, b, n)
            k2 = molien_dim(d, a, b, n) if n <= 5 else None
            row.append((n, k1, k2))
        print(f"  ab{d} (a,b)={(a, b)}: " + " ".join(f"n={n}:{k1}/{k2}" for n, k1, k2 in row if k1 or k2))

# ---------------------------------------------------------------- (2) heat-kernel operator identities
print("(2) operator identities on Omega*(R) (x) Omega^{0,*}(C^2), generic coefficients")
tau, z1, z2, w1, w2 = sp.symbols("tau z1 z2 zb1 zb2")
X = [tau, w1, w2]  # derivative variable for form direction v
Zh = [z1, z2]


def nb(T, v):
    return bin(T & ((1 << v) - 1)).count("1")


def iota(F, v):
    R = {}
    for T, f in F.items():
        if (T >> v) & 1:
            s = -1 if nb(T, v) % 2 else 1
            R[T & ~(1 << v)] = R.get(T & ~(1 << v), 0) + s * f
    return R


def wedge(v, F):
    R = {}
    for T, f in F.items():
        if not (T >> v) & 1:
            s = -1 if nb(T, v) % 2 else 1
            R[T | (1 << v)] = R.get(T | (1 << v), 0) + s * f
    return R


def dif(F, var):
    return {T: sp.diff(f, var) for T, f in F.items()}


def plus(*Fs, c=None):
    R = {}
    for i, F in enumerate(Fs):
        k = 1 if c is None else c[i]
        for T, f in F.items():
            R[T] = R.get(T, 0) + k * f
    return R


def Dop(F):
    return plus(*[wedge(v, dif(F, X[v])) for v in range(3)])


def Dstar(F):
    return plus(dif(iota(F, 0), tau), dif(iota(F, 1), z1), dif(iota(F, 2), z2), c=[-1, -4, -4])


def iszero(F):
    return all(sp.simplify(f) == 0 for f in F.values())


fs = {T: sp.Function(f"f{T}")(tau, z1, z2, w1, w2) for T in range(8)}
F = fs
ok = {}
for v in range(3):
    ok[f"[D*,iota_{v}]=0"] = iszero(plus(Dstar(iota(F, v)), iota(Dstar(F), v)))
    ok[f"[D,iota_{v}]=d_{v}"] = iszero(plus(Dop(iota(F, v)), iota(Dop(F), v), dif(F, X[v]), c=[1, 1, -1]))
H = plus(Dop(Dstar(F)), Dstar(Dop(F)))
lap = {T: -(sp.diff(f, tau, 2) + 4 * sp.diff(f, z1, w1) + 4 * sp.diff(f, z2, w2)) for T, f in F.items()}
ok["[D,D*]=-(d_tau^2+4 sum d_z d_zb)"] = iszero(plus(H, lap, c=[1, -1]))
for v in range(3):
    HF = lambda G: plus(Dop(Dstar(G)), Dstar(Dop(G)))
    ok[f"[[D,D*],iota_{v}]=0"] = iszero(plus(HF(iota(F, v)), iota(HF(F), v), c=[1, -1]))
for k, v in ok.items():
    print(f"  {k}: {v}")

# X x X: kernel k(x-y) prod_v (dx_v - dy_v); bits 0..2 dx, 3..5 dy
print("(2b) X x X: (iota^x_v + iota^y_v) annihilates K = k(x-y) (dx-dy)^3 and D*_x K")
xs = sp.symbols("tx z1x z2x w1x w2x")
ys = sp.symbols("ty z1y z2y w1y w2y")
kf = sp.Function("k")(*[a - b for a, b in zip(xs, ys)])


def nb6(T, v):
    return bin(T & ((1 << v) - 1)).count("1")


def wedge6(v, F):
    R = {}
    for T, f in F.items():
        if not (T >> v) & 1:
            s = -1 if nb6(T, v) % 2 else 1
            R[T | (1 << v)] = R.get(T | (1 << v), 0) + s * f
    return R


def iota6(F, v):
    R = {}
    for T, f in F.items():
        if (T >> v) & 1:
            s = -1 if nb6(T, v) % 2 else 1
            R[T & ~(1 << v)] = R.get(T & ~(1 << v), 0) + s * f
    return R


def plus6(*Fs, c=None):
    R = {}
    for i, Fm in enumerate(Fs):
        k = 1 if c is None else c[i]
        for T, f in Fm.items():
            R[T] = R.get(T, 0) + k * f
    return R


# build prod_v (dx_v - dy_v) left to right
K = {0: kf}
for v in range(3):
    K = plus6(wedge6(v, K), wedge6(v + 3, K), c=[1, -1])
    # wedge on the left; order (dx_2-dy_2)(dx_1-dy_1)(dx_0-dy_0) k: any fixed order suffices
xvars = [xs[0], xs[3], xs[4]]
Dsx = plus6({T: sp.diff(f, xs[0]) for T, f in iota6(K, 0).items()},
            {T: sp.diff(f, xs[1]) for T, f in iota6(K, 1).items()},
            {T: sp.diff(f, xs[2]) for T, f in iota6(K, 2).items()}, c=[-1, -4, -4])
for name, G in [("K", K), ("D*_x K", Dsx)]:
    for v in range(3):
        R = plus6(iota6(G, v), iota6(G, v + 3))
        print(f"  (iota^x_{v} + iota^y_{v}) {name} = 0: {all(sp.simplify(f) == 0 for f in R.values())}")
        R2 = plus6(iota6(G, v), iota6(G, v + 3), c=[1, -1])
        print(f"  control (iota^x_{v} - iota^y_{v}) {name} = 0: {all(sp.simplify(f) == 0 for f in R2.values())}")

# ---------------------------------------------------------------- (3) characters
print("(3) characters of psi_S(phi), phi in C_g^0_(a,b), w = -1_{U-S}; I_k: U(1)^2 (l1 l2)^{k-1}, dilation lam^{2k-2}")
sols_u = []
sols_ud = []
for k in range(0, 6):
    for S in range(8):
        w = [0 if (S >> v) & 1 else -1 for v in range(3)]
        for a in range(0, 2 * k + 4):
            for b in range(0, 2 * k + 4):
                u1 = (a - 1 - w[1] == k - 1) and (b - 1 - w[2] == k - 1)
                dil = (w[0] + a - 1 + b - 1 + w[1] + w[2] == 2 * k - 2)
                if u1:
                    sols_u.append((k, S, a, b, dil))
for k, S, a, b, dil in sols_u:
    if S != 7 or (a, b) != (k, k) or not dil:
        print(f"  k={k} S={S:03b} (bits tau,1b,2b = 1,2,4) (a,b)=({a},{b}) passes U(1)^2; passes dilation: {dil}")
print("  all U(1)^2+dilation solutions:", sorted(set((k, S, a, b) for k, S, a, b, dil in sols_u if dil)))
