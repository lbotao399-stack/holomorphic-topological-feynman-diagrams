"""C3 first product obstruction: third Koszul bracket of Q_2 on three gauge-invariant factors.

Phi3_{Q2}(tr b1 g^I, tr b2 g^J, tr b3 g^K) = Q2(product) [lower terms vanish: Q2 needs three contracted b/beta letters]
 = -i g^2 tau^3/(48 (32 pi^2)^2) sum_sigma (C_{x_s1} C_{x_s2} C_{x_s3})_{DE} B^{DE},  (C_x)_{DR} = x^A C_{ADR},
B^{DE} = (d_{+a}d_{+b} c^D)(d_+^a d_+^b c^E) = -B^{ED}   (LLC eq. 2-28 with spectators).
"""
from itertools import permutations, product
import numpy as np
import sympy as sp

# ---- antisymmetry and non-vanishing of B^{DE}: Grassmann algebra on u^D_{ab}, a,b in {+,-} symmetric ----
def grassmann_B():
    # generators u[D][k], k = 0:(++), 1:(+-), 2:(--); represent monomials as sorted tuples with sign
    def mul(x, y):
        out = {}
        for m1, c1 in x.items():
            for m2, c2 in y.items():
                if set(m1) & set(m2):
                    continue
                seq = list(m1 + m2)
                sign = 1
                for i in range(len(seq)):          # bubble sort, count transpositions
                    for j in range(len(seq) - 1 - i):
                        if seq[j] > seq[j + 1]:
                            seq[j], seq[j + 1] = seq[j + 1], seq[j]
                            sign = -sign
                key = tuple(seq)
                out[key] = out.get(key, 0) + sign * c1 * c2
        return {k: v for k, v in out.items() if v}

    def add(x, y, cy=1):
        out = dict(x)
        for k, v in y.items():
            out[k] = out.get(k, 0) + cy * v
        return {k: v for k, v in out.items() if v}

    idx = {('+', '+'): 0, ('+', '-'): 1, ('-', '+'): 1, ('-', '-'): 2}
    eps_up = {('+', '-'): 1, ('-', '+'): -1}
    u = lambda D, a, b: {((D, idx[(a, b)]),): 1}

    def B(D, E):
        tot = {}
        for a, b, c, d in product('+-', repeat=4):
            e = eps_up.get((a, c), 0) * eps_up.get((b, d), 0)
            if e:
                tot = add(tot, mul(u(D, a, b), u(E, c, d)), e)
        return tot
    b12, b21, b11 = B(1, 2), B(2, 1), B(1, 1)
    return add(b12, b21) == {}, b11 == {}, b12


anti, diag0, b12 = grassmann_B()
print('B^{DE} = -B^{ED}:', anti, ';  B^{DD} = 0:', diag0, ';  B^{12} =', b12)

# ---- SU(2), t = sigma/2: C_{ABC} = eps_{ABC}, tau = 1/2 ----
x = sp.symbols('x1:4')            # components of the spectator gamma (same flavour in all three factors)
y = sp.symbols('y1:4')
z = sp.symbols('z1:4')
Bs = {(1, 2): sp.Symbol('B12'), (1, 3): sp.Symbol('B13'), (2, 3): sp.Symbol('B23')}
eps = lambda i, j, k: sp.LeviCivita(i, j, k)


def Cx(v):
    return sp.Matrix(3, 3, lambda D, R: sum(v[A] * eps(A, D, R) for A in range(3)))


def T_su2(vs):
    T = sp.zeros(3, 3)
    for p in permutations(range(3)):
        T += Cx(vs[p[0]]) * Cx(vs[p[1]]) * Cx(vs[p[2]])
    return T


def contract(T):
    tot = 0
    for D in range(3):
        for E in range(3):
            if D == E:
                continue
            key = (min(D, E) + 1, max(D, E) + 1)
            tot += T[D, E] * Bs[key] * (1 if D < E else -1)
    return sp.factor(sp.expand(tot))


g, pi = sp.symbols('g pi', positive=True)
tau = sp.Rational(1, 2)
pref = -sp.I * g**2 * tau**3 / (48 * (32 * pi**2)**2)
same = contract(T_su2([x, x, x]))
print('SU(2), I=J=K:  sum_DE T_DE B^DE =', same)
print('   Phi3 =', sp.factor(pref * same))
mixed = contract(T_su2([x, y, z]))
print('SU(2), three different spectators: nonzero =', sp.expand(mixed) != 0)


# ---- SU(N) by matrices: Phi3 colour tensor = (i/tau) sum_sigma tr(ad_{x_s3} ad_{x_s2} ad_{x_s1}(t_D) t_E) ----
def large_N_trace(N, seed=1):
    rng = np.random.default_rng(seed)
    def tl():
        m = rng.integers(-2, 3, size=(N, N)).astype(object)
        m[N - 1, N - 1] = -sum(m[i, i] for i in range(N - 1))
        return sp.Matrix(m.tolist())
    X1, X2, X3 = tl(), tl(), tl()
    Dm, Em = tl(), tl()
    tot = 0
    for p in permutations((X1, X2, X3)):
        W = Dm
        for Xk in p:
            W = Xk * W - W * Xk
        tot += (W * Em).trace()
    return tot


for N in (2, 3, 4, 5):
    print(f'SU({N}) random traceless x1,x2,x3,D,E: sum_sigma tr(ad ad ad (D) E) =', large_N_trace(N))
