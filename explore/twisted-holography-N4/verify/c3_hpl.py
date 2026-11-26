"""C3 toy check of Q(1+Delta s) = (1+Delta s) Q0 on random filtered complexes over Q (truncated at hbar^K)."""
import random
import sympy as sp

random.seed(5)
h, a, K = 2, 3, 3
n = h + 2 * a


def rmat(r, c, lo=-3, hi=3):
    return sp.Matrix(r, c, lambda i, j: random.randint(lo, hi))


# standard form V = H + A + B, Q0: B -> A (invertible M), s: A -> B
while True:
    M = rmat(a, a)
    if M.det() != 0:
        break
Z = sp.zeros
Q0s = sp.BlockMatrix([[Z(h, h), Z(h, a), Z(h, a)], [Z(a, h), Z(a, a), M], [Z(a, h), Z(a, a), Z(a, a)]]).as_explicit()
ss = sp.BlockMatrix([[Z(h, h), Z(h, a), Z(h, a)], [Z(a, h), Z(a, a), Z(a, a)], [Z(a, h), M.inv(), Z(a, a)]]).as_explicit()
iota_s = sp.BlockMatrix([[sp.eye(h)], [Z(a, h)], [Z(a, h)]]).as_explicit()
pi_s = iota_s.T
while True:
    P = rmat(n, n)
    if P.det() != 0:
        break
Pi = P.inv()
Q0, s, iota, pi = P * Q0s * Pi, P * ss * Pi, P * iota_s, pi_s * Pi
assert Q0 * Q0 == Z(n, n)
assert Q0 * s + s * Q0 == sp.eye(n) - iota * pi
assert pi * iota == sp.eye(h)


def solve_order(R, require_iota):
    """random X with Q0 X + X Q0 = R (and X iota = 0 if required)."""
    xs = sp.symbols(f'x0:{n*n}')
    X = sp.Matrix(n, n, xs)
    eqs = list(Q0 * X + X * Q0 - R)
    if require_iota:
        eqs += list(X * iota)
    sol = sp.linsolve(eqs, xs)
    (sol,) = sol
    free = sorted(set().union(*[e.free_symbols for e in sol]), key=str)
    sub = {f: random.randint(-2, 2) for f in free}
    return sp.Matrix(n, n, [e.subs(sub) for e in sol])


def build(require_iota):
    Qs = [Q0]
    for k in range(1, K + 1):
        R = -sum((Qs[i] * Qs[k - i] for i in range(1, k)), Z(n, n))
        Qs.append(solve_order(R, require_iota))
    return Qs


def check(Qs, label):
    t = sp.symbols('hbar')
    Q = sum((Qs[k] * t**k for k in range(K + 1)), Z(n, n))
    D = Q - Q0
    trunc = lambda Mx: Mx.applyfunc(lambda e: sp.expand(sum(sp.expand(e).coeff(t, k) * t**k for k in range(K + 1))))
    sq = trunc(Q * Q)
    lhs = trunc(Q * (sp.eye(n) + D * s) - (sp.eye(n) + D * s) * Q0)
    two_loop = Qs[2] - ((Qs[2] * s) * Q0 - Q0 * (Qs[2] * s) - Qs[1] * Qs[1] * s)
    print(f'[{label}] Q^2 = 0 mod hbar^{K+1}: {sq == Z(n, n)};  Delta iota = 0: {trunc(D * iota) == Z(n, h)};'
          f'  Q(1+Ds)-(1+Ds)Q0 = 0: {lhs == Z(n, n)};  Q2 = [Q2 s,Q0] - Q1^2 s: {two_loop == Z(n, n)}')
    if lhs == Z(n, n):
        g = sp.eye(n) + D * s
        # inverse as a truncated series: g = 1 + O(hbar)
        X = D * s
        ginv = sum(((-X)**k for k in range(1, K + 1)), sp.eye(n))
        conj = trunc(trunc(ginv * Q) * g)
        print(f'   g^-1 Q g = Q0 mod hbar^{K+1}: {conj == Q0}')
    # induced differential on H: pi Delta iota + pi Delta s Delta iota + ...
    return lhs == Z(n, n)


ok1 = check(build(True), 'random Q with Q iota = 0')
# negative control: Q = g Q0 g^-1 with g acting nontrivially on H (Q^2 = 0 but Delta iota != 0)
t = sp.symbols('hbar')
X1 = rmat(n, n)
g = sp.eye(n) + t * X1
ginv = sum(((-t * X1)**k for k in range(1, K + 1)), sp.eye(n))
Qc = (g * Q0 * ginv).applyfunc(sp.expand)
Qs = [Qc.applyfunc(lambda e: e.coeff(t, k)) for k in range(K + 1)]
Qs[0] = Q0
ok2 = check(Qs, 'control: conjugate by g with g iota != iota')
