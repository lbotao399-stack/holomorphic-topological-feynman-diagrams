"""All N >= 9 (formal sl_N, kappa = tr): B = d(C^7_inv) does not depend on N.  With random functionals l_1..l_4
vanishing on B (mod p), the class matrix M(N) = [l_i(o'), l_i(Q55), l_i({j5,j23}), l_i(Q2323)] has entries
polynomial in N (o': degree <= 2) resp. in 1/N (Q: Q0 - Q1/N).  f(N) = N^3 det M(N) is a polynomial of degree <= 5
(interpolated at 6 values, verified at the others).  f(N) != 0 mod p  =>  [o'] not in span{[Q]} at N;
g(N) = l_1(o'(N)) (degree <= 2) != 0  =>  [o'] != 0.  Roots of f, g in F_p via gcd with x^p - x."""
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement
from common import Log, PRIMES
import colour as col
import trace_complex as tc
import cartan_lib as cl

log = Log(__file__)
H = 4
rng = np.random.default_rng(91)
plan8 = tc.DiffPlan(7, H, 90, rng)
subs = list(combinations(range(6), 4))
gens = []
for lam in tc.fpf_types(7):
    for I, J in combinations_with_replacement(subs, 2):
        g = tc.Cochain(7, H)
        g.add_term(tc.canonical(lam), {(I, J): 1}, 1)
        if g.data:
            gens.append(g)
T5, T23 = tc.sym_trace_dict((5,)), tc.sym_trace_dict((2, 3))
Tmix = {k: T5.get(k, 0) + T23.get(k, 0) for k in set(T5) | set(T23)}


def matmul_mod(A, B, p):
    A = np.asarray(A, dtype=np.int64) % p
    B = np.asarray(B, dtype=np.int64) % p
    lo, hi = B & 0xFFFF, B >> 16
    return ((A @ lo) % p + ((A @ hi) % p) * 65536) % p


def nullspace_mod(B, p):
    M = np.array(B, dtype=np.int64) % p
    nr, nc = M.shape
    piv = []
    r = 0
    for c in range(nc):
        if r == nr:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        M[[r, i]] = M[[i, r]]
        M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        f = M[:, c].copy(); f[r] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(f[nzr], M[r]) % p) % p
        piv.append(c); r += 1
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fc in free:
        v = np.zeros(nc, dtype=np.int64); v[fc] = 1
        for i, c in enumerate(piv):
            v[c] = (-M[i, fc]) % p
        basis.append(v)
    return np.array(basis)


# ---- polynomials mod p (lists of coefficients, low degree first)
def pmod(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def pdivmod(a, b, p):
    a = pmod(a, p); b = pmod(b, p)
    q = [0] * max(len(a) - len(b) + 1, 1)
    inv = pow(b[-1], p - 2, p)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        k = len(a) - len(b)
        q[k] = c
        for i, x in enumerate(b):
            a[i + k] = (a[i + k] - c * x) % p
        a = pmod(a, p)
    return pmod(q, p), a


def pmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return pmod(out, p)


def pgcd(a, b, p):
    a, b = pmod(a, p), pmod(b, p)
    while b:
        a, b = b, pdivmod(a, b, p)[1]
    return a


def n_roots_mod(f, p):
    """number of distinct roots of f in F_p: deg gcd(f, x^p - x)."""
    f = pmod(f, p)
    if len(f) <= 1:
        return 0 if f else None
    # x^p mod f by square and multiply
    res, base, e = [1], [0, 1], p
    while e:
        if e & 1:
            res = pdivmod(pmul(res, base, p), f, p)[1]
        base = pdivmod(pmul(base, base, p), f, p)[1]
        e >>= 1
    xpx = pmod([(res[i] if i < len(res) else 0) - (1 if i == 1 else 0) for i in range(max(len(res), 2))], p)
    g = pgcd(f, xpx, p)
    return len(g) - 1


def interp(xs, ys, p):
    """Lagrange interpolation mod p -> coefficient list."""
    n = len(xs)
    coef = [0] * n
    for i in range(n):
        num, den = [1], 1
        for j in range(n):
            if j != i:
                num = pmul(num, [(-xs[j]) % p, 1], p)
                den = den * (xs[i] - xs[j]) % p
        c = ys[i] * pow(den, p - 2, p) % p
        for k, x in enumerate(num):
            coef[k] = (coef[k] + c * x) % p
    return pmod(coef, p)


def peval(c, x, p):
    return sum(ci * pow(x, i, p) for i, ci in enumerate(c)) % p


NVALS = list(range(9, 23))
for p in PRIMES[:2]:
    B = np.stack([plan8.apply(g, p) for g in gens])
    K = nullspace_mod(B, p)
    assert not np.any(matmul_mod(B, K.T, p))
    R = rng.integers(0, p, size=(4, K.shape[0]))
    ell = np.array([[sum(int(R[i, k]) * int(K[k, c]) for k in range(K.shape[0])) % p for c in range(K.shape[1])] for i in range(4)], dtype=np.int64)
    fv, gv = [], []
    for N in NVALS:
        o = None
        for shape in cl.SHAPES:
            om = tc.Cochain(8, H)
            for perm, w in col.colour_perm_dict(shape, N).items():
                om.add_term(perm, {(I, J): e for I, J, e in cl.DTH[shape]}, cl.SIGMA * w)
            om = om.scaled(cl.KC[shape] / cl.AUT[shape])
            o = om if o is None else o.plus(om)
        vo = tc.coordinates(o, plan8.pts, plan8.tabs, p)
        q55 = tc.coordinates(tc.half_bracket_sym(T5, H, 120, N), plan8.pts, plan8.tabs, p)
        q2323 = tc.coordinates(tc.half_bracket_sym(T23, H, 120, N), plan8.pts, plan8.tabs, p)
        qmix = tc.coordinates(tc.half_bracket_sym(Tmix, H, 120, N), plan8.pts, plan8.tabs, p)
        q523 = (qmix - q55 - q2323) % p
        cols = [vo, q55 * N % p, q523 * N % p, q2323 * N % p]
        M = [[int(sum(int(a) * int(b) for a, b in zip(ell[i], c)) % p) for c in cols] for i in range(4)]
        import sympy as sp
        d = int(sp.Matrix(M).det() % p)
        fv.append(d)
        gv.append(M[0][0])
    fpoly = interp(NVALS[:6], fv[:6], p)
    gpoly = interp(NVALS[:3], gv[:3], p)
    okf = all(peval(fpoly, N, p) == v for N, v in zip(NVALS, fv))
    okg = all(peval(gpoly, N, p) == v for N, v in zip(NVALS, gv))
    log('p = %d: f(N) = N^3 det M(N): degree %d (interpolated at N=9..14, verified at N=9..22: %s), roots in F_p: %s'
        % (p, len(fpoly) - 1, okf, n_roots_mod(fpoly, p)))
    log('p = %d: g(N) = l_1(o\'(N)): degree %d (verified at N=9..22: %s), roots in F_p: %s; f(9..22) all nonzero: %s'
        % (p, len(gpoly) - 1, okg, n_roots_mod(gpoly, p), all(fv)))
log.close()
