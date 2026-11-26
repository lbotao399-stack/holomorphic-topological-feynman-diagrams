"""Part (a) verification at random points, exact modulo primes (independent constructions vlie.py, engine vpoint.py).

For each algebra: two random points X = (X0..X3) (residues mod p), all 6 QB coefficients.
  (1) F133 + F223 - (1/2) beta = 0                                  (Killing form)
  (2) F223 - x beta - a G_kk = 0      with the explorer's (x, a)
  (3) minors beta[i] X[j] - X[i] beta[j] over the 12 values, X in {F223, F133, -3F133+2F223, -3F133-2F223, G_kk}:
      a nonzero minor mod p certifies X not in C beta (a rational identity X = c beta forces all minors to vanish).
Optional: d F_223 = 0 at a random point of 5 arguments (numeric momenta), for the keys given after '--dF'.
"""
import sys, time, itertools
import numpy as np
from fractions import Fraction as Fr
import vlie, vpoint

CLAIM = {  # (x, a) from simple-all REPORT 2.1 and logs/b_analysis.log
    'sl3': (Fr(-2, 3), Fr(-1, 12)), 'sl4': (Fr(-3, 8), Fr(-3, 64)), 'sl5': (Fr(-6, 25), Fr(-3, 100)),
    'sl6': (Fr(-1, 6), Fr(-1, 48)), 'sl7': (Fr(-6, 49), Fr(-3, 196)), 'so5': (Fr(-1), Fr(-1, 18)), 'sp4': (Fr(-1), Fr(-1, 18)),
    'so6': (Fr(-3, 8), Fr(-3, 64)),
    'so7': (Fr(-3, 25), Fr(-9, 250)), 'sp6': (Fr(-21, 32), Fr(-15, 512)), 'so8': (Fr(0), Fr(-1, 36)),
    'G2': (Fr(-7, 24), Fr(-5, 96)), 'F4': (Fr(4, 27), Fr(-5, 324)), 'so9': (Fr(3, 49), Fr(-15, 686)),
    'sp8': (Fr(-12, 25), Fr(-9, 500)), 'so10': (Fr(3, 32), Fr(-9, 512)), 'so11': (Fr(1, 9), Fr(-7, 486)),
    'sp10': (Fr(-3, 8), Fr(-7, 576)), 'so12': (Fr(3, 25), Fr(-3, 250)), 'E6': (Fr(5, 24), Fr(-1, 96)),
    'E7': (Fr(7, 27), Fr(-1, 162)), 'E8': (Fr(22, 75), Fr(-1, 300)),
}


def r(x, p):
    x = Fr(x)
    return x.numerator % p * pow(x.denominator, p - 2, p) % p


def first_minor(Bv, Xv, p):
    n = len(Bv)
    for i in range(n):
        if Bv[i] % p:
            for j in range(n):
                m = (Bv[i] * Xv[j] - Xv[i] * Bv[j]) % p
                if m:
                    return (i, j, m)
            return None
    return 'beta=0'


def run(key, p, npts=2, seed=11):
    t0 = time.time()
    L = vlie.build(key)
    E = vpoint.Engine(L, p)
    rng = np.random.default_rng(seed + p)
    vals = {'F133': [], 'F223': [], 'beta': [], 'Gkk': []}
    for _ in range(npts):
        X = [rng.integers(0, p, L.d).astype(np.float64) for _ in range(4)]
        F1, F2 = E.F(X)
        vals['F133'] += F1; vals['F223'] += F2
        vals['beta'] += E.beta(X); vals['Gkk'] += E.Gkk(X)
    h = r(Fr(1, 2), p)
    B = vals['beta']
    s1 = all((a + b - h * c) % p == 0 for a, b, c in zip(vals['F133'], vals['F223'], B))
    out = ['%-16s p=%d  d=%d  %d points x 6 coefficients' % (L.name, p, L.d, npts)]
    out.append('   F133 + F223 = (1/2) beta (Killing): %s' % s1)
    if key in CLAIM:
        x, a = CLAIM[key]
        s2 = all((f - r(x, p) * b - r(a, p) * g) % p == 0 for f, b, g in zip(vals['F223'], B, vals['Gkk']))
        out.append('   F223 = x beta + a G_kk with (x, a) = (%s, %s): %s' % (x, a, s2))
    if key == 'sl2':
        out.append('   F223 = -3 beta: %s ; G_kk = 8 beta: %s' % (all((f + 3 * b) % p == 0 for f, b in zip(vals['F223'], B)),
                                                         all((g - 8 * b) % p == 0 for g, b in zip(vals['Gkk'], B))))
    combos = {'F223': vals['F223'], 'F133': vals['F133'],
              '-3F133+2F223': [(-3 * a + 2 * b) % p for a, b in zip(vals['F133'], vals['F223'])],
              '-3F133-2F223': [(-3 * a - 2 * b) % p for a, b in zip(vals['F133'], vals['F223'])],
              'G_kk': vals['Gkk']}
    for name, Xv in combos.items():
        m = first_minor(B, Xv, p)
        out.append('   %-13s in C beta: %s' % (name, 'no (minor[%d,%d] = %d mod p)' % m if isinstance(m, tuple) else 'all minors 0 (%s)' % m))
    out.append('   time %.1fs' % (time.time() - t0))
    return '\n'.join(out)


def dF(key, p, seed=5):
    """d F_G at a random point (5 arguments, numeric momenta): sum_{i<j} (-1)^{i+j} F([x_i,x_j], x_k, x_l, x_m)."""
    t0 = time.time()
    L = vlie.build(key)
    E = vpoint.Engine(L, p)
    rng = np.random.default_rng(seed)
    X = [rng.integers(0, p, L.d).astype(np.float64) for _ in range(5)]
    P = [rng.integers(0, p, 2) for _ in range(4)]
    P.append(np.array([(-sum(q[0] for q in P)) % p, (-sum(q[1] for q in P)) % p]))
    Bk = lambda a, b: (int(a[0]) * int(b[1]) - int(a[1]) * int(b[0])) % p

    def Fnum(Y, Q):
        u, v, w = Bk(Q[0], Q[1]), Bk(Q[0], Q[2]), Bk(Q[1], Q[2])
        mono = [u * u, v * v, w * w, u * v, u * w, v * w]
        F1, F2 = E.F(Y)
        return (sum(m * f for m, f in zip(mono, F1)) % p, sum(m * f for m, f in zip(mono, F2)) % p)
    tot1 = tot2 = 0
    for i, j in itertools.combinations(range(5), 2):
        rest = [k for k in range(5) if k not in (i, j)]
        Y = [E.bracket(X[i], X[j])] + [X[k] for k in rest]
        Q = [(P[i] + P[j]) % p] + [P[k] for k in rest]
        f1, f2 = Fnum(Y, Q)
        sg = (-1) ** (i + j)
        tot1 = (tot1 + sg * f1) % p; tot2 = (tot2 + sg * f2) % p
    return '%-16s p=%d  d F_133 = %d, d F_223 = %d (mod p) at a random point  (%.1fs)' % (L.name, p, tot1, tot2, time.time() - t0)


if __name__ == '__main__':
    args = sys.argv[1:]
    p = 1048573
    if args and args[0].startswith('--p='):
        p = int(args[0][4:]); args = args[1:]
    if args and args[0] == '--dF':
        for key in args[1:]:
            print(dF(key, p), flush=True)
    else:
        for key in args:
            print(run(key, p), flush=True)
