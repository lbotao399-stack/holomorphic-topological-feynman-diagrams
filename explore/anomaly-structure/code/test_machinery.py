"""Tests of lie.py and cochain.py: Lie algebra data, d o d = 0, transposition conventions,
the Kac-Moody 2-cocycle (H = 1) is closed and not exact."""
import sys, random, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import lie, cochain as cc

rng = np.random.default_rng(5)
for g in [lie.sl(2), lie.gl(2), lie.sl(3), lie.gl(3)]:
    print(g.name, 'dim', g.dim, 'centre dim', g.center_dim(), 'kinv scale', g.kinv_scale)

# transposition convention
W = rng.integers(-5, 5, size=(2, 3, 4, 5))
sig = (2, 0, 3, 1)
B = np.transpose(W, cc.transpose_for(sig))
ok = all(B[a] == W[tuple(a[s] for s in sig)] for a in itertools.product(*[range(n) for n in B.shape]))
print('transpose convention', ok, B.shape)

# d o d = 0 on random cochains
for g, H, k, degs in [(lie.sl(2), 1, 2, (2,)), (lie.sl(2), 2, 2, (1, 1)), (lie.gl(2), 2, 2, (1, 1)), (lie.sl(3), 1, 2, (1,))]:
    mom = cc.Momenta(6, H)
    T = rng.integers(-3, 4, size=(g.dim,) * k)
    P = sum(int(rng.integers(-3, 4)) * mom.mono_expr(m, k + 1, [mom.p[i] for i in range(k)])
            for m in cc.all_monos(mom, k + 1, degs))
    w = cc.alternate(mom, P, T, k)
    dw = cc.d_CE(mom, g, w, k)
    ddw = cc.d_CE(mom, g, dw, k + 1)
    print('d d = 0:', g.name, 'H', H, 'k', k, cc.is_zero(ddw), 'dw nonzero', not cc.is_zero(dw))

# alternating property: alternate(alternate) = k! alternate  (P depends on all k momenta)
g = lie.sl(2); mom = cc.Momenta(6, 2); k = 3
T = rng.integers(-3, 4, size=(g.dim,) * k)
P = mom.p[0][0] * mom.p[1][1] + 2 * mom.p[2][0] * mom.p[0][1]
w = cc.alternate(mom, P, T, k)
# rebuild from canonical form as a P with all momenta and re-alternate
monos = sorted(w)
ok = True
for m in monos:
    e = mom.mono_expr(m, k, [mom.p[i] for i in range(k - 1)])
    w2 = cc.alternate(mom, e, w[m], k)
    # w2 should equal 6 * (contribution of the single term m (x) w[m]) alternated: sum over m of w2 = 6 w
w_sum = {}
for m in monos:
    e = mom.mono_expr(m, k, [mom.p[i] for i in range(k - 1)])
    for M, Tm in cc.alternate(mom, e, w[m], k).items():
        w_sum[M] = w_sum.get(M, 0) + Tm
print('alternate(alternate) = 3! alternate:', all(np.array_equal(w_sum.get(M, 0) * 1, 6 * w.get(M, 0)) for M in set(w) | set(w_sum)))

# Kac-Moody cocycle for H = 1: omega(X e^{pz}, Y e^{qz}) = p kappa(X, Y): closed, not exact
for g in [lie.sl(2), lie.sl(3)]:
    mom = cc.Momenta(5, 1)
    w = cc.alternate(mom, mom.p[0][0], g.kappa, 2)   # alternation doubles it
    print('KM cocycle', g.name, 'closed', cc.is_zero(cc.d_CE(mom, g, w, 2)))
    # exactness: d of 1-cochains of degree 1: eta(X e^{pz}) = p phi(X); canonical form of a 1-cochain
    # at k = 1 has no free momentum (p_1 = 0): image of d from degree-1 one-cochains is zero.
    mono1 = cc.all_monos(mom, 1, (1,))
    print('  canonical degree-1 monomials of 1-cochains:', mono1)
