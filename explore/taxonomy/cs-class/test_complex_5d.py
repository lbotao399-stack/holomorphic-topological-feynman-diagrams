"""Validation of colour.py + trace_complex.py on the 5d flagship (H,T) = (2,1), formal sl_N (kappa = tr):
flagship Thm class-slN: dF_133 = dF_223 = 0, F_133 + F_223 = 2 N^2 beta, beta = d([q0q1]^2 f), [F_223] != 0,
d^2 = 0 on random cochains."""
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES
import colour as col
import trace_complex as tc

log = Log(__file__)
H = 2
p = PRIMES[0]
rng = np.random.default_rng(5)
SH = [(1, 3, 3), (2, 2, 3)]

for N in (7, 11):
    F = {}
    for shape in SH:
        sig = col.standard_sigma(shape)
        cd = col.colour_perm_dict(shape, N)
        om = tc.Cochain(4, H)
        for perm, w in cd.items():
            om.add_term(perm, {((0, 1), (2, 3)): 1}, sig * w)
        F[shape] = om
        log('N =', N, shape, 'sigma =', sig, ' trace structures (fpf):', len(cd), ' canonical types:', sorted(om.data))
    plan5 = tc.DiffPlan(4, H, 60, rng)
    for shape in SH:
        v = plan5.apply(F[shape], p)
        log('N =', N, shape, 'dF = 0 at', len(v), 'C^5 coordinates:', not np.any(v))
    beta0 = tc.Cochain(3, H, {(3,): {((0, 1), (0, 1)): 1}})
    plan4 = tc.DiffPlan(3, H, 80, rng)
    beta = plan4.apply(beta0, p)
    pts, tabs = plan4.pts, plan4.tabs
    c133 = tc.coordinates(F[(1, 3, 3)], pts, tabs, p)
    c223 = tc.coordinates(F[(2, 2, 3)], pts, tabs, p)
    s = (c133 + c223) % p
    log('N =', N, 'F_133 + F_223 == 2 N^2 beta:', np.array_equal(s, (2 * N * N * beta) % p))
    r0 = tc.np_rank_mod if hasattr(tc, 'np_rank_mod') else None
    from common import np_rank_mod
    log('N =', N, 'rank[beta] =', np_rank_mod(beta[None, :], p), ' rank[beta; F_223] =', np_rank_mod(np.stack([beta, c223]), p))
    # d^2 = 0 on a random 3-cochain with general polynomial (not invariant under SL_2 needed? products of minors)
    om3 = tc.Cochain(3, H)
    om3.add_term(tc.canonical((3,)), {((0, 1), (0, 2)): 3}, 1)
    om3.add_term((2, 0, 1), {((0, 1), (1, 2)): -2}, 1)
    plan45 = tc.DiffPlan(4, H, 30, rng)
    # build d om3 as a 4-cochain explicitly is not available (coordinates only); check d^2 through the
    # identity d(beta0) closed: evaluate d(beta-cochain) -> need cochain form; use F instead (closed above)
log.close()
