"""Validation of trace_complex.half_bracket_sym on the 5d flagship Thm quad (formal sl_N, kappa = tr):
1/2 {D,D}_sl = F_223/(12 N^2) + d eta_N, eta_N = (N^2+4)/(2N^2) [p1p2]^2 tr(1[23]), D = [p1p2] tr(1{23})."""
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, np_rank_mod
import colour as col
import trace_complex as tc

log = Log(__file__)
H = 2
rng = np.random.default_rng(9)
for N in (5, 7, 11):
    for p in PRIMES[:2]:
        F223 = tc.Cochain(4, H)
        for perm, w in col.colour_perm_dict((2, 2, 3), N).items():
            F223.add_term(perm, {((0, 1), (2, 3)): 1}, -w)
        Tdict = {(1, 2, 0): 1, (2, 0, 1): 1}
        DD = tc.half_bracket_sym(Tdict, H, 1, N)
        plan4 = tc.DiffPlan(3, H, 80, rng)
        beta = plan4.apply(tc.Cochain(3, H, {(3,): {((0, 1), (0, 1)): 1}}), p)
        cDD = tc.coordinates(DD, plan4.pts, plan4.tabs, p)
        cF = tc.coordinates(F223, plan4.pts, plan4.tabs, p)
        a = Fr(1, 12 * N * N); b = Fr(N * N + 4, 2 * N * N)
        am = a.numerator * pow(a.denominator, p - 2, p) % p
        bm = b.numerator * pow(b.denominator, p - 2, p) % p
        res = (cDD - am * cF - bm * beta) % p
        plan5 = tc.DiffPlan(4, H, 40, rng)
        dDD = plan5.apply(DD, p)
        log('N =', N, 'p =', p, ' 1/2{D,D} - F_223/(12N^2) - (N^2+4)/(2N^2) beta = 0:', not np.any(res),
            ' d(1/2{D,D}) = 0:', not np.any(dDD))
log.close()
