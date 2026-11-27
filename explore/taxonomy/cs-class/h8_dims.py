"""Dimensions of the SL_H x SL_N-invariant translation-coinvariant CE complex of sl_N[[z_1..z_H]] (formal N, trace
basis) at multidegree (2,..,2): dim C^k (rank of spanning generators at coordinates), rank d_k, dim H^k.
H = 2 (flagship weight (2,2), k = 3,4,5) as a control; H = 4 (9d, k = 6,7,8)."""
import sys
import numpy as np
from itertools import combinations, combinations_with_replacement
from common import Log, PRIMES, np_rank_mod
import trace_complex as tc

log = Log(__file__)
p = PRIMES[0]


def generators(k, H):
    subs = list(combinations(range(k - 1), H))
    gens = []
    for lam in tc.fpf_types(k):
        for I, J in combinations_with_replacement(subs, 2):
            g = tc.Cochain(k, H)
            g.add_term(tc.canonical(lam), {(I, J): 1}, 1)
            if g.data:
                gens.append(g)
    return gens


def dims(H, ks, M):
    rng = np.random.default_rng(4)
    res = {}
    for k in ks:
        gens = generators(k, H)
        pl_here = tc.DiffPlan(k - 1, H, M, rng)          # coordinates of C^k
        pl_next = tc.DiffPlan(k, H, M, rng)              # d : C^k -> C^{k+1}
        Cm = np.stack([tc.coordinates(g, pl_here.pts, pl_here.tabs, p) for g in gens])
        dC = np_rank_mod(Cm, p)
        Dm = np.stack([pl_next.apply(g, p) for g in gens])
        rd = np_rank_mod(Dm, p)
        res[k] = (len(gens), dC, rd)
        log('H = %d, k = %d: generators %d, dim C^k = %d, rank d_k = %d (coordinates %d / %d)'
            % (H, k, len(gens), dC, rd, Cm.shape[1], Dm.shape[1]))
    for k in ks[1:-1]:
        if k - 1 in res and k in res:
            dC, rd = res[k][1], res[k][2]
            log('H = %d: dim Z^%d = %d, dim B^%d = %d, dim H^%d = %d' % (H, k, dC - rd, k, res[k - 1][2], k, dC - rd - res[k - 1][2]))
    return res


which = sys.argv[1:] or ['2', '4']
if '2' in which:
    dims(2, [3, 4, 5], 120)
if '4' in which:
    dims(4, [6, 7, 8], 160)
log.close()
