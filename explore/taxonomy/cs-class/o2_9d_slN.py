"""9d minimal twist Map(C^4 x R_dR, Bg), (H,T) = (4,1), two-loop obstruction for formal sl_N (kappa = tr),
trace basis (faithful for N >= 9).  Cochains of sl_N[[z_1..z_4]] restricted to sum p = 0, SL_4-invariant
(products of two 4-minors), multidegree (2,2,2,2).

  F_G  = alt(D_theta,G  W^G),  W^G = sigma_G C_std,G              (G = theta(1,5,5),(2,4,5),(3,3,5),(3,4,4))
  o'   = sum_G |Aut G|^{-1} alt((K_G/pi) W^G),  K_G/pi = s_G c_1,G D_theta,G
  o_2  = U'/(1024 pi^8) o'            (G = c_{4,1}^{-2} K up to the sign U', c_{4,1} = 2 v_4 pi^{9/2}, v_4 = 16)
  j_T  = alt([p_1..p_4] T) = 120 [p_1..p_4] T,  T in (Sym^5 g^*)^g:  T_5 = Sym tr X^5,  T_23 = Sym tr X^2 tr X^3
Checks: d F_G = 0 (C^9 coordinates), d j_T = 0 (C^6), d(1/2{j,j}) = 0; B^8 = d(C^7_inv); classes."""
import sys
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, np_rank_mod, np_solve_mod, ratrec, crt_pair
import colour as col
import trace_complex as tc

log = Log(__file__)
H = 4
SHAPES = [(1, 5, 5), (2, 4, 5), (3, 3, 5), (3, 4, 4)]
KC = {(1, 5, 5): Fr(1, 5184), (2, 4, 5): Fr(-1, 9072), (3, 3, 5): Fr(1, 10080), (3, 4, 4): Fr(-1, 22680)}   # s c_1
AUT = {(1, 5, 5): 4, (2, 4, 5): 2, (3, 3, 5): 4, (3, 4, 4): 4}
DTH = {(1, 5, 5): {((0, 1, 2, 3), (4, 5, 6, 7)): 1}, (2, 4, 5): {((0, 1, 2, 3), (4, 5, 6, 7)): 1},
       (3, 3, 5): {((0, 1, 2, 3), (4, 5, 6, 7)): 1},
       (3, 4, 4): {((0, 1, 2, 3), (4, 5, 6, 7)): 1, ((0, 1, 2, 4), (3, 5, 6, 7)): -1, ((0, 1, 3, 4), (2, 5, 6, 7)): 1}}
NS = [int(x) for x in sys.argv[1:]] or [9, 10, 11, 12, 13, 14, 15]
M7, M8, M9 = 70, 110, 30
rng = np.random.default_rng(2026)

plan6 = tc.DiffPlan(5, H, 40, rng)            # d : C^5 -> C^6
plan5 = tc.DiffPlan(4, H, 40, rng)            # coordinates of C^5
plan6_src_pts, plan6_src_tabs = plan5.pts, plan5.tabs
plan8 = tc.DiffPlan(7, H, M8, rng)            # d : C^7 -> C^8 (coordinates of C^8)
plan9 = tc.DiffPlan(8, H, M9, rng)            # d : C^8 -> C^9
# C^7 generators: alt([I][J] T_{pi_lambda}), I, J 4-subsets of {0..5}
from itertools import combinations, combinations_with_replacement
subs = list(combinations(range(6), 4))
gens = []
for lam in tc.fpf_types(7):
    for I, J in combinations_with_replacement(subs, 2):
        om = tc.Cochain(7, H)
        om.add_term(tc.canonical(lam), {(I, J): 1}, 1)
        if om.data:
            gens.append(om)
log('C^7 generators (spanning):', len(gens), ' C^8 coordinates:', sum(v.shape[0] for v in plan8.pts.values()),
    ' C^9 coordinates:', sum(v.shape[0] for v in plan9.pts.values()))
T5 = tc.sym_trace_dict((5,))
T23 = tc.sym_trace_dict((2, 3))
results = {}
for N in NS:
    out = {}
    for p in PRIMES[:2]:
        B = np.stack([plan8.apply(g, p) for g in gens])
        rB = np_rank_mod(B, p)
        F = {}
        for shape in SHAPES:
            sig = col.standard_sigma(shape)
            om = tc.Cochain(8, H)
            for perm, w in col.colour_perm_dict(shape, N).items():
                om.add_term(perm, DTH[shape], sig * w)
            F[shape] = om
        o = F[SHAPES[0]].scaled(KC[SHAPES[0]] / AUT[SHAPES[0]])
        for shape in SHAPES[1:]:
            o = o.plus(F[shape], KC[shape] / AUT[shape])
        closed = {s: not np.any(plan9.apply(F[s], p)) for s in SHAPES}
        cF = {s: tc.coordinates(F[s], plan8.pts, plan8.tabs, p) for s in SHAPES}
        co = tc.coordinates(o, plan8.pts, plan8.tabs, p)
        j5 = tc.Cochain(5, H); j23 = tc.Cochain(5, H)
        for succ, c in T5.items():
            j5.add_term(succ, {((), (1, 2, 3, 4)): 1}, c)
        for succ, c in T23.items():
            j23.add_term(succ, {((), (1, 2, 3, 4)): 1}, c)
        closed_j = [not np.any(plan6.apply(j, p)) for j in (j5, j23)]
        nonzero_j = [bool(np.any(tc.coordinates(j, plan6_src_pts, plan6_src_tabs, p))) for j in (j5, j23)]
        Q55 = tc.half_bracket_sym(T5, H, 120, N)
        Q2323 = tc.half_bracket_sym(T23, H, 120, N)
        Tmix = {k: T5.get(k, 0) + T23.get(k, 0) for k in set(T5) | set(T23)}
        Qmix = tc.half_bracket_sym(Tmix, H, 120, N)         # = Q55 + Q2323 + {j5, j23}
        cQ55 = tc.coordinates(Q55, plan8.pts, plan8.tabs, p)
        cQ2323 = tc.coordinates(Q2323, plan8.pts, plan8.tabs, p)
        cQ523 = (tc.coordinates(Qmix, plan8.pts, plan8.tabs, p) - cQ55 - cQ2323) % p
        closedQ = all(not np.any(plan9.apply(Q, p)) for Q in (Q55, Q2323, Qmix))
        ranks = {
            'B': rB,
            '[B;o]': np_rank_mod(np.vstack([B, co[None]]), p),
            '[B;F_G]': [np_rank_mod(np.vstack([B, cF[s][None]]), p) for s in SHAPES],
            '[B;F_all]': np_rank_mod(np.vstack([B] + [cF[s][None] for s in SHAPES]), p),
            '[B;Q55]': np_rank_mod(np.vstack([B, cQ55[None]]), p),
            '[B;Q2323]': np_rank_mod(np.vstack([B, cQ2323[None]]), p),
            '[B;Q523]': np_rank_mod(np.vstack([B, cQ523[None]]), p),
            '[B;Q*]': np_rank_mod(np.vstack([B, cQ55[None], cQ2323[None], cQ523[None]]), p),
            '[B;Q*;o]': np_rank_mod(np.vstack([B, cQ55[None], cQ2323[None], cQ523[None], co[None]]), p),
        }
        # solve o = x1 Q55 + x2 Q523 + x3 Q2323 + B y  (mod p)
        A = np.vstack([cQ55[None], cQ523[None], cQ2323[None], B]).T
        ok, x, _ = np_solve_mod(A, co, p)
        # F_G classes in terms of each other: F_G = sum_H a_H Q + B
        out[p] = dict(closed_j=closed_j, nonzero_j=nonzero_j, closed=closed, closedQ=closedQ, ranks=ranks, sol=(ok, None if x is None else [int(v) for v in x[:3]]))
    p0, p1 = PRIMES[:2]
    r = out[p0]
    sol = None
    if r['sol'][0] and out[p1]['sol'][0]:
        sol = []
        for a0, a1 in zip(r['sol'][1], out[p1]['sol'][1]):
            v, m = crt_pair(a0, p0, a1, p1)
            sol.append(ratrec(v, m))
    log('N =', N, ' closed F_G:', r['closed'], ' d j_5 = d j_23 = 0:', r['closed_j'], ' j nonzero:', r['nonzero_j'], ' closed 1/2{j,j}:', r['closedQ'])
    log('   ranks:', r['ranks'], ' (second prime equal:', r['ranks'] == out[p1]['ranks'], ')')
    log('   o\' = x55 Q55 + x523 {j5,j23} + x2323 Q2323 mod B:', r['sol'][0], ' x =', sol,
        ' (solution unique only if [B;Q*] has full rank +3)')
    results[N] = (r, sol)
log.close()
