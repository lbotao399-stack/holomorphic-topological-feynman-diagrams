"""Formal sl_N, 9d weight (2,2,2,2): H^8_inv = 4 (h8_dims).  Trace-product representatives
R_kkkk = alt([0123][4567] tr01 tr23 tr45 tr67), R_44 = alt([0123][4567] tr(0123) tr(4567)),
R_x = alt([0123][4567] tr04 tr15 tr26 tr37); decomposition of [F_155] and of the bracket classes."""
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement
from common import Log, PRIMES, np_rank_mod, np_solve_mod, ratrec, crt_pair
import colour as col
import trace_complex as tc
import cartan_lib as cl

log = Log(__file__)
H = 4
rng = np.random.default_rng(5)
plan8 = tc.DiffPlan(7, H, 90, rng)
subs = list(combinations(range(6), 4))
gens = []
for lam in tc.fpf_types(7):
    for I, J in combinations_with_replacement(subs, 2):
        g = tc.Cochain(7, H); g.add_term(tc.canonical(lam), {(I, J): 1}, 1)
        if g.data:
            gens.append(g)
P8 = {((0, 1, 2, 3), (4, 5, 6, 7)): 1}
REPS = {'kkkk': (1, 0, 3, 2, 5, 4, 7, 6), 't44': (1, 2, 3, 0, 5, 6, 7, 4), 'x': (4, 5, 6, 7, 0, 1, 2, 3)}
T5, T23 = tc.sym_trace_dict((5,)), tc.sym_trace_dict((2, 3))
Tmix = {k: T5.get(k, 0) + T23.get(k, 0) for k in set(T5) | set(T23)}
for N in (9, 12, 15):
    sols = {}
    for p in PRIMES[:2]:
        B = np.stack([plan8.apply(g, p) for g in gens])
        rB = np_rank_mod(B, p)
        F = tc.Cochain(8, H)
        for perm, w in col.colour_perm_dict((1, 5, 5), N).items():
            F.add_term(perm, P8, cl.SIGMA * w)
        vF = tc.coordinates(F, plan8.pts, plan8.tabs, p)
        vR = {}
        for k, succ in REPS.items():
            om = tc.Cochain(8, H); om.add_term(succ, P8, 1)
            vR[k] = tc.coordinates(om, plan8.pts, plan8.tabs, p)
        q55 = tc.coordinates(tc.half_bracket_sym(T5, H, 120, N), plan8.pts, plan8.tabs, p)
        q2323 = tc.coordinates(tc.half_bracket_sym(T23, H, 120, N), plan8.pts, plan8.tabs, p)
        qm = (tc.coordinates(tc.half_bracket_sym(Tmix, H, 120, N), plan8.pts, plan8.tabs, p) - q55 - q2323) % p
        vR['Q55'], vR['Q523'], vR['Q2323'] = q55, qm, q2323
        names = list(vR)
        r_all = np_rank_mod(np.vstack([B] + [vR[k][None] for k in names]), p) - rB
        vR['F155'] = vF
        allnames = list(vR)
        # greedy basis preferring trace-product representatives
        order = ['kkkk', 't44', 'x', 'Q2323', 'Q523', 'Q55', 'F155']
        basis, rcur = [], rB
        for k in order:
            rr = np_rank_mod(np.vstack([B] + [vR[j][None] for j in basis + [k]]), p)
            if rr > rcur:
                basis.append(k); rcur = rr
        sols.setdefault('basis', basis)
        out = {}
        for target in [t for t in ['F155', 'Q55', 'Q523', 'Q2323', 'kkkk', 't44', 'x'] if t not in basis]:
            v = vF if target == 'F155' else vR[target]
            A = np.vstack([vR[k][None] for k in basis] + [B]).T
            ok, x, _ = np_solve_mod(A, v, p)
            out[target] = [int(a) for a in x[:len(basis)]] if ok else None
        sols[p] = (r_all, None, out)
    p0, p1 = PRIMES[:2]
    r_all, rb, out0 = sols[p0]
    rb = None
    out1 = sols[p1][2]
    rec = {}
    for t in out0:
        if out0[t] is None or out1[t] is None:
            rec[t] = None
        else:
            rec[t] = [str(ratrec(*crt_pair(a, p0, b, p1))) for a, b in zip(out0[t], out1[t])]
    log('N = %d: rank of {kkkk,t44,x,Q55,Q523,Q2323} mod B = %d' % (N, r_all))
    log('   H^8 basis (greedy, trace products first):', sols['basis'], ' coefficients mod B:', rec)
log.close()
