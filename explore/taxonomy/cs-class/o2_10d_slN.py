"""10d holomorphic twist = hCS on C^5, (H,T) = (5,0), pure open two-loop theta sum for formal sl_N (kappa = tr):
  F_G = alt(Q_G(Lbar) W^G),  o'_10 = sum_G |Aut G|^{-1} alt((K_G) W^G),  K_G = s_G c_0,G Q_G(Lbar)
Q(Lbar) expanded into products of two 5-minors of the eight interior momenta (Lbar = sum mu_v p_v).
Checks: d F_G and d o'_10 at C^9 coordinates (closedness); the type-I README states d != 0 for so(N)."""
import sys
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations
from common import Log, PRIMES, np_rank_mod, parity
import colour as col
import trace_complex as tc
import twoloop

log = Log(__file__)
H = 5
SHAPES = [(1, 5, 5), (2, 4, 5), (3, 3, 5), (3, 4, 4)]
KC0 = {(1, 5, 5): Fr(1, 17280), (2, 4, 5): Fr(-1, 17280), (3, 3, 5): Fr(1, 17280), (3, 4, 4): Fr(-1, 25920)}
AUT = {(1, 5, 5): 4, (2, 4, 5): 2, (3, 3, 5): 4, (3, 4, 4): 4}


def q_poly(shape):
    """Q(Lbar) as {(I,J): coef}, letters 0..7 = interior vertices 1..8."""
    th = twoloop.Theta(*twoloop.theta_graph(*shape))
    mu = twoloop.mu_weights(th)
    A, B, C = th.interior['A'], th.interior['B'], th.interior['C']
    # columns: ('vt'), ('vb'), A (top), B (both), C (bottom); 10 rows: top 0..4, bottom 5..9
    cols = [('vt', None), ('vb', None)] + [('A', u) for u in A] + [('B', u) for u in B] + [('C', u) for u in C]
    out = {}
    for S in combinations(range(10), 5):
        types = [cols[c][0] for c in S]
        if any(t in ('vb', 'C') for t in types):
            continue
        Sc = [c for c in range(10) if c not in S]
        if any(cols[c][0] in ('vt', 'A') for c in Sc):
            continue
        sign = (-1) ** (sum(range(5)) + sum(S))
        # top minor: columns S (vt -> Lbar), bottom: Sc (vb -> Lbar)
        for k, mk in mu.items():
            if mk == 0 or k in (th.N, th.Np):
                continue
            for l, ml in mu.items():
                if ml == 0 or l in (th.N, th.Np):
                    continue
                top = [k if cols[c][0] == 'vt' else cols[c][1] for c in S]
                bot = [l if cols[c][0] == 'vb' else cols[c][1] for c in Sc]
                if len(set(top)) < 5 or len(set(bot)) < 5:
                    continue
                top0 = [x - 1 for x in top]; bot0 = [x - 1 for x in bot]
                s1 = parity(top0); s2 = parity(bot0)
                I, J = tuple(sorted(top0)), tuple(sorted(bot0))
                key = (I, J) if I <= J else (J, I)
                out[key] = out.get(key, 0) + sign * s1 * s2 * mk * ml
    return {k: v for k, v in out.items() if v}, th


# verify the minor expansion against twoloop.Q_form(Lbar) at random sources
rng0 = np.random.default_rng(3)
QP = {}
for shape in SHAPES:
    P, th = q_poly(shape)
    QP[shape] = P
    for _ in range(3):
        lam = [[Fr(int(x)) for x in rng0.integers(-4, 5, 5)] for _ in range(th.nv)]
        lam[th.N] = [Fr(0)] * 5
        lam[th.Np] = [-sum(lam[v][i] for v in range(th.nv) if v != th.Np) for i in range(5)]
        ref = twoloop.Q_form(th, lam, twoloop.Lbar(th, lam))
        from common import det_int
        val = sum(c * det_int([lam[i + 1] for i in I]) * det_int([lam[j + 1] for j in J]) for (I, J), c in P.items())
        assert val == ref, (shape, val, ref)
    log(shape, 'Q(Lbar) minor expansion: %d terms, equal to twoloop.Q_form at 3 random source sets' % len(P))

NS = [int(x) for x in sys.argv[1:]] or [9, 11]
rng = np.random.default_rng(10)
plan9 = tc.DiffPlan(8, H, 40, rng)
for N in NS:
    for p in PRIMES[:1]:
        F = {}
        for shape in SHAPES:
            sig = col.standard_sigma(shape)
            om = tc.Cochain(8, H)
            for perm, w in col.colour_perm_dict(shape, N).items():
                om.add_term(perm, QP[shape], sig * w)
            F[shape] = om
        o = F[SHAPES[0]].scaled(KC0[SHAPES[0]] / AUT[SHAPES[0]])
        for shape in SHAPES[1:]:
            o = o.plus(F[shape], KC0[shape] / AUT[shape])
        dF = {s: int(np.count_nonzero(plan9.apply(F[s], p))) for s in SHAPES}
        do = plan9.apply(o, p)
        log('N =', N, ' nonzero C^9 coordinates of d F_G:', dF, ' of d o\'_10:', int(np.count_nonzero(do)), 'of', len(do))
log.close()
