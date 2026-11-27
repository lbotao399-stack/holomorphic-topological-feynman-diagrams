"""Formal sl_N (N = 9..17): class ratios [F_G] = r_G(N) [F_155] mod B = d(C^7_inv), [o'] = c(N) [F_155],
rational reconstruction in N; simple closed representatives (four-kappa and trace-product cochains) tested."""
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement
import sympy as sp
from common import Log, PRIMES, np_rank_mod, np_solve_mod, ratrec, crt_pair
import colour as col
import trace_complex as tc
import cartan_lib as cl

log = Log(__file__)
H = 4
rng = np.random.default_rng(77)
plan8 = tc.DiffPlan(7, H, 90, rng)
plan9 = tc.DiffPlan(8, H, 30, rng)
subs = list(combinations(range(6), 4))
gens = []
for lam in tc.fpf_types(7):
    for I, J in combinations_with_replacement(subs, 2):
        g = tc.Cochain(7, H)
        g.add_term(tc.canonical(lam), {(I, J): 1}, 1)
        if g.data:
            gens.append(g)
P8 = {((0, 1, 2, 3), (4, 5, 6, 7)): 1}


def simple_reps(N):
    """candidate closed cochains alt([0123][4567] T) with simple T (traceless sl_N)."""
    reps = {}
    # four kappa: tr(01)tr(23)tr(45)tr(67)
    reps['kkkk'] = [((1, 0, 3, 2, 5, 4, 7, 6), 1)]
    reps['tr(0123)tr(4567)'] = [((1, 2, 3, 0, 5, 6, 7, 4), 1)]
    reps['tr(01234567)'] = [((1, 2, 3, 4, 5, 6, 7, 0), 1)]
    reps['tr(04)tr(15)tr(26)tr(37)'] = [((4, 5, 6, 7, 0, 1, 2, 3), 1)]
    out = {}
    for name, terms in reps.items():
        om = tc.Cochain(8, H)
        for succ, c in terms:
            om.add_term(succ, P8, c)
        out[name] = om
    return out


rows = {}
for N in range(9, 18):
    res = {}
    for p in PRIMES[:2]:
        B = np.stack([plan8.apply(g, p) for g in gens])
        rB = np_rank_mod(B, p)
        F = {}
        for shape in cl.SHAPES:
            om = tc.Cochain(8, H)
            for perm, w in col.colour_perm_dict(shape, N).items():
                om.add_term(perm, {(I, J): e for I, J, e in cl.DTH[shape]}, cl.SIGMA * w)
            F[shape] = om
        cF = {s: tc.coordinates(F[s], plan8.pts, plan8.tabs, p) for s in cl.SHAPES}
        base = cF[(1, 5, 5)]
        ratios = {}
        A = np.vstack([base[None], B]).T
        for s in cl.SHAPES:
            ok, x, _ = np_solve_mod(A, cF[s], p)
            ratios[s] = int(x[0]) if ok else None
        reps = simple_reps(N)
        repinfo = {}
        for name, om in reps.items():
            cl_ = not np.any(plan9.apply(om, p))
            v = tc.coordinates(om, plan8.pts, plan8.tabs, p)
            nz = np_rank_mod(np.vstack([B, v[None]]), p) > rB
            ok, x, _ = np_solve_mod(np.vstack([v[None], B]).T, base, p) if cl_ and nz else (False, None, 0)
            repinfo[name] = (cl_, nz, int(x[0]) if ok else None)
        res[p] = (rB, ratios, repinfo)
    p0, p1 = PRIMES[:2]
    rat = {}
    for s in cl.SHAPES:
        a0, a1 = res[p0][1][s], res[p1][1][s]
        v, m = crt_pair(a0, p0, a1, p1)
        rat[s] = ratrec(v, m)
    c = sum(cl.KC[s] / cl.AUT[s] * rat[s] for s in cl.SHAPES)
    reps_q = {}
    for name in res[p0][2]:
        cl0, nz0, x0 = res[p0][2][name]
        x1 = res[p1][2][name][2]
        reps_q[name] = (cl0, nz0, ratrec(*crt_pair(x0, p0, x1, p1)) if x0 is not None and x1 is not None else None)
    rows[N] = (rat, c, reps_q)
    log('N =', N, ' rank B =', res[p0][0], ' r_G = [F_G]/[F_155]:', {''.join(map(str, s)): str(rat[s]) for s in cl.SHAPES},
        ' c(N) = [o\']/[F_155] =', c)
    log('     simple representatives (closed, non-exact, [F_155]/[rep]):', {k: (v[0], v[1], str(v[2])) for k, v in reps_q.items()})
# rational-function fit in N
Nsym = sp.Symbol('N')
for s in cl.SHAPES + ['c']:
    pts = [(N, rows[N][0][s] if s != 'c' else rows[N][1]) for N in sorted(rows)]
    for deg in range(0, 4):
        # fit numerator degree deg_n, denominator degree deg (monic-free): try polynomial times N^-k forms via sympy
        pass
    vals = [sp.Rational(v.numerator, v.denominator) for _, v in pts]
    f = sp.interpolate([(sp.Integer(N), v) for (N, _), v in zip(pts, vals)], Nsym)
    f = sp.factor(sp.simplify(f))
    # check rational form: multiply by N^k
    best = None
    for k in range(0, 5):
        g = sp.interpolate([(sp.Integer(N), v * sp.Integer(N) ** k) for (N, _), v in zip(pts, vals)], Nsym)
        if sp.degree(g, Nsym) <= len(pts) - 3:
            best = sp.factor(g / Nsym ** k)
            break
    log('fit', s if s == 'c' else ''.join(map(str, s)), ':', best if best is not None else f)
log.close()
