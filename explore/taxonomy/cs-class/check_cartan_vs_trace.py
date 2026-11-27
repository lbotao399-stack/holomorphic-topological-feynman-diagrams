"""Cross-check of the Cartan detector (root-space contraction, set-partition reduction of the alternation) against
the trace-basis cochain o' of o2_9d_slN (ribbon colour, add_term alternation) evaluated on diagonal sl_N inputs:
o'(h_v e^{p_v z}) = sum_rho c_rho(p) T_rho(h), all fixed-point-free rho in S_8, exact mod p."""
import itertools
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, det_int
import colour as col
import trace_complex as tc
import cartan_lib as cl
import lie

log = Log(__file__)
p = PRIMES[0]
H = 4
for N in (3, 4):
    L = lie.build('A%d' % (N - 1))
    detc = cl.Detector(L, p)
    rng = np.random.default_rng(N)
    o = None
    for shape in cl.SHAPES:
        om = tc.Cochain(8, H)
        for perm, w in col.colour_perm_dict(shape, N).items():
            om.add_term(perm, {(I, J): e for I, J, e in cl.DTH[shape]}, cl.SIGMA * w)
        om = om.scaled(cl.KC[shape] / cl.AUT[shape])
        o = om if o is None else o.plus(om)
    rhos = [r for r in itertools.permutations(range(8)) if all(r[i] != i for i in range(8))]
    for t in range(2):
        hc, mom = detc.sample(rng)
        v_det = detc.o_prime(hc, mom)
        diag = []
        for v in range(8):
            d = [0] * N
            for i in range(N - 1):
                d[i] += int(hc[v, i]); d[i + 1] -= int(hc[v, i])
            diag.append(d)
        minors = {(): 1}
        for I in itertools.combinations(range(8), 4):
            minors[I] = det_int([[int(mom[i][k]) for k in range(4)] for i in I])
        total = 0
        cache = {}
        for rho in rhos:
            lam, sig0 = tc.conjugator(rho)
            P = o.data.get(lam)
            if not P:
                continue
            Q = tc.poly_pullback(P, sig0, H)
            cv = sum(cl.fmod(c, p) * (minors[I] * minors[J] % p) for (I, J), c in Q.items()) * tc.sgn(sig0) % p
            Tv = 1
            for cyc in tc.cycles(rho):
                Tv *= sum(np.prod([diag[a][i] for a in cyc]) for i in range(N))
            total = (total + cv * (Tv % p)) % p
        log('sl_%d sample %d: Cartan detector = %d, trace-basis cochain = %d, equal: %s' % (N, t, v_det, total, v_det == total))
        assert v_det == total
log('## 1/2{j_T,j_T} on Cartan cycles: cartan_lib.half_bracket_tensor vs trace-basis half_bracket_sym')
for N, Tname in ((4, (2, 3)), (5, (5,)), (5, (2, 3))):
    L = lie.build('A%d' % (N - 1))
    W = [[(1 if k == i else 0) - (1 if k == i + 1 else 0) for i in range(N - 1)] for k in range(N)]
    if Tname == (5,):
        Pf = lambda x, W=W: sum(sum(w[i] * x[i] for i in range(len(x))) ** 5 for w in W)
    else:
        Pf = lambda x, W=W: sum(sum(w[i] * x[i] for i in range(len(x))) ** 2 for w in W) * sum(sum(w[i] * x[i] for i in range(len(x))) ** 3 for w in W)
    Tt = cl.sym_tensor_from_poly(Pf, N - 1, p)
    Q = tc.half_bracket_sym(tc.sym_trace_dict(Tname), H, 120, N)
    detc = cl.Detector(L, p)
    rng = np.random.default_rng(50 + N)
    rhos = [r for r in itertools.permutations(range(8)) if all(r[i] != i for i in range(8))]
    for t in range(2):
        hc, mom = detc.sample(rng)
        v_c = cl.half_bracket_tensor(L, Tt, hc, mom, p)
        diag = []
        for v in range(8):
            d = [0] * N
            for i in range(N - 1):
                d[i] += int(hc[v, i]); d[i + 1] -= int(hc[v, i])
            diag.append(d)
        minors = {(): 1}
        for I in itertools.combinations(range(8), 4):
            minors[I] = det_int([[int(mom[i][k]) for k in range(4)] for i in I])
        total = 0
        for rho in rhos:
            lam, sig0 = tc.conjugator(rho)
            P = Q.data.get(lam)
            if not P:
                continue
            Qp = tc.poly_pullback(P, sig0, H)
            cv = sum(cl.fmod(c, p) * (minors[I] * minors[J] % p) for (I, J), c in Qp.items()) * tc.sgn(sig0) % p
            Tv = 1
            for cyc in tc.cycles(rho):
                Tv *= sum(np.prod([diag[a][i] for a in cyc]) for i in range(N))
            total = (total + cv * (Tv % p)) % p
        log('sl_%d T=%s sample %d: Cartan = %d, trace basis = %d, equal: %s' % (N, Tname, t, v_c, total, v_c == total))
        assert v_c == total
log.close()
