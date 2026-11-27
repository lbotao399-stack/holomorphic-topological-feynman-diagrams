"""Tadpole-type part of Q_2^hat (O bivalent; one-loop (D+1)-cycle not through O) at (H,1):
additivity over multiplets and proportionality to the pure-gauge part.
5d (H=2): SU(2) (sl2) pure, + Nf fundamental hypers (BF R/g, Nf = 1, 2), + adjoint (5d N=2), + half-hyper
(U//g); 3d (H=1), 7d (H=3): pure and adjoint.  Monomials: all quadratic monomials with total derivative
order <= maxd.  Fit T(theory) = x * T(pure) on every coefficient.  Log: tadpole_trace.log"""
import sys, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng, theories as T

log = open("tadpole_trace.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

def tad_parts(th, maxd, shapes=None):
    E = hteng.Engine(th, shapes=shapes)
    tad = [k for k, s in enumerate(E.shapes) if not s.name.endswith("@N")]
    H = th.H
    gl = [L for L in th.letter_pool() if L[0] in ('c', 'B')]
    pool = [(n, i, a) for (n, i, _) in gl for a in itertools.product(range(maxd + 1), repeat=H) if sum(a) <= maxd]
    res = {}
    for combo in itertools.combinations_with_replacement(pool, 2):
        if sum(sum(L[2]) for L in combo) > maxd:
            continue
        O = th.mono(*combo)
        if not O:
            continue
        (m, c), = O.items()
        res[m] = E.Q2(O, tad)
    return res, E.shapes

def ratio(Tth, Tpure):
    xs = set()
    nz = 0
    for m in set(Tth) | set(Tpure):
        a, b = Tth.get(m, {}), Tpure.get(m, {})
        for mm in set(a) | set(b):
            u, v = a.get(mm, 0), b.get(mm, 0)
            if v == 0 and u == 0:
                continue
            nz += 1
            xs.add(None if v == 0 else Fr(u) / v)
    return xs, nz

for H, maxd in [(1, 1), (2, 1), (3, 0)]:
    t0 = time.time()
    if H == 2:
        cases = [
            ("pure", hteng.Theory(2, 3, T.F_SL2)),
            ("+1 fund hyper (BF R/g)", hteng.Theory(2, 3, T.F_SL2, bf=[(2, T.FUND_SL2)])),
            ("+2 fund hypers (BF R/g)", hteng.Theory(2, 3, T.F_SL2, bf=[(2, T.FUND_SL2)] * 2)),
            ("+adjoint (5d N=2, BF g/g)", hteng.Theory(2, 3, T.F_SL2, bf=[(3, T.adj_rho(T.F_SL2, 3))])),
            ("+half-hyper fund (U//g)", T.sym_theory(2, 'sl2', [T.half_hyper_sl2()], "half")),
            ("+1 fund hyper (U//g, U = T^*fund)", T.sym_theory(2, 'sl2', [T.full_hyper_sym(T.FUND_SL2, 2)], "full")),
        ]
    else:
        cases = [("pure", hteng.Theory(H, 3, T.F_SL2)),
                 ("+adjoint", hteng.Theory(H, 3, T.F_SL2, bf=[(3, T.adj_rho(T.F_SL2, 3))])),
                 ("+fund", hteng.Theory(H, 3, T.F_SL2, bf=[(2, T.FUND_SL2)]))]
    shapes = hteng.theta_shapes(H)
    base = None
    for name, th in cases:
        Tt, _ = tad_parts(th, maxd, shapes)
        nnz = sum(1 for v in Tt.values() if v)
        if base is None:
            base = Tt
            out("H=%d SU(2) %s: tadpole-type part nonzero on %d of %d quadratic (c,B) monomials (derivatives <= %d)"
                % (H, name, nnz, len(Tt), maxd))
            for m, v in sorted(Tt.items())[:200]:
                if v and sum(sum(L[2]) for L in m) == 0:
                    out("    Q2^hat_tad(%s) = %s" % (" ".join(hteng.lname(L) for L in m), hteng.pstr(v, 12)))
            continue
        xs, nz = ratio(Tt, base)
        out("H=%d SU(2) %s: tadpole-type part nonzero on %d monomials; T / T_pure on %d coefficients: %s"
            % (H, name, nnz, nz, sorted(xs, key=lambda x: (x is None, x))))
    out("   [%.0fs]" % (time.time() - t0))
