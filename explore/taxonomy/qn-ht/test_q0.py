"""Q_0^2 = 0 for each theory and each choice of the interaction signs (sB, sM, sU) (and eta signs);
records which choices are consistent.  Usage: python3 test_q0.py   Log: test_q0.log"""
import sys, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng, theories as T

log = open("test_q0.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

def q0sq_fail(th, maxdeg, maxd):
    E = hteng.Engine(th, shapes=[])
    H = th.H
    pool = []
    for (n, i, _) in th.letter_pool():
        for a in itertools.product(range(maxd + 1), repeat=H):
            if sum(a) <= maxd:
                pool.append((n, i, a))
    nm = nf = nz = 0
    for deg in range(1, maxdeg + 1):
        for combo in itertools.combinations_with_replacement(pool, deg):
            if sum(sum(L[2]) for L in combo) > maxd:
                continue
            O = th.mono(*combo)
            if not O:
                continue
            nm += 1
            q = E.Q0(O)
            nz += bool(q)
            if E.Q0(q):
                nf += 1
    return nm, nz, nf

def mk(kind, H, sg):
    sB, sM, sU = sg
    if kind == "BF adj":
        return hteng.Theory(H, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))], signs=(sB, sM, 1, sU))
    if kind == "BF sl2 fund":
        return hteng.Theory(H, 3, T.F_SL2, bf=[(2, T.FUND_SL2)], signs=(sB, sM, 1, sU))
    if kind == "BF u1 q=1,2":
        return hteng.Theory(H, 1, {}, bf=T.u1_reps([1, 2]), signs=(sB, sM, 1, sU))
    if kind == "pure so3":
        return hteng.Theory(H, 3, T.F_SO3, signs=(sB, sM, 1, sU))
    if kind == "U//g sl2 half-hyper":
        return T.sym_theory(H, 'sl2', [T.half_hyper_sl2()], "half", signs=(sB, sM, 1, sU))
    if kind == "U//g sl2 full hyper":
        return T.sym_theory(H, 'sl2', [T.full_hyper_sym(T.FUND_SL2, 2)], "full", signs=(sB, sM, 1, sU))

for H, kinds in [(1, ["BF adj", "BF sl2 fund", "BF u1 q=1,2"]),
                 (2, ["pure so3", "BF adj", "BF sl2 fund", "BF u1 q=1,2", "U//g sl2 half-hyper", "U//g sl2 full hyper"]),
                 (3, ["pure so3", "BF adj", "BF sl2 fund"])]:
    for kind in kinds:
        res = []
        for sg in itertools.product([1, -1], repeat=3):
            if "U//g" not in kind and sg[2] == -1:
                continue
            if "U//g" in kind and sg[1] == -1:
                continue
            if kind == "pure so3" and sg != (1, 1, 1):
                continue
            th = mk(kind, H, sg)
            t0 = time.time()
            nm, nz, nf = q0sq_fail(th, 3, 1 if H < 3 else 1)
            res.append((sg, nm, nz, nf))
            out("H=%d D=%d %-22s signs (sB,sM,sU)=%s: %d monomials (deg<=3, derivatives<=1), Q0 != 0 on %d, Q0^2 != 0 on %d  [%.1fs]"
                % (H, H + 1, kind, sg, nm, nz, nf, time.time() - t0))
