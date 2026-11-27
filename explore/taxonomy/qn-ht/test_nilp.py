"""{Q_0, Q_2} = 0 on all monomials of bounded degree and total derivative order (or a random sample),
for the theories of THEORIES.  Also: fan-type (O trivalent) and tadpole-type (O bivalent) parts
separately; nontriviality (Q0 Q2 and Q2 Q0 both != 0); negative control (one shape weight doubled).
Usage: python3 test_nilp.py <theory> <maxdeg> <maxd> [nsample seed]   Log: test_nilp_<theory>_<maxdeg>_<maxd>.log"""
import sys, time, itertools, random
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng, theories as T

THEORIES = {
    # 3d (H=1)
    "3d_sqed1": lambda: T.sqed_3d(1),
    "3d_sqed2": lambda: T.sqed_3d(2),
    "3d_su2N4_1": lambda: T.su2_N4_3d(1),
    "3d_su2N8": lambda: T.su2_N8_3d(),
    # 5d (H=2)
    "5d_pure_so3": lambda: hteng.Theory(2, 3, T.F_SO3, name="5d N=1 pure SU(2)"),
    "5d_N2_so3": lambda: hteng.Theory(2, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))], name="5d N=2 SU(2) (BF g/g)"),
    "5d_sl2_fund": lambda: hteng.Theory(2, 3, T.F_SL2, bf=[(2, T.FUND_SL2)], name="5d N=1 SU(2) + 1 hyper (BF R/g, R = fund)"),
    "5d_sl2_2fund": lambda: hteng.Theory(2, 3, T.F_SL2, bf=[(2, T.FUND_SL2), (2, T.FUND_SL2)], name="5d N=1 SU(2) + 2 hypers (BF R/g)"),
    "5d_u1_q1": lambda: hteng.Theory(2, 1, {}, bf=T.u1_reps([1]), name="5d N=1 U(1) + 1 hyper q=1 (BF R/g)"),
    "5d_u1_q12": lambda: hteng.Theory(2, 1, {}, bf=T.u1_reps([1, 2]), name="5d N=1 U(1) + hypers q=1,2 (BF R/g)"),
    "5d_half_sl2": lambda: T.sym_theory(2, 'sl2', [T.half_hyper_sl2()], "5d N=1 SU(2) + half-hyper fund (U//g)"),
    "5d_fullsym_sl2": lambda: T.sym_theory(2, 'sl2', [T.full_hyper_sym(T.FUND_SL2, 2)], "5d N=1 SU(2) + 1 hyper (U//g, U = T^*fund)"),
    # 7d (H=3)
    "7d_so3": lambda: hteng.Theory(3, 3, T.F_SO3, name="7d N=1 SU(2) (BF Bg)"),
    "7d_sl2": lambda: hteng.Theory(3, 3, T.F_SL2, name="7d N=1 SU(2) (BF Bg, sl2 basis)"),
}


def run(name, maxdeg, maxd, nsample=0, seed=1, logf=None):
    th = THEORIES[name]()
    H = th.H
    def out(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        if logf:
            logf.write(s + "\n"); logf.flush()
    t0 = time.time()
    E = hteng.Engine(th)
    out("theory:", th.name, "(H,T) = (%d,1), D = %d; interaction terms %d; shapes:" % (H, H + 1, len(th.I)),
        ", ".join("%s |Aut|=%d" % (s.name, s.aut) for s in E.shapes), "[setup %.1fs]" % (time.time() - t0))
    fan = [k for k, s in enumerate(E.shapes) if s.name.endswith("@N")]
    tad = [k for k, s in enumerate(E.shapes) if not s.name.endswith("@N")]
    pool = []
    for (n, i, _) in th.letter_pool():
        for a in itertools.product(range(maxd + 1), repeat=H):
            if sum(a) <= maxd:
                pool.append((n, i, a))
    monos = []
    for deg in range(2, maxdeg + 1):
        for combo in itertools.combinations_with_replacement(pool, deg):
            if sum(sum(L[2]) for L in combo) > maxd:
                continue
            O = th.mono(*combo)
            if O:
                monos.append(O)
    if nsample and len(monos) > nsample:
        monos = random.Random(seed).sample(monos, nsample)
    nq2 = nnt = nbad = nbadF = nbadT = ntad = 0
    ctrl_bad = {}
    for O in monos:
        q2 = E.Q2(O)
        a = E.Q0(q2)
        b = E.Q2(E.Q0(O))
        C = hteng.psum(a, b)
        nq2 += bool(q2)
        nnt += bool(a) and bool(b)
        if C:
            nbad += 1
            if nbad <= 3:
                out("  FAIL on", hteng.pstr(O), ":", hteng.pstr(C, 8))
        qf = E.Q2(O, fan); qt = E.Q2(O, tad)
        ntad += bool(qt)
        if hteng.psum(E.Q0(qf), E.Q2(E.Q0(O), fan)):
            nbadF += 1
        if hteng.psum(E.Q0(qt), E.Q2(E.Q0(O), tad)):
            nbadT += 1
        # negative control: shape k weight doubled
        for k in range(len(E.shapes)):
            qk = E.Q2(O, [k])
            Ck = hteng.psum(C, E.Q0(qk), E.Q2(E.Q0(O), [k]))
            if Ck:
                ctrl_bad[k] = ctrl_bad.get(k, 0) + 1
    out("%d monomials (deg 2..%d, total derivatives <= %d%s): Q2^hat != 0 on %d (tadpole-type part != 0 on %d); "
        "Q0 Q2 and Q2 Q0 both != 0 on %d; {Q0,Q2} != 0 on %d" % (len(monos), maxdeg, maxd,
        ", random sample seed %d" % seed if nsample else "", nq2, ntad, nnt, nbad))
    out("   fan-type part alone: {Q0,.} != 0 on %d; tadpole-type part alone: %d" % (nbadF, nbadT))
    out("   control (weight of one shape doubled): failures", {E.shapes[k].name: v for k, v in sorted(ctrl_bad.items())})
    out("   [%.0fs]" % (time.time() - t0))
    return nbad


if __name__ == "__main__":
    name, maxdeg, maxd = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nsample = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    f = open("test_nilp_%s_%d_%d%s.log" % (name, maxdeg, maxd, "_s%d" % seed if nsample else ""), "w")
    run(name, maxdeg, maxd, nsample, seed, f)
