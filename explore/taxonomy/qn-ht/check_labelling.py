"""Labelling independence of Q_2^hat at H = 2, 3: every rooted theta relabelled by a random connected
search order (new vertex labels, edge order (-v,w), orientations small -> large, kernel sign refitted
against htk); Q2^hat compared on all monomials of degree <= 3, derivatives <= maxd.  At D odd the rule's
vertex Koszul signs (I odd) and the kernel's sgn(pi)^D must compensate.  Log: check_labelling.log"""
import sys, itertools, time
sys.path.insert(0, ".")
import hteng, theories as T

log = open("check_labelling.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

cases = [("5d SU(2) + 1 fund hyper (BF R/g)", lambda: hteng.Theory(2, 3, T.F_SL2, bf=[(2, T.FUND_SL2)]), 1),
         ("5d SU(2) + half-hyper (U//g)", lambda: T.sym_theory(2, 'sl2', [T.half_hyper_sl2()], "half"), 1),
         ("7d SU(2) (BF Bg)", lambda: hteng.Theory(3, 3, T.F_SO3), 0),
         ("3d SU(2) N=4 Nf=1", lambda: T.su2_N4_3d(1), 1)]
for name, mk, maxd in cases:
    th = mk()
    E0 = hteng.Engine(th)
    for oseed in (101, 202):
        t0 = time.time()
        E1 = hteng.Engine(th, shapes=hteng.theta_shapes(th.H, order_seed=oseed))
        diff_lab = sum(1 for a, b in zip(E0.shapes, E1.shapes) if a.edges != b.edges)
        H = th.H
        pool = [(n, i, a) for (n, i, _) in th.letter_pool() for a in itertools.product(range(maxd + 1), repeat=H) if sum(a) <= maxd]
        n = bad = nz = 0
        for deg in (2, 3):
            for combo in itertools.combinations_with_replacement(pool, deg):
                if sum(sum(L[2]) for L in combo) > maxd:
                    continue
                O = th.mono(*combo)
                if not O:
                    continue
                n += 1
                q0, q1 = E0.Q2(O), E1.Q2(O)
                nz += bool(q0)
                if hteng.psum(q0, {m: -c for m, c in q1.items()}):
                    bad += 1
        out("%s, relabelling seed %d (%d of %d rooted shapes with different labelled edge lists): %d monomials, "
            "Q2^hat != 0 on %d, mismatches %d  [%.0fs]" % (name, oseed, diff_lab, len(E0.shapes), n, nz, bad, time.time() - t0))
