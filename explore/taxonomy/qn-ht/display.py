"""Explicit data: (1) z = 0 kernels K(lam;0)/pi of all rooted thetas at H = 1,2,3 in the hteng labelling;
(2) Q_0 on undifferentiated letters; (3) Q_2^hat on undifferentiated monomials of degree 2, 3 (fan-type and
tadpole-type parts) for the theories of THEORIES.  Log: display.log"""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import sympy as sp
import hteng, thetak, theories as T, specs

log = open("display.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

out("## (1) z = 0 kernels K(lam;0)/pi, labelled graphs of hteng (root = nv-1, edges in contraction order)")
for H in (1, 2, 3):
    for sh in hteng.theta_shapes(H):
        Tt = sh.theta
        mom = Tt.moment(tuple((0,) * H for _ in Tt.root_edges))
        lam = {u: [sp.Symbol("l%d%s" % (u, "xyz"[i] if H > 1 else "")) for i in range(H)] for u in Tt.surv}
        expr = 0
        for beta, c in mom.items():
            t = sp.Rational(c.numerator, c.denominator)
            for k, u in enumerate(Tt.surv):
                for i in range(H):
                    t *= lam[u][i] ** beta[k][i]
            expr += t
        out("H=%d %-18s |Aut_O|=%d edges %s survivors %s: K/pi = %s" % (H, sh.name, sh.aut, sh.edges, sh.surv, sp.factor(expr)))

THEORIES = [
    ("3d N=4 SQED Nf=1", lambda: T.sqed_3d(1)),
    ("3d N=8 SU(2) (so3)", lambda: T.su2_N8_3d()),
    ("5d N=1 U(1) + hyper q (BF R/g), q = 1", lambda: hteng.Theory(2, 1, {}, bf=T.u1_reps([1]))),
    ("5d N=1 SU(2) pure (so3)", lambda: hteng.Theory(2, 3, T.F_SO3)),
    ("5d N=2 SU(2) (so3, BF g/g)", lambda: hteng.Theory(2, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))])),
    ("7d N=1 SU(2) (so3, BF Bg)", lambda: hteng.Theory(3, 3, T.F_SO3)),
]
for name, mk in THEORIES:
    th = mk()
    E = hteng.Engine(th)
    out("\n## %s, (H,T) = (%d,1); interaction I = %s" % (name, th.H, hteng.pstr(th.I, 60)))
    pool = th.letter_pool()
    for L in pool:
        q = E.Q0(th.mono(L))
        out("  Q0 %s = %s" % (hteng.lname(L), hteng.pstr(q, 30)))
    fan = [k for k, s in enumerate(E.shapes) if s.name.endswith("@N")]
    tad = [k for k, s in enumerate(E.shapes) if not s.name.endswith("@N")]
    shown = 0
    for deg in (2, 3):
        for combo in itertools.combinations_with_replacement(pool, deg):
            O = th.mono(*combo)
            if not O:
                continue
            qf, qt = E.Q2(O, fan), E.Q2(O, tad)
            if not qf and not qt:
                continue
            if shown >= 14:
                break
            shown += 1
            out("  Q2^hat(%s):" % hteng.pstr(O))
            if qf:
                out("     fan-type     = %s" % hteng.pstr(qf, 16))
            if qt:
                out("     tadpole-type = %s" % hteng.pstr(qt, 16))
