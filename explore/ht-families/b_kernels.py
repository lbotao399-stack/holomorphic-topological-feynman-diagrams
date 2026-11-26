"""(b) K(lambda; 0) at n = 4, (H,T) = (1,1): the five series-parallel 2-Laman graphs with 6 vertices.

 B1  engine kernel with symbolic sources at all non-root vertices (plan 'large'), cached in b_symbolic.txt
 B2  my exact pointwise check: engine ratio density g(xi) == D_brute(t(xi)) |det[t; dt/dxi]| (sign and square),
     D_brute from hf_brute (own Kirchhoff solution, BFS tree), t(xi) from my own ratio map
 B3  engine kernel at integer sources == symbolic kernel evaluated there
 B4  QMC integral of the float brute-force density over the ratio coordinates (own map), versus B3
"""
import sys, random, time, json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import hf_brute as hb
import hf_numeric as hn

SYM = {
    0: "pi**2*l0*l1*l2*l3/5",
    1: "-pi**2*l1*l4*(3*l0 - 3*l2 + 7*l3 + 5*l4)*(l0 + l1 + l2 + l3 + l4)/90",
    2: "0",
    3: "-pi**2*l2*l3*(8*l0**2 + 20*l0*l1 + 18*l0*l2 + 20*l0*l3 + 20*l0*l4 + 12*l1**2 + 22*l1*l2 + 22*l1*l3 + 20*l1*l4 + 10*l2**2 + 21*l2*l3 + 20*l2*l4 + 10*l3**2 + 18*l3*l4 + 8*l4**2)/180",
    4: "-pi**2*l1*(l0 + l1 + l2 + l3 + l4)*(l0*l1 + 2*l0*l3 - 9*l1**2 - 17*l1*l2 - 19*l1*l3 - l1*l4 - 8*l2**2 - 18*l2*l3 - 10*l3**2 - 2*l3*l4)/180",
}


def plan_of(red):
    return [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]


def main(which=range(5), qmc_m=15, nrep=16):
    gs = [Gx for Gx in eg.enumerate_tight(2, 4) if eg.is_sp(Gx)]
    rng = random.Random(7)
    lsyms = sp.symbols("l0:5")
    for k in which:
        G, _ = htk.from_networkx(gs[k], None)
        degs = sorted(G.degree(v) for v in range(G.nv))
        print("=" * 100)
        print("graph %d: degrees %s, edges %s, root %d" % (k, degs, G.edges, G.root))
        red = spred.sp_reduce(G, prefer="large")
        plan = plan_of(red)
        print("  plan (p, q, chains):", [(p, q, ch) for p, q, ch, s in plan])
        Ksym = sp.sympify(SYM[k], locals={str(s): s for s in lsyms})
        print("  B1  K_sym =", Ksym)
        Gm = hb.G(G.nv, G.edges)
        C = Gm.cycles()
        rm = hn.RatioMap(G.ne, plan)
        assert rm.names == red.bary and rm.coords == red.coords
        for trial in range(2):
            lam = [[Fr(rng.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]))] for _ in range(G.nv - 1)]
            t0 = time.time()
            dens = red.density(lam, T=1, check=0)
            Knum = dens.integrate()
            # B2: pointwise, own brute force at own t(xi)
            ok = True
            for _ in range(4):
                pt = red.random_point(rng)
                t_own = [Fr(1)] * G.ne
                for e in range(G.ne):
                    for b, ex in rm.texp[e].items():
                        t_own[e] *= Fr(pt[b]) ** ex
                assert t_own == red.tmap(pt)
                # exact Jacobian det[t; dt/dxi]
                rows = [t_own]
                for c in rm.coords:
                    blk = next(bl for bl in rm.blocks if c in bl)
                    row = []
                    for e in range(G.ne):
                        s = Fr(rm.texp[e].get(c, 0)) / Fr(pt[c]) - Fr(rm.texp[e].get(blk[0], 0)) / Fr(pt[blk[0]])
                        row.append(t_own[e] * s)
                    rows.append(row)
                J = hb.fdet(rows)
                coeff, U = hb.density(Gm, t_own, lam, T=1, C=C)
                lhs_c = coeff * abs(J)                      # D |J| = lhs_c * U^{-1/2}
                vals = tuple(Fr(pt[c]) for c in red.coords)
                rat, rad = dens.evaluate(vals)              # g = rat * sqrt(rad)
                same = ((lhs_c > 0) - (lhs_c < 0) == (rat > 0) - (rat < 0)) and lhs_c ** 2 / U == rat ** 2 * rad
                ok &= same
            Kval = Ksym.subs({s: int(l[0]) for s, l in zip(lsyms, lam)})
            Ks = Knum.sympy()
            print("  lambda = %s" % [int(l[0]) for l in lam])
            print("    B2  pointwise g == D_brute |J| at 4 rational points: %s" % ok)
            print("    B3  K (engine, numeric sources) = %s ; K_sym(lambda) = %s ; equal: %s  (%.1fs)" % (
                Knum, Kval, sp.simplify(Ks - Kval) == 0, time.time() - t0))
            # B4: QMC
            t1 = time.time()
            lamf = np.array([[float(l[0])] for l in lam])
            ex = dens.exponents()
            alpha = {b: max(0.5, float(ex.get(b, 0)) + 1) for b in rm.names}
            mean, err, _ = hn.qmc_integral(G.nv, G.edges, C, plan, lamf, T=1, m=qmc_m, nrep=nrep, seed=1000 + 100 * k + 10 * trial, alpha=alpha)
            exact = float(sp.N(Kval))
            rel = abs(mean - exact) / abs(exact) if exact != 0 else float("nan")
            print("    B4  QMC (2^%d x %d scrambled Sobol, Dirichlet(alpha_b+1) blocks): %.8g +- %.2g ; exact %.10g ; |diff| = %.2g ; rel = %.2g  (%.1fs)" % (
                qmc_m, nrep, mean, err, exact, abs(mean - exact), rel, time.time() - t1))
            t1 = time.time()
            mean, err, _ = hn.qmc_integral(G.nv, G.edges, C, plan, lamf, T=1, m=14, nrep=32, seed=5000 + 100 * k + 10 * trial, alpha=None)
            rel = abs(mean - exact) / abs(exact) if exact != 0 else float("nan")
            print("    B4' QMC (2^14 x 32 scrambled Sobol, Dirichlet(1/2) blocks, no engine input): %.8g +- %.2g ; |diff| = %.2g ; rel = %.2g  (%.1fs)" % (
                mean, err, abs(mean - exact), rel, time.time() - t1))


if __name__ == "__main__":
    which = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else range(5)
    main(which)
