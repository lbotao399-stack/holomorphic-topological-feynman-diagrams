"""t1_vanishing.py -- the full reduced integrand at general (H, T), exact at rational t:

    Omega = U^{-H} prod_{a=1..T} Tf^a[kappa^a] ^ prod_{i=1..H} prod_r psi^i_r = U^{-H-T/2} e^{sum_a P_a} Om,
    Om = Phi^1 ^ ... ^ Phi^T ^ prod_{i,r} psi^i_r  (degree n(H+T) = N-1),  D(t) = Om_{hat 0} / t_0 (times U, e^P factors).

 (V1) kappa = 0, T >= 2: Om = 0 for every admissible graph with D = H+T <= 4, n <= 3 (all splittings H+T = D).
 (V2) kappa = 0, T = 1, n odd: Om = 0.
 (V3) controls: T = 1, n even, kappa = 0: Om != 0;  T >= 1 with kappa != 0 (n odd at T = 1; T = 2): Om != 0 in general.
"""
import random, time
from fractions import Fraction as Fr
import tb
from tb import TopoData, wedge_all, one_form, currents
import graphs as gr

rng = random.Random(1618)


def full_Om(G, t, H, T, lam, kappas):
    D = TopoData(G, t)
    n = D.n
    C = D.C
    forms = [D.topo_form(kappas[a])[0] for a in range(T)]
    ys = [currents(G, D.t, [lam[v][i] for v in range(G.nv - 1)], D.Linv) for i in range(H)]
    psis = [one_form([C[r][e] * ys[i][e] for e in range(G.ne)]) for i in range(H) for r in range(n)]
    Om = wedge_all(forms + psis)
    return Om


def main():
    t0 = time.time()
    print(__doc__)
    stats = {}
    ok = True
    for D_ in (1, 2, 3, 4):
        for n in (1, 2, 3):
            for G in gr.admissible(D_, n):
                for T in range(1, D_ + 1):
                    H = D_ - T
                    for withk in (False, True):
                        t = tb.rand_t(G.ne, rng)
                        lam = [[Fr(rng.randint(-5, 5)) for _ in range(H)] for _ in range(G.nv - 1)]
                        kap = [[Fr(rng.randint(-4, 4), rng.randint(1, 3)) if withk else Fr(0) for _ in range(G.nv - 1)] for _ in range(T)]
                        Om = full_Om(G, t, H, T, lam, kap)
                        key = (D_, n, H, T, "kappa" if withk else "kappa=0")
                        s = stats.setdefault(key, [0, 0])
                        s[0] += 1
                        s[1] += (len(Om) > 0)
                        if not withk and (T >= 2 or n % 2 == 1) and Om:
                            ok = False
                            print("FAIL", key, G.name)
    print("(D, n, H, T, sources): graphs, number with Om != 0")
    for key in sorted(stats):
        D_, n, H, T, k = key
        tag = ""
        if k == "kappa=0":
            tag = "  [expected 0]" if (T >= 2 or n % 2) else "  [expected all]"
        print("   D=%d n=%d H=%d T=%d %-8s: %3d graphs, Om != 0 for %3d%s" % (D_, n, H, T, k, stats[key][0], stats[key][1], tag))
    print("\nV1, V2 hold: %s ; time %.1fs" % (ok, time.time() - t0))


if __name__ == "__main__":
    main()
