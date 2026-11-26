"""t1_normalization.py -- alpha_Gamma (GKW) versus rho_n (htk); four independent evaluations, exact.

 (a) Phi_GKW : spanning-tree expansion of the fiber integral of prod_e Theta_e (GKW Higher.tex l.1451-1457,
               (3.53)-(3.54); x = tau/2, s_e = eta_e); fiber integral = (2 sqrt pi)^{-n} U^{-1/2} Phi_GKW.
 (b) Pf      : Pf(-1/2 dM M^{-1} dM)                   (rho_n = U^{-1/2} Pf).
 (c) Wick    : < prod_r psi_r >, position-space currents (S3 eq-ww).
 (d) Haf     : (-2)^{-n/2} sum_S det(C_S) Haf(K_SS) dt_S,  K = C^T M^{-1} C.
 Claims checked:  (b) = (c) = (d);  (a) = (-1)^n det(G)^{-1} (b),  G = (P; C), P B = 1;
 GKW bitriangle formula (3.56) (Higher.tex l.1500-1508):  alpha = (8 pi)^{-1} U^{-3/2} Br  <=>  2 U Phi_GKW = Br.
"""
import random, sys, time
from fractions import Fraction as Fr
import tb
from tb import TopoData, Graph, det_G_sign, fadd, fscale
import graphs as gr

rng = random.Random(20261122)


def check_graph(G, npts=2):
    res = []
    C = G.cycle_matrix()
    dG = det_G_sign(G, C)
    for _ in range(npts):
        t = tb.rand_t(G.ne, rng)
        D = TopoData(G, t, C)
        pf = D.rho_pfaffian()
        wick, U, P = D.topo_form()
        hf = D.rho_hafnian()
        gk = D.gkw_fiber()
        s = Fr((-1) ** D.n) / dG
        ok_bc = (pf == wick)
        ok_bd = (pf == hf)
        ok_a = (gk == fscale(pf, s))
        res.append((ok_bc, ok_bd, ok_a))
    return dG, res


def main():
    t0 = time.time()
    print(__doc__)
    # ---- GKW bitriangle (3.56)
    G = gr.bitriangle_gkw()
    e = {"12": 0, "23": 1, "24": 2, "31": 3, "43": 4}
    allok = True
    for _ in range(4):
        t = tb.rand_t(G.ne, rng)
        T = {k: t[v] for k, v in e.items()}
        D = TopoData(G, t)
        gk = D.gkw_fiber()
        U = D.U
        Ugkw = T["31"] * (T["24"] + T["43"]) + T["12"] * (T["23"] + T["24"] + T["43"]) + T["23"] * (T["24"] + T["31"] + T["43"])

        def dd(a, b):
            m = (1 << e[a]) | (1 << e[b])
            sgn = 1 if e[a] < e[b] else -1
            return {m: Fr(sgn)}
        Br = {}
        for (a, b, c) in [("12", "23", T["24"] + T["43"]), ("31", "23", T["24"] + T["43"]),
                          ("23", "24", T["12"] + T["31"]), ("23", "43", T["12"] + T["31"]),
                          ("12", "24", -T["23"]), ("12", "43", -T["23"]), ("31", "24", -T["23"]), ("31", "43", -T["23"])]:
            Br = fadd(Br, dd(a, b), c)
        ok = (U == Ugkw) and (fscale(gk, 2 * U) == Br)
        okm = (U == Ugkw) and (fscale(gk, -2 * U) == Br)
        print("bitriangle GKW(3.56): U = GKW bracket: %s ; 2 U Phi_GKW = Br: %s ; 2 U Phi_GKW = -Br: %s" % (U == Ugkw, ok, okm))
        allok &= ok
    dG = det_G_sign(G, G.cycle_matrix())
    print("bitriangle: det G = %s, sgn det[B^T;C] > 0 by construction" % dG)
    # ---- all admissible graphs D <= 4, n <= 3 and the D = 2,3 graphs with n = 4, thetas
    print("\n(b)=(c) Pf = Wick, (b)=(d) Pf = Haf, (a) Phi_GKW = (-1)^n det(G)^{-1} Pf;  2 random rational t each")
    tot = 0
    bad = 0
    for D in (1, 2, 3, 4):
        for n in (1, 2, 3):
            Gs = gr.admissible(D, n)
            okall = True
            dGs = set()
            for G in Gs:
                dG, res = check_graph(G)
                dGs.add(dG)
                tot += 1
                if not all(all(r) for r in res):
                    okall = False
                    bad += 1
                    print("   FAIL", G.name, res)
            print("D=%d n=%d: %3d graphs, all three identities: %s, det G values %s" % (D, n, len(Gs), okall, sorted(dGs)))
    for D, n in ((1, 4), (2, 4), (1, 6)):
        Gs = gr.admissible(D, n)
        okall = True
        for G in Gs:
            dG, res = check_graph(G, 1)
            tot += 1
            if not all(all(r) for r in res):
                okall = False
                bad += 1
        print("D=%d n=%d: %3d graphs, all three identities: %s" % (D, n, len(Gs), okall))
    # reversed / permuted / non-Kruskal cycle bases: identity (a) with det G of that basis
    print("\nnon-fundamental bases C' = A C (random unimodular A), graphs D=3 n=2:")
    for G in gr.admissible(3, 2):
        C = G.cycle_matrix()
        A = [[1, 1], [0, 1]] if rng.random() < 0.5 else [[0, 1], [1, 0]]
        C2 = [[sum(A[i][k] * C[k][e] for k in range(2)) for e in range(G.ne)] for i in range(2)]
        t = tb.rand_t(G.ne, rng)
        D = TopoData(G, t, C2)
        dG = det_G_sign(G, C2)
        print("   %s A=%s det G=%s: Phi_GKW = det(G)^{-1} Pf: %s" % (G.name, A, dG, D.gkw_fiber() == fscale(D.rho_pfaffian(), Fr(1) / dG)))
    print("\ntotal graphs %d, failures %d, bitriangle (3.56) %s, time %.1fs" % (tot, bad, allok, time.time() - t0))


if __name__ == "__main__":
    main()
