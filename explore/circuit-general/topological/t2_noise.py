"""t2_noise.py -- chain move before integrating out the chain noise, with current insertions F (exact).

Decomposition (claim): position-space currents of chain i (walk orientation)
    w_ij = r_i w_* + nu_ij,   w_* = current of the merged edge in Gamma',
    nu_i independent of Gamma' and of nu_l (l != i), Gaussian with
    <nu_ij> = Q_ij - mhat_i,   <dnu_ij dnu_ij'> = (1/2) (delta_jj' / t_ij - 1 / R_i);
    psi(lift of r) = psi'_r + C'_{r*} chi_1,   psi(chain_i - chain_1) = chi_i - chi_1,
    chi_i = -t_red w_* dlog r_i + R_i sum_j nu_ij dx_ij.
Check: < F prod_Gamma psi >_Gamma (brute force, Gamma's own Gaussian) = E[ F prod psi ] in the decomposed
variables (Wick engine), for F = 1, w_e, w_e w_f, w_e^2 w_f (chain and non-chain edges), random t, kappa.
Also the three-chain moves of the theta graphs with trivalent root, D = 2, 3, 4.
"""
import random, time
from fractions import Fraction as Fr
import tb
from tb import Graph, TopoData, fadd, fscale
from wick import Dual, dlog, padd, pmul, pconst, pvar, pscale, block_gauss
import graphs as gr
from moves import find_chains, chain_interior, chain_move_check

rng = random.Random(99)


def noise_check(G, p, q, chains, t, kappa, F_edges_list):
    k = len(chains)
    interior = [chain_interior(G, ch, p) for ch in chains]
    allint = set(v for vs in interior for v in vs)
    keep_v = [v for v in range(G.nv) if v not in allint]
    vmap = {v: i for i, v in enumerate(keep_v)}
    ch_e = set(e for ch in chains for e, s in ch)
    keep_e = [e for e in range(G.ne) if e not in ch_e]
    emap = {e: i for i, e in enumerate(keep_e)}
    edges2 = [(vmap[G.edges[e][0]], vmap[G.edges[e][1]]) for e in keep_e] + [(vmap[p], vmap[q])]
    G2 = Graph(len(keep_v), edges2)
    es = len(edges2) - 1
    kap = {v: (Fr(kappa[v]) if v < G.root else Fr(0)) for v in range(G.nv)}
    T = [Dual.var(t[e], e) for e in range(G.ne)]
    R = [sum((T[e] for e, s in ch), Dual(0)) for ch in chains]
    tred = 1 / sum((1 / Ri for Ri in R), Dual(0))
    r = [tred / Ri for Ri in R]
    kred = {vmap[v]: kap[v] for v in keep_v}
    blocks_mean, blocks_cov, xs = [], [], []
    for i, ch in enumerate(chains):
        a = len(ch)
        x = [T[e] / R[i] for e, s in ch]
        xs.append(x)
        kk = [kap[u] for u in interior[i]]
        sig = [sum((x[j].v for j in range(m + 1, a)), Fr(0)) for m in range(a - 1)]
        mh = sum((sig[m] * kk[m] for m in range(a - 1)), Fr(0))
        Q = [sum(kk[:j], Fr(0)) for j in range(a)]
        kred[vmap[p]] += mh
        for m, u in enumerate(interior[i]):
            kred[vmap[q]] += (1 - sig[m]) * kap[u]
        blocks_mean.append([Q[j] - mh for j in range(a)])
        blocks_cov.append([[(Fr(1, 2) * ((1 / Fr(t[ch[j][0]])) if j == jj else 0) - Fr(1, 2) / R[i].v) for jj in range(a)] for j in range(a)])
    kred_list = [kred[v] for v in range(G2.nv - 1)]
    t2 = [Fr(t[e]) for e in keep_e] + [tred.v]
    C2 = G2.cycle_matrix()
    Clift = []
    for row in C2:
        new = [0] * G.ne
        for e in keep_e:
            new[e] = row[emap[e]]
        for e, s in chains[0]:
            new[e] += row[es] * s
        Clift.append(new)
    for i in range(1, k):
        new = [0] * G.ne
        for e, s in chains[i]:
            new[e] += s
        for e, s in chains[0]:
            new[e] -= s
        Clift.append(new)
    D1 = TopoData(G, t, Clift)
    D2 = TopoData(G2, t2, C2)
    wb2 = D2.wbar(kred_list)
    N2 = G2.ne
    gauss, offs = block_gauss([(wb2, [row[:] for row in D2.cov])] + list(zip(blocks_mean, blocks_cov)))
    nu_index = {}
    for i, ch in enumerate(chains):
        for j, (e, s) in enumerate(ch):
            nu_index[e] = (i, j, s, offs[1 + i] + j)
    img = {emap[e]: {1 << e: Fr(1)} for e in keep_e}
    img[es] = tred.form()

    def w_poly(e):
        """position-space current of Gamma edge e (edge orientation) as a degree-1 polynomial"""
        if e in emap:
            return pvar(emap[e], {0: Fr(1)})
        i, j, s, idx = nu_index[e]
        P = padd(pvar(es, {0: r[i].v}), pvar(idx, {0: Fr(1)}))
        return pscale(P, s)
    chi = []
    for i, ch in enumerate(chains):
        P = pvar(es, fscale(dlog(r[i]), -tred.v))
        for j, (e, s) in enumerate(ch):
            P = padd(P, pvar(offs[1 + i] + j, fscale(xs[i][j].form(), R[i].v)))
        chi.append(P)
    prod = pconst({0: Fr(1)})
    for row in C2:
        fac = {}
        for e in range(N2):
            if row[e]:
                fac = padd(fac, pvar(e, img[e]), row[e])
        fac = padd(fac, chi[0], row[es])
        prod = pmul(prod, fac)
    for i in range(1, k):
        prod = pmul(prod, padd(chi[i], chi[0], -1))
    results = []
    for Fe in F_edges_list:
        Fp = pconst({0: Fr(1)})
        for e in Fe:
            Fp = pmul(Fp, w_poly(e))
        rhs = gauss.E(pmul(Fp, prod))
        lhs, _, _ = D1.topo_form(kappa, {tuple(Fe): Fr(1)})
        results.append(lhs == rhs)
    return results


def main():
    t0 = time.time()
    print(__doc__)
    cases = []
    for D_ in (2, 3, 4):
        for n in (2, 3):
            for G in gr.admissible(D_, n):
                for pp in range(G.nv):
                    for qq in range(pp + 1, G.nv):
                        chs = find_chains(G, pp, qq)
                        if len(chs) < 1 or all(len(c) == 1 for c in chs):
                            continue
                        ints = [chain_interior(G, ch, pp) for ch in chs]
                        if any(G.root in vs for vs in ints):
                            continue
                        cases.append((G, pp, qq, chs))
    for pqr in ((1, 2, 2), (1, 3, 3), (2, 2, 3), (1, 4, 4), (2, 3, 4), (3, 3, 3)):
        G = gr.theta(*pqr)
        cases.append((G, 0, G.root, find_chains(G, 0, G.root)))
    nok, ntot = 0, 0
    for (G, pp, qq, chs) in cases[:60]:
        t = tb.rand_t(G.ne, rng)
        kap = [Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(G.nv - 1)]
        che = [e for ch in chs for e, s in ch]
        oth = [e for e in range(G.ne) if e not in che]
        e1 = rng.choice(che)
        e2 = rng.choice(che)
        Fs = [(), (e1,), (e1, e2), (e1, e1, e2)]
        if oth:
            Fs.append((rng.choice(oth), e1))
        res = noise_check(G, pp, qq, chs, t, kap, Fs)
        ntot += len(res)
        nok += sum(res)
        if not all(res):
            print("FAIL", G.name, pp, qq, [len(c) for c in chs], res)
    print("noise decomposition with insertions: %d/%d (graph, move, insertion) cases exact" % (nok, ntot))
    # theta three-chain moves with the integrated formula
    print("\nthree-chain moves of theta graphs (Gamma' = one edge), integrated formula, kappa random / 0:")
    for pqr in ((1, 2, 2), (1, 3, 3), (2, 2, 3), (1, 4, 4), (2, 3, 4), (3, 3, 3)):
        G = gr.theta(*pqr)
        chs = find_chains(G, 0, G.root)
        outs = []
        for zero in (False, True):
            t = tb.rand_t(G.ne, rng)
            kap = [Fr(0) if zero else Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(G.nv - 1)]
            res = chain_move_check(G, 0, G.root, chs, t, kap)
            outs.append(all(v for k_, v in res.items() if not k_.startswith("_")))
        print("   theta%s: kappa != 0: %s, kappa = 0: %s" % (pqr, outs[0], outs[1]))
    print("\ntime %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
