"""moves.py -- the local moves acting on the topological factor (one direction), with sources kappa.

Notation (CONVENTIONS sec 5): chain move with k chains p -> q, chain i has edges (i,j), j = 1..a_i
(walk orientation from p), interior vertices u_{i,1..a_i-1} with sources kappa_{i,k};
  R_i = sum_j t_ij, x_ij = t_ij/R_i, 1/t_red = sum_i 1/R_i, r_i = t_red/R_i,
  sigma_ik = sum_{j>k} x_ij, mhat_i = sum_k sigma_ik kappa_ik, Q_ij = sum_{k<j} kappa_ik.
Claims (checked here):
  kappa_red_p = kappa_p + sum_i mhat_i, kappa_red_q = kappa_q + sum_{i,k} (1 - sigma_ik) kappa_ik;
  wbar_ij = r_i wbar_red - mhat_i + Q_ij;   P = P_red + sum_i R_i Var_{x_i}(Q_i);   U = (prod_i R_i / t_red) U_red;
  < prod_Gamma psi > = < prod psi' ^ bd(chibar_1..chibar_k) > + (-1)^{n'-1} < d_zeta(prod psi') ^ chibar_1...chibar_k >,
  chibar_i = R_i (d mhat_i - w_* d r_i),  bd(c_1..c_k) = sum_j (-1)^{j-1} c_1..^c_j..c_k,
  basis of Gamma: lifts (through chain 1) of the basis of Gamma', then chain_i - chain_1, i = 2..k.
"""
from fractions import Fraction as Fr
import tb
from tb import Graph, TopoData, fadd, fscale, wedge, pullback
from wick import Dual, dlog, padd, pmul, pconst, pvar, pscale, Gauss


def find_chains(G, p, q):
    adj = {v: [] for v in range(G.nv)}
    for e, (a, b) in enumerate(G.edges):
        adj[a].append((b, e, +1))
        adj[b].append((a, e, -1))
    deg = {v: len(adj[v]) for v in adj}
    chains = []
    used = set()
    for (u, e, s) in adj[p]:
        if e in used:
            continue
        ch = [(e, s)]
        prev, cur = p, u
        path_e = {e}
        ok = True
        while cur not in (p, q):
            if deg[cur] != 2:
                ok = False
                break
            nxt = [(w, f, s2) for (w, f, s2) in adj[cur] if f not in path_e]
            if len(nxt) != 1:
                ok = False
                break
            w, f, s2 = nxt[0]
            ch.append((f, s2))
            path_e.add(f)
            prev, cur = cur, w
        if ok and cur == q:
            chains.append(ch)
            used |= path_e
    return chains


def chain_interior(G, ch, p):
    """interior vertices of a chain in walk order"""
    vs = []
    cur = p
    for e, s in ch[:-1]:
        a, b = G.edges[e]
        cur = b if s == +1 else a
        vs.append(cur)
    return vs


def chain_move_check(G, p, q, chains, t, kappa):
    """kappa: list over the non-root vertices of G.  Returns dict of booleans."""
    k = len(chains)
    interior = [chain_interior(G, ch, p) for ch in chains]
    allint = set(v for vs in interior for v in vs)
    assert G.root not in allint
    keep_v = [v for v in range(G.nv) if v not in allint]
    vmap = {v: i for i, v in enumerate(keep_v)}
    ch_e = set(e for ch in chains for e, s in ch)
    keep_e = [e for e in range(G.ne) if e not in ch_e]
    emap = {e: i for i, e in enumerate(keep_e)}
    edges2 = [(vmap[G.edges[e][0]], vmap[G.edges[e][1]]) for e in keep_e] + [(vmap[p], vmap[q])]
    G2 = Graph(len(keep_v), edges2, name=G.name + "_red")
    es = len(edges2) - 1
    kap = {v: (Fr(kappa[v]) if v < G.root else Fr(0)) for v in range(G.nv)}
    # ratio coordinates as exact duals in the coordinates t_e of Gamma
    T = [Dual.var(t[e], e) for e in range(G.ne)]
    R = [sum((T[e] for e, s in ch), Dual(0)) for ch in chains]
    tred = 1 / sum((1 / Ri for Ri in R), Dual(0))
    r = [tred / Ri for Ri in R]
    mhat, Qs, varloc = [], [], Fr(0)
    for i, ch in enumerate(chains):
        a = len(ch)
        x = [T[e] / R[i] for e, s in ch]
        kk = [kap[u] for u in interior[i]]                      # kappa_{i,1..a-1}
        sig = [sum((x[j] for j in range(kk_ + 1, a)), Dual(0)) for kk_ in range(a - 1)]   # sigma_{i,k}, k = 1..a-1
        mh = sum((sig[m] * kk[m] for m in range(a - 1)), Dual(0))
        Q = [sum(kk[:j], Fr(0)) for j in range(a)]             # Q_{i,j}, j = 1..a  (sum_{k<j})
        mhat.append(mh)
        Qs.append((Q, sig))
        varloc += R[i].v * (sum((x[j].v * Q[j] ** 2 for j in range(a)), Fr(0)) - sum((x[j].v * Q[j] for j in range(a)), Fr(0)) ** 2)
    # reduced sources
    kred = {vmap[v]: kap[v] for v in keep_v}
    kred[vmap[p]] += sum((m.v for m in mhat), Fr(0))
    for i, ch in enumerate(chains):
        Q, sig = Qs[i]
        for m, u in enumerate(interior[i]):
            kred[vmap[q]] += (1 - sig[m].v) * kap[u]
    kred_list = [kred[v] for v in range(G2.nv - 1)]
    t2 = [Fr(t[e]) for e in keep_e] + [tred.v]
    # bases
    C2 = G2.cycle_matrix()
    n2 = len(C2)
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
    PhiG, UG, PG = D1.topo_form(kappa)
    D2 = TopoData(G2, t2, C2)
    res = {}
    # classical
    wb1 = D1.wbar(kappa)
    wb2 = D2.wbar(kred_list)
    ok = all(wb1[e] == wb2[emap[e]] for e in keep_e)
    for i, ch in enumerate(chains):
        Q, sig = Qs[i]
        for j, (e, s) in enumerate(ch):
            ok &= (s * wb1[e] == r[i].v * wb2[es] - mhat[i].v + Q[j])
    res["wbar reconstruction"] = ok
    res["kappa transfer (P)"] = (PG == tb.power(G2, t2, wb2) + varloc)
    prodR = Fr(1)
    for Ri in R:
        prodR *= Ri.v
    res["U = prod R / t_red U'"] = (UG == prodR / tred.v * D2.U)
    # forms
    images = {emap[e]: {1 << e: Fr(1)} for e in keep_e}
    images[es] = tred.form()
    chis = []
    for i in range(k):
        alpha = fscale(mhat[i].form(), R[i].v)
        beta = fscale(r[i].form(), -R[i].v)
        P = {}
        if alpha:
            P[()] = alpha
        if beta:
            P[(0,)] = beta
        chis.append(P)
    bd = {}
    for j in range(k):
        prod = pconst({0: Fr(1)})
        for i in range(k):
            if i != j:
                prod = pmul(prod, chis[i])
        bd = padd(bd, prod, 1 if j % 2 == 0 else -1)
    full = pconst({0: Fr(1)})
    for i in range(k):
        full = pmul(full, chis[i])
    rhs = {}
    for mono, form in bd.items():
        F = {tuple([es] * len(mono)): Fr(1)}
        Phi2, _, _ = D2.topo_form(kred_list, F)
        rhs = fadd(rhs, wedge(pullback(Phi2, images), form))
    sgn = -1 if (n2 - 1) % 2 else 1
    for mono, form in full.items():
        F = {tuple([es] * len(mono)): Fr(1)}
        Phi2, _, _ = D2.topo_form_cut(kred_list, F, (es,))
        rhs = fadd(rhs, wedge(pullback(Phi2, images), form), sgn)
    res["<prod psi> = chain formula"] = (PhiG == rhs)
    res["_ncomp"] = len(PhiG)
    return res


# ------------------------------------------------------------------ star-triangle

def star_triangle_check(G, v, u, w, x, a, b, c, d, f, t, kappa):
    """G.edges[a]=(v,u), [b]=(v,w), [c]=(v,x), [d]=(u,x), [f]=(w,x); no edge u-w; v != root."""
    assert G.edges[a] == (v, u) and G.edges[b] == (v, w) and G.edges[c] == (v, x)
    assert G.edges[d] == (u, x) and G.edges[f] == (w, x) and v != G.root
    keep_v = [y for y in range(G.nv) if y != v]
    vmap = {y: i for i, y in enumerate(keep_v)}
    loc = {a, b, c, d, f}
    keep_e = [e for e in range(G.ne) if e not in loc]
    emap = {e: i for i, e in enumerate(keep_e)}
    edges2 = [(vmap[G.edges[e][0]], vmap[G.edges[e][1]]) for e in keep_e] + [(vmap[u], vmap[w]), (vmap[u], vmap[x]), (vmap[w], vmap[x])]
    G2 = Graph(len(keep_v), edges2, name=G.name + "_YD")
    iU, iX, iW = len(edges2) - 3, len(edges2) - 2, len(edges2) - 1
    T = [Dual.var(t[e], e) for e in range(G.ne)]
    ta, tb_, tc, td, tf = T[a], T[b], T[c], T[d], T[f]
    Sig = ta * tb_ + tb_ * tc + tc * ta
    TU, TX, TW = Sig / tc, Sig / tb_, Sig / ta          # triangle uw, ux, wx
    tX = 1 / (1 / td + 1 / TX)
    tW = 1 / (1 / tf + 1 / TW)
    r1, r2 = tX / td, tW / tf
    S = TU + TX + TW
    Lvv = 1 / ta + 1 / tb_ + 1 / tc
    kap = {y: (Fr(kappa[y]) if y < G.root else Fr(0)) for y in range(G.nv)}
    kv = kap[v]
    kred = {vmap[y]: kap[y] for y in keep_v}
    kred[vmap[u]] += kv * (1 / ta).v / Lvv.v
    kred[vmap[w]] += kv * (1 / tb_).v / Lvv.v
    kred[vmap[x]] += kv * (1 / tc).v / Lvv.v
    kred_list = [kred[y] for y in range(G2.nv - 1)]
    t2 = [Fr(t[e]) for e in keep_e] + [TU.v, tX.v, tW.v]
    # Gamma' basis with Delta = U + W - X last
    C2 = G2.cycle_matrix()
    _, cot2 = G2.fundamental_cycles()
    delta = [0] * G2.ne
    delta[iU], delta[iW], delta[iX] = 1, 1, -1
    coeffs = [delta[ch] for ch in cot2]
    r0 = next(i for i, cf in enumerate(coeffs) if cf != 0)
    C2b = [row for i, row in enumerate(C2) if i != r0] + [delta]
    n2 = len(C2b)
    # lifted basis of Gamma: non-Delta rows lifted, then c1 = a + d - c, c2 = b + f - c
    Cg = []
    for row in C2b[:-1]:
        new = [0] * G.ne
        for e in keep_e:
            new[e] = row[emap[e]]
        new[a] += -row[iU]
        new[b] += row[iU]
        new[d] += row[iX]
        new[f] += row[iW]
        Cg.append(new)
    c1 = [0] * G.ne
    c1[a], c1[d], c1[c] = 1, 1, -1
    c2 = [0] * G.ne
    c2[b], c2[f], c2[c] = 1, 1, -1
    Cg += [c1, c2]
    D1 = TopoData(G, t, Cg)
    PhiG, UG, PG = D1.topo_form(kappa)
    D2 = TopoData(G2, t2, C2b)
    res = {}
    # classical
    wb1 = D1.wbar(kappa)
    wb2 = D2.wbar(kred_list)
    epsbar = kv / Lvv.v
    ok = all(wb1[e] == wb2[emap[e]] for e in keep_e)
    ok &= wb1[d] == r1.v * wb2[iX] and wb1[f] == r2.v * wb2[iW]
    ok &= wb1[a] == -wb2[iU] - (1 - r1.v) * wb2[iX] + epsbar / ta.v
    ok &= wb1[b] == wb2[iU] - (1 - r2.v) * wb2[iW] + epsbar / tb_.v
    ok &= wb1[c] == (1 - r1.v) * wb2[iX] + (1 - r2.v) * wb2[iW] + epsbar / tc.v
    res["wbar reconstruction"] = ok
    res["P = P' + kappa_v^2/Lam_vv"] = (PG == tb.power(G2, t2, wb2) + kv ** 2 / Lvv.v)
    res["U = t_c/(r1 r2) U'"] = (UG == tc.v / (r1.v * r2.v) * D2.U)
    # Wick engine: variables = Gamma' edge currents (0..N2-1), eps_v (N2)
    N2 = G2.ne
    cov = [[D2.cov[i][j] for j in range(N2)] + [Fr(0)] for i in range(N2)] + [[Fr(0)] * N2 + [1 / (2 * Lvv.v)]]
    gauss = Gauss(list(wb2) + [epsbar], cov)
    ev = N2
    img = {emap[e]: {1 << e: Fr(1)} for e in keep_e}
    img[iU], img[iX], img[iW] = TU.form(), tX.form(), tW.form()

    def psi2(row):
        P = {}
        for e in range(N2):
            if row[e]:
                P = padd(P, pvar(e, img[e]), row[e])
        return P
    psiD = padd(padd(pvar(iU, TU.form()), pvar(iW, fscale(TW.form(), 1 - r2.v))), pvar(iX, fscale(TX.form(), 1 - r1.v)), -1)
    muU = padd(pscale(psiD, -(TU.v / S.v)), pvar(ev, dlog(TW / TX)))
    muX = pvar(iX, fscale(r1.form(), -tX.v / r1.v))
    muW = pvar(iW, fscale(r2.form(), -tW.v / r2.v))
    pc1 = padd(padd(pscale(psiD, -(TX.v / S.v)), pvar(iX, fscale(r1.form(), -tX.v / (r1.v * (1 - r1.v))))), pvar(ev, dlog(TU / TW)))
    pc2 = padd(padd(pscale(psiD, TW.v / S.v), pvar(iW, fscale(r2.form(), -tW.v / (r2.v * (1 - r2.v))))), pvar(ev, dlog(TU / TX)))
    prod = pconst({0: Fr(1)})
    for row in C2b[:-1]:
        fac = psi2(row)
        fac = padd(fac, muU, row[iU])
        fac = padd(fac, muX, row[iX])
        fac = padd(fac, muW, row[iW])
        prod = pmul(prod, fac)
    prod = pmul(pmul(prod, pc1), pc2)
    rhs = gauss.E(prod)
    res["<prod psi> = star-triangle formula"] = (PhiG == rhs)
    res["_ncomp"] = len(PhiG)
    return res
