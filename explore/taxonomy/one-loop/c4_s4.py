"""(S4) T = 1: odd-loop kernels vanish for every H; topological sources kappa only enter with positive weight.

Topological factor (one real direction), CONVENTIONS sec. 2 / topological REPORT sec. 0:
   T_Gamma[kappa] = U^{-1/2} e^{P} < prod_{r=1}^n psi_r >,  psi_r = sum_e C_re w_e dt_e,
   <w_e> = wbar_e (Kirchhoff currents of the sources kappa), <dw_e dw_f> = (delta_ef/t_e - K_ef)/2, K = C^T M^-1 C,
   P = sum_e t_e wbar_e^2.
Coefficient of dt_S (S an n-set of edges, increasing): G_S = det C[:,S] * E[prod_{f in S} w_f] (Isserlis with mean).
Checks, exact Fractions at random rational (t, kappa):
 (1) kappa = 0, n even: G_S = [Pf(-1/2 dM M^-1 dM)]_S  (Isserlis route = Pfaffian route; S3 prop-rho);
 (2) kappa = 0, n odd: G = 0;   kappa != 0, n odd: G != 0 (control);
 (3) scaling: P(4t;k) = P(t;2k) and 2^n G_S(4t;k) = G_S(t;2k); kappa-degree-j part: G^(j)(4t) = 2^(j-n) G^(j)(t),
     i.e. phi_c^* (U^{-1/2} G^(j) dt_S) = c^{j/2} (U^{-1/2} G^(j) dt_S): only j = 0 has weight 0;
 (4) T = 2: product of two directions with independent sources, kappa-bidegree (0,0) = Pf^Pf = 0, total != 0.
Graphs: (H,1), H = 1..6 (D = H+1): n = 1 polygons; n = 3 all admissible graphs D = 2..6; n = 5 at D = 2 (htk
catalog); even controls n = 2 (D = 2..7), n = 4 (D = 2, 3).
"""
import sys
import random
import time
from fractions import Fraction as Fr
from itertools import combinations

from common import Log, det, inverse, wedge, pfaffian, fadd, P as PROJ
import c4graphs as cg

sys.path.insert(0, str(PROJ / "code"))
from htk import catalog  # noqa: E402

log = Log("c4_s4.log")
rng = random.Random(4242)


def data(nv, E, t, kap):
    """wbar, P, covariance Sigma, cycle matrix C, M"""
    Ne = len(E)
    B = cg.incidence(nv, E)
    m = nv - 1
    Lam = [[sum(Fr(B[e][u] * B[e][v]) / t[e] for e in range(Ne)) for v in range(m)] for u in range(m)]
    Li = inverse(Lam)
    phi = [sum(Li[u][v] * kap[v] for v in range(m)) for u in range(m)]
    wb = [sum(B[e][v] * phi[v] for v in range(m)) / t[e] for e in range(Ne)]
    C = cg.cycle_matrix(nv, E)
    n = len(C)
    assert all(sum(B[e][v] * wb[e] for e in range(Ne)) == kap[v] for v in range(m))          # KCL
    assert all(sum(C[r][e] * t[e] * wb[e] for e in range(Ne)) == 0 for r in range(n))       # KVL
    Pv = sum(t[e] * wb[e] ** 2 for e in range(Ne))
    assert Pv == sum(kap[u] * phi[u] for u in range(m))
    M = [[sum(C[a][e] * C[b][e] * t[e] for e in range(Ne)) for b in range(n)] for a in range(n)]
    Mi = inverse(M)
    K = [[sum(C[a][e] * Mi[a][b] * C[b][f] for a in range(n) for b in range(n)) for f in range(Ne)] for e in range(Ne)]
    Sig = [[(Fr(1) / t[e] if e == f else 0) * Fr(1, 2) - K[e][f] / 2 for f in range(Ne)] for e in range(Ne)]
    return wb, Pv, Sig, C, M, Mi


def isserlis_by_degree(S, wb, Sig):
    """E[prod_{f in S} w_f] split by number j of unmatched (mean) factors: dict j -> value"""
    memo = {}

    def rec(T):
        if not T:
            return {0: Fr(1)}
        if T in memo:
            return memo[T]
        a, rest = T[0], T[1:]
        out = {}
        for j, v in rec(rest).items():
            out[j + 1] = out.get(j + 1, 0) + wb[a] * v
        for i, b in enumerate(rest):
            if Sig[a][b]:
                sub = rest[:i] + rest[i + 1:]
                for j, v in rec(sub).items():
                    out[j] = out.get(j, 0) + Sig[a][b] * v
        memo[T] = out
        return out
    return rec(tuple(S))


def G_components(nv, E, t, kap):
    wb, Pv, Sig, C, M, Mi = data(nv, E, t, kap)
    n = len(C)
    comps = {}
    for S in combinations(range(len(E)), n):
        dC = det([[C[r][f] for f in S] for r in range(n)])
        if dC == 0:
            continue
        for j, v in isserlis_by_degree(S, wb, Sig).items():
            if v:
                comps.setdefault(j, {})[S] = dC * v
    return comps, Pv, C, M, Mi


def total(comps):
    out = {}
    for j, d in comps.items():
        for S, v in d.items():
            out[S] = out.get(S, 0) + v
    return {S: v for S, v in out.items() if v != 0}


def pf_route(nv, E, t, C, M, Mi):
    n = len(C)
    Ne = len(E)
    dM = [[{1 << e: C[a][e] * C[b][e] for e in range(Ne) if C[a][e] * C[b][e]} for b in range(n)] for a in range(n)]
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Mi[a][b] and dM[r][a] and dM[b][s]:
                        acc = fadd(acc, wedge(dM[r][a], dM[b][s]), Mi[a][b])
            A[r][s] = {k: -v / 2 for k, v in acc.items()}
    Pf = pfaffian(A)
    out = {}
    for mask, v in Pf.items():
        S = tuple(e for e in range(Ne) if mask >> e & 1)
        out[S] = v
    return out, Pf


def check_graph(nv, E):
    Ne = len(E)
    m = nv - 1
    t = [Fr(rng.randint(1, 40), rng.randint(1, 6)) for _ in range(Ne)]
    kap = [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(m)]
    if all(k == 0 for k in kap):
        kap[0] = Fr(1)
    zero = [Fr(0)] * m
    c0, _, C, M, Mi = G_components(nv, E, t, zero)
    n = len(C)
    G0 = total(c0)
    res = {}
    if n % 2 == 0:
        pfr, _ = pf_route(nv, E, t, C, M, Mi)
        res["pf_route"] = (G0 == {S: v for S, v in pfr.items() if v != 0}) and len(G0) > 0
    else:
        res["k0_zero"] = (len(G0) == 0)
    ck, Pk, _, _, _ = G_components(nv, E, t, kap)
    Gk = total(ck)
    res["k_nonzero"] = len(Gk) > 0
    t4 = [4 * x for x in t]
    c4, P4, _, _, _ = G_components(nv, E, t4, kap)
    c2, P2, _, _, _ = G_components(nv, E, t, [2 * k for k in kap])
    G4, G2 = total(c4), total(c2)
    res["P_scaling"] = (P4 == P2)
    res["G_scaling"] = ({S: 2 ** n * v for S, v in G4.items()} == G2)
    okj = True
    for j in set(ck) | set(c4):
        a = c4.get(j, {})
        b = ck.get(j, {})
        okj &= ({S: v for S, v in a.items()} == {S: Fr(2) ** (j - n) * v for S, v in b.items()})
    res["degree_weights"] = okj
    res["degrees"] = sorted(ck)
    return n, res


def run(label, graphs, H):
    t0 = time.time()
    agg = {}
    degs = set()
    for nv, E in graphs:
        n, res = check_graph(nv, E)
        for k, v in res.items():
            if k == "degrees":
                degs.update(v)
                continue
            agg[k] = agg.get(k, True) and v
    log("  %-34s H=%d T=1: %3d graphs  %s  kappa-degrees %s  (%.1fs)"
        % (label, H, len(graphs), "  ".join("%s=%s" % kv for kv in agg.items()), sorted(degs), time.time() - t0))
    return agg


log("odd loops (S4) and even controls; (H,1), D = H+1")
allok = True
for H in range(1, 7):
    D = H + 1
    a = run("n=1 polygon C_%d" % (D + 1), cg.admissible(D, 1), H)
    allok &= all(a.values())
for D in range(2, 7):
    gs = cg.admissible(D, 3)
    a = run("n=3 all admissible (D=%d)" % D, gs, D - 1)
    allok &= all(a.values())
gs5 = [(G.nv, [tuple(e) for e in G.edges]) for G in catalog.tight_graphs(2, 5)]
a = run("n=5 all admissible (D=2, htk catalog)", gs5, 1)
allok &= all(a.values())
for D in range(2, 8):
    a = run("n=2 thetas (D=%d)" % D, cg.admissible(D, 2), D - 1)
    allok &= all(a.values())
gs4 = [(G.nv, [tuple(e) for e in G.edges]) for G in catalog.tight_graphs(2, 4)]
a = run("n=4 all admissible (D=2, htk catalog)", gs4, 1)
allok &= all(a.values())
gs43 = cg.admissible(3, 4)
a = run("n=4 all admissible (D=3)", gs43, 2)
allok &= all(a.values())
log("  admissible counts n=3, D=2..6:", [len(cg.admissible(D, 3)) for D in range(2, 7)],
    "; n=4 D=3:", len(gs43), "; n=5 D=2:", len(gs5))


# cubic-skeleton subcount at D = 5, n = 3 (atlas: bulk3(5,3) = 83 graphs)
def max_deg(nv, E):
    d = [0] * nv
    for x, y in E:
        d[x] += 1
        d[y] += 1
    return max(d)


log("  D=5 n=3 admissible with max degree 3:", sum(1 for nv, E in cg.admissible(5, 3) if max_deg(nv, E) == 3),
    "(atlas REPORT bulk3(5,3): 83)")

log("T = 2 (S3 with sources): two directions, independent kappa^1, kappa^2")
ok2 = True
cnt = 0
for D in (2, 3, 4, 5):
    for nv, E in cg.admissible(D, 2)[:3]:
        Ne = len(E)
        m = nv - 1
        t = [Fr(rng.randint(1, 40), rng.randint(1, 6)) for _ in range(Ne)]
        k1 = [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(m)]
        k2 = [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(m)]
        cA, _, _, _, _ = G_components(nv, E, t, k1)
        cB, _, _, _, _ = G_components(nv, E, t, k2)

        def form(d):
            return {sum(1 << e for e in S): v for S, v in d.items()}
        prod00 = wedge(form(cA.get(0, {})), form(cB.get(0, {})))
        prodall = wedge(form(total(cA)), form(total(cB)))
        ok2 &= (len(prod00) == 0) and (len(prodall) > 0)
        cnt += 1
log("  n=2 graphs at D = 2..5 (%d): bidegree (0,0) = Pf^Pf = 0 and full product != 0: %s" % (cnt, ok2))
log("all checks:", allok and ok2)
log.close()
