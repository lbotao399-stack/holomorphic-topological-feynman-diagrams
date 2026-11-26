"""(a) rho_n on series-parallel networks.  All checks exact (Fractions).

 A1  recursion == U^{-1/2} Pf(-1/2 dM M^{-1} dM) for C_rec, all plans of the five SP graphs (d=2,n=4),
     thetas, and D=3, n=4 bananas / two-level graphs
 A2  rho_n ^ rho_n = 0: (i) on P_n (n = 4 all of Sym_4; n = 6, 8 on random 2n-planes), (ii) pulled back
     to Schwinger space of the graphs of A1; (iii) det(B^2) = 0, tr B^{2j} = 0 for B = dM M^{-1}
 A3  graph 2 (degrees [2,2,2,4,4,4]), small plan: rho_4 = s U^{-1/2} prod_e (t_e dw_e/(w_e w'_e)) ^ (1/4T) dlog(t_9 t_10 t_11/T^3)
"""
import sys, random, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import hf_brute as hb
from hf_brute import wedge, fadd, merge_sign, pfaff, fdet, finv
from hf_recursion import Recursion

rng = random.Random(2026)


def plan_of(red):
    return [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]


def top_pair(f, g, gens):
    """coefficient of the top form (gens in order) of f ^ g"""
    tot = Fr(0)
    full = tuple(gens)
    for S, v in f.items():
        rest = tuple(x for x in gens if x not in S)
        if rest in g:
            tot += merge_sign(S, rest) * v * g[rest]
    return tot


def graphs_A1():
    out = []
    for k, Gx in enumerate([Gx for Gx in eg.enumerate_tight(2, 4) if eg.is_sp(Gx)]):
        G, _ = htk.from_networkx(Gx, None)
        for prefer in ("large", "small"):
            out.append(("d2n4 graph%d %s" % (k, prefer), G.nv, G.edges, plan_of(spred.sp_reduce(G, prefer=prefer))))
    for pqr in [(1, 2, 2), (1, 3, 3), (2, 2, 3)]:
        G, _ = htk.theta(*pqr)
        out.append(("theta%s" % (pqr,), G.nv, G.edges, plan_of(spred.sp_reduce(G, prefer="large"))))
    # D = 3, n = 4 banana (1,3,3,3,3) and (2,2,3,3,3)
    for lens in [(1, 3, 3, 3, 3), (2, 2, 3, 3, 3)]:
        G = banana_graph(lens)
        out.append(("banana%s" % (lens,), G.nv, G.edges, plan_of(spred.sp_reduce(G, prefer="large"))))
    return out


def banana_graph(lens):
    """two terminals 0 and 1 (root = last vertex label is re-assigned: terminal q = nv-1)"""
    nv = 2 + sum(a - 1 for a in lens)
    p, q = 0, nv - 1
    edges, v = [], 1
    for a in lens:
        prev = p
        for j in range(a):
            nxt = q if j == a - 1 else v
            if j < a - 1:
                v += 1
            edges.append((prev, nxt))
            prev = nxt
    return htk.Graph(nv, edges)


def A1(cases):
    print("A1  recursion versus Pfaffian (C_rec), exact, 2 random rational points each")
    for name, nv, edges, plan in cases:
        rec = Recursion(nv, edges, plan)
        C = rec.cycle_basis()
        Gm = hb.G(nv, edges)
        res = []
        for _ in range(2):
            t = [Fr(rng.randint(1, 30), rng.randint(1, 12)) for _ in edges]
            rho_r, U_r = rec.run(t)
            pf, U = hb.rho_form(Gm, t, C)
            res.append(rho_r == pf and U_r == U)
        print("   %-28s n=%d  equal: %s" % (name, len(C), res))


def A2_universal():
    print("A2(i)  rho_n ^ rho_n on P_n (exact, random rational positive M)")
    for n, dim in [(2, None), (4, None), (6, 12), (8, 16)]:
        # random positive definite rational M
        X = [[Fr(rng.randint(-5, 5)) for _ in range(n)] for _ in range(n)]
        M = [[sum(X[i][k] * X[j][k] for k in range(n)) + (n if i == j else 0) for j in range(n)] for i in range(n)]
        Mi = finv(M)
        pairs = [(i, j) for i in range(n) for j in range(i, n)]
        if dim is None:
            gens = list(range(len(pairs)))
            dM = [[{} for _ in range(n)] for _ in range(n)]
            for g, (i, j) in enumerate(pairs):
                dM[i][j] = {(g,): Fr(1)}
                dM[j][i] = {(g,): Fr(1)}
        else:   # pullback to a random dim-plane: dM = sum_a S_a dxi_a
            gens = list(range(dim))
            Ss = []
            for a in range(dim):
                S = [[Fr(0)] * n for _ in range(n)]
                for (i, j) in pairs:
                    S[i][j] = S[j][i] = Fr(rng.randint(-3, 3))
                Ss.append(S)
            dM = [[{(a,): Ss[a][i][j] for a in range(dim) if Ss[a][i][j] != 0} for j in range(n)] for i in range(n)]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                acc = {}
                for a in range(n):
                    for b in range(n):
                        if Mi[a][b] != 0:
                            acc = fadd(acc, wedge(dM[r][a], dM[s][b]), -Mi[a][b] / 2)
                A[r][s] = acc
                A[s][r] = {k: -v for k, v in acc.items()}
        t0 = time.time()
        P = pfaff(A, list(range(n)))
        if 2 * n <= len(gens):
            if 2 * n == len(gens):
                val = top_pair(P, P, gens)
                nz = 0 if val == 0 else 1
            else:
                W = wedge(P, P)
                nz = len(W)
        else:
            nz = 0
        # B = dM M^{-1}: tr B^{2j} and det(B^2) (B^2 has commuting 2-form entries)
        B = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                acc = {}
                for a in range(n):
                    if Mi[a][j] != 0:
                        acc = fadd(acc, dM[i][a], Mi[a][j])
                B[i][j] = acc
        B2 = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                acc = {}
                for a in range(n):
                    acc = fadd(acc, wedge(B[i][a], B[a][j]))
                B2[i][j] = acc
        trs = []
        Pw = B2
        for jj in range(1, (min(n, 3) if n < 8 else 2) + 1):
            tr = {}
            for i in range(n):
                tr = fadd(tr, Pw[i][i])
            trs.append(len(tr))
            nxt = [[{} for _ in range(n)] for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    acc = {}
                    for a in range(n):
                        acc = fadd(acc, wedge(Pw[i][a], B2[a][j]))
                    nxt[i][j] = acc
            Pw = nxt
        print("   n=%d (%s): #components of rho^rho = %d;  #components of tr(B^2)^j, j=1..%d: %s  (%.1fs)" % (
            n, "Sym_%d" % n if dim is None else "random %d-plane" % dim, nz, len(trs), trs, time.time() - t0))


def A2_graphs(cases):
    print("A2(ii)  rho_n ^ rho_n pulled back to Schwinger space (gauge slice t_{N-1} = 1), exact")
    for name, nv, edges, plan in cases:
        Gm = hb.G(nv, edges)
        C = Gm.cycles()
        n, N = len(C), len(edges)
        if 2 * n > N - 1 or n % 2:
            continue
        t = [Fr(rng.randint(1, 30), rng.randint(1, 12)) for _ in edges]
        gens = list(range(N - 1))
        P, U = hb.rho_form(Gm, t, C, gens)
        W = wedge(P, P)
        print("   %-28s n=%d, N-1=%d: #components of Pf ^ Pf = %d" % (name, n, N - 1, len(W)))


def A3():
    print("A3  graph 2, small plan: rho_4 = s U^{-1/2} prod_m [t_m dw_m1 / (w_m1 w_m2)] ^ (1/(4T)) dlog(t_9 t_10 t_11 / T^3)")
    edges = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)]
    G = htk.Graph(6, edges)
    red = spred.sp_reduce(G, prefer="small")
    plan = plan_of(red)
    print("   plan", [(p, q, ch, sg) for p, q, ch, sg in plan])
    rec = Recursion(6, edges, plan)
    C = rec.cycle_basis()
    Gm = hb.G(6, edges)
    signs = set()
    for _ in range(4):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 12)) for _ in edges]
        rho_r, U = rec.run(t)
        pf, U2 = hb.rho_form(Gm, t, C)
        assert rho_r == pf
        tt, dt = rec.tt, rec.dt
        form = {(): Fr(1)}
        for m in range(3):
            info = rec.local[m]
            w1, w2 = info["w"]
            tm = info["tred"]
            # dw_1 = d(t_m / R_1) as a form in the dt_e
            R1 = info["R"][0]
            ch1 = plan[m][2][0]
            dR1 = {}
            for e in ch1:
                for kk, v in dt[e].items():
                    dR1[kk] = dR1.get(kk, 0) + v
            dw1 = fadd({kk: v / R1 for kk, v in dt[G.ne + m].items()}, {kk: -tm * v / R1 ** 2 for kk, v in dR1.items()})
            form = wedge(form, {kk: tm / (w1 * w2) * v for kk, v in dw1.items()})
        es = [G.ne, G.ne + 1, G.ne + 2]
        T = sum(tt[e] for e in es)
        dT = {}
        for e in es:
            for kk, v in dt[e].items():
                dT[kk] = dT.get(kk, 0) + v
        dlog = {}
        for e in es:
            dlog = fadd(dlog, {kk: v / tt[e] for kk, v in dt[e].items()})
        dlog = fadd(dlog, {kk: -3 * v / T for kk, v in dT.items()})
        form = wedge(form, {kk: v / (4 * T) for kk, v in dlog.items()})
        ratio = set(rho_r.get(kk, 0) / v for kk, v in form.items()) | set(
            (Fr(0) if kk not in form else None) for kk in rho_r if kk not in form)
        signs |= ratio
    print("   ratio rec / formula over all components and 4 points:", signs)


if __name__ == "__main__":
    cases = graphs_A1()
    A1(cases)
    A2_universal()
    A2_graphs(cases)
    A3()
