"""Exact pointwise checks: C0 currents, C1 step identities, C2 densities (incl. relabellings, reversals,
edge permutations, interior roots), C4 displacements, C5/C1(5) chart, C7 basis sign, C8 relabel sign."""
import sys, time, random
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.dont_write_bytecode = True
from vcore import *

rng = random.Random(20261120)
t0 = time.time()
GRAPHS = [(0, H, pqr) for H in range(1, 7) for pqr in theta_list(H)] + \
         [(1, H, pqr) for H in range(0, 6) for pqr in theta_list(H + 1)]
print("graphs:", len(GRAPHS))


def wedge2_on(a, b, u, v):
    return matvec(a, u) * matvec(b, v) - matvec(a, v) * matvec(b, u)


def form_eval_2_then_1(two_forms, ones, vecs):
    """(sum_k c_k a_k^b_k) ^ ones[0] ^ ... evaluated on vecs; two_forms = [(c, a, b)]."""
    m = len(vecs)
    FV = [[matvec(F, v) for v in vecs] for F in ones]
    tot = Fr(0)
    for k in range(m):
        for l in range(k + 1, m):
            w = sum(c * wedge2_on(a, b, vecs[k], vecs[l]) for c, a, b in two_forms)
            if w == 0:
                continue
            rest = [x for x in range(m) if x not in (k, l)]
            tot += sg(k + l + 3) * w * det([[row[x] for x in rest] for row in FV])
    return tot


ok_all = True
for (T, H, pqr) in GRAPHS:
    th = theta_std(*pqr)
    ne, nv = th.ne, th.nv
    D = H + T
    p, q, r = pqr
    res = {}
    # C3 rows are cycles
    res["C3 cycles"] = all(is_cycle_basis_row(nv, th.edges, row) for row in th.C3())
    # ---------------- C2 density, C0 currents (root N', sources at all non-root vertices)
    okd, okc = True, True
    for _ in range(3):
        t = rand_t(ne, rng)
        lam = rand_sources(nv, th.Np, H, rng)
        bf, U, y = density_bf(th, t, lam, H, T, rng=rng)
        cl, U2 = density_closed(th, t, lam, H, T)
        okd &= (bf == cl and U == U2)
        if H:
            R = Rsum(th, t)
            sinv = sum(1 / R[P] for P in "ABC")
            L = L_t(th, t, lam, H)
            sig = sigma_t(th, t)
            for P in "ABC":
                wP = (1 / R[P]) / sinv
                mhat = [sum(s * Fr(lam[u][i]) for s, u in zip(sig[P], th.interior[P])) for i in range(H)]
                for j, (e, s) in enumerate(th.paths[P]):
                    Qj = [sum(Fr(lam[u][i]) for u in th.interior[P][:j]) for i in range(H)]
                    pred = [wP * L[i] - mhat[i] + Qj[i] for i in range(H)]
                    okc &= all(s * y[e][i] == pred[i] for i in range(H))
    res["C2 density (root N')"] = okd
    res["C0 currents"] = okc
    # ---------------- C2 under relabelling / reversal / edge permutation, root at a random vertex
    okt = True
    for _ in range(2):
        vperm = list(range(nv)); rng.shuffle(vperm)
        eperm = list(range(ne)); rng.shuffle(eperm)
        rev = [e for e in range(ne) if rng.random() < 0.5]
        th2 = transform(th, vperm, eperm, rev)
        root = rng.randrange(nv)
        t = rand_t(ne, rng)
        lam = rand_sources(nv, root, H, rng)
        bf, U, _ = density_bf(th2, t, lam, H, T, rng=rng, root=root)
        cl, _ = density_closed(th2, t, lam, H, T)
        okt &= (bf == cl)
    res["C2 relabel/reverse/permute/any root"] = okt
    # ---------------- C7: sign det[B^T; C3] and pointwise density in the htk basis
    Bf = incidence_full(nv, th.edges)
    rows = [[Bf[e][v] for e in range(ne)] for v in range(nv) if v != th.Np] + th.C3()
    sdet = 1 if det(rows) > 0 else -1
    res["sgn det[B^T;C3] = (-1)^(q+1)"] = (sdet == sg(q + 1))
    import importlib
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
    import htk
    G, _ = htk.graph.theta(*pqr)
    assert G.edges == th.edges
    Ch = G.cycle_matrix()
    t = rand_t(ne, rng)
    lam = rand_sources(nv, th.Np, H, rng)
    bf_h, _, _ = density_bf(th, t, lam, H, T, C=Ch, rng=rng)
    bf_3, _, _ = density_bf(th, t, lam, H, T, rng=rng)
    res["htk-basis density = (-1)^((q+1)D) D[C3]"] = (bf_h == sg((q + 1) * D) * bf_3)
    # ---------------- C8: relabel rho <-> N' : s s' = -1
    oks = True
    for P in "ABC":
        for rho in th.interior[P]:
            perm = list(range(nv)); perm[rho], perm[th.Np] = th.Np, rho
            th3 = transform(th, perm)
            Bf3 = incidence_full(nv, th3.edges)
            rows3 = [[Bf3[e][v] for e in range(ne)] for v in range(nv) if v != nv - 1] + th3.C3()
            s3 = 1 if det(rows3) > 0 else -1
            oks &= (s3 * sdet == -1)
    res["C8 s*s' = -1 (all interior rho)"] = oks
    # ---------------- C4: displacement map, root N' and root at first interior vertex of each path
    if H:
        okx = True
        roots = [th.Np] + [th.interior[P][0] for P in "ABC" if th.interior[P]]
        for root in roots:
            t = rand_t(ne, rng)
            lam = rand_sources(nv, root, H, rng)
            y = kirchhoff(nv, th.edges, t, lam, root=root)
            z = [[Fr(rng.randint(-7, 7), rng.randint(1, 5)) for _ in range(H)] for _ in range(ne)]
            lhs = -sum(y[e][i] * z[e][i] for e in range(ne) for i in range(H))
            R = Rsum(th, t); sinv = sum(1 / R[P] for P in "ABC")
            w = {P: (1 / R[P]) / sinv for P in "ABC"}
            sig = sigma_t(th, t)
            zp = {P: [[s * z[e][i] for i in range(H)] for e, s in th.paths[P]] for P in "ABC"}
            Z = {P: [sum(zz[i] for zz in zp[P]) for i in range(H)] for P in "ABC"}
            W = [sum(w[P] * Z[P][i] for P in "ABC") for i in range(H)]
            X = {th.N: [-x for x in W], th.Np: [Fr(0)] * H}
            for P in "ABC":
                for k, u in enumerate(th.interior[P], start=1):
                    tail = [sum(zz[i] for zz in zp[P][k:]) for i in range(H)]
                    X[u] = [sig[P][k - 1] * (Z[P][i] - W[i]) - tail[i] for i in range(H)]
            rhs = sum(lam[v][i] * (X[v][i] - X[root][i]) for v in range(nv) if v != root for i in range(H))
            okx &= (lhs == rhs)
        res["C4 X_v"] = okx
    # ---------------- C1(5)/C3: chart Jacobian and chart density
    okj, okg = True, True
    for _ in range(2):
        wa, wb = Fr(rng.randint(1, 30), 97), Fr(rng.randint(1, 30), 97)
        w = {"A": wa, "B": wb, "C": 1 - wa - wb}
        xs = {}
        for P in "ABC":
            raw = [Fr(rng.randint(1, 20)) for _ in range(th.n[P])]
            S = sum(raw)
            xs[P] = [a / S for a in raw]
        t = [None] * ne
        for P in "ABC":
            for (e, _), xj in zip(th.paths[P], xs[P]):
                t[e] = xj / w[P]
        # d t / d xi, xi = (w_A, w_B, x^A_2.., x^B_2.., x^C_2..), w_C = 1 - w_A - w_B, x^P_1 = 1 - sum_{j>=2}
        J = []
        for a in "AB":
            v = [Fr(0)] * ne
            for P in "ABC":
                dw = 1 if P == a else (-1 if P == "C" else 0)
                for (e, _), xj in zip(th.paths[P], xs[P]):
                    v[e] = -xj / w[P] ** 2 * dw
            J.append(v)
        for P in "ABC":
            es = [e for e, _ in th.paths[P]]
            for j in range(1, len(es)):
                v = [Fr(0)] * ne
                v[es[j]] = 1 / w[P]
                v[es[0]] = -1 / w[P]
                J.append(v)
        jac = det([t] + J)
        pred = sg(q - 1)
        for P in "ABC":
            pred /= w[P] ** (th.n[P] + 1)
        okj &= (jac == pred)
        lam = rand_sources(nv, th.Np, H, rng)
        bf, U, _ = density_bf(th, t, lam, H, T, rng=rng)
        if T == 0:
            L = [Fr(x) for x in lam[th.N]]
            for P in "ABC":
                for k, u in enumerate(th.interior[P], start=1):
                    sgm = sum(xs[P][k:])
                    L = [a + sgm * Fr(b) for a, b in zip(L, lam[u])]
            gcl = sg(H * (H + 1) // 2 + 1) * Qform(th, lam, H, L)
            for P in "ABC":
                gcl *= w[P] ** (H - th.n[P])
            okg &= (bf * abs(jac) == gcl)
        else:
            # g = (-1)^{H(H+1)/2} (D_theta/2) prod w^{H+1/2-n_P}; bf = sqrt(U) D, U = 1/(w_A w_B w_C)
            g_times_sqrtU = bf * abs(jac)
            gcl_sq = (D_theta(th, lam, H) / 2) ** 2
            for P in "ABC":
                gcl_sq *= w[P] ** (2 * H + 1 - 2 * th.n[P])
            okg &= (g_times_sqrtU ** 2 / U == gcl_sq) and ((g_times_sqrtU > 0) == (sg(H * (H + 1) // 2) * D_theta(th, lam, H) > 0) or D_theta(th, lam, H) == 0)
            okg &= (U == 1 / (w["A"] * w["B"] * w["C"]))
    res["C1(5) chart jacobian"] = okj
    res["C3 chart density"] = okg
    # ---------------- C1 step identities (1), (2: E = (-1)^{H+q} Q), (3), (4)
    t = rand_t(ne, rng)
    lam = rand_sources(nv, th.Np, H, rng)
    y = kirchhoff(nv, th.edges, t, lam) if H else None
    R = Rsum(th, t)
    sinv = sum(1 / R[P] for P in "ABC")
    sig = sigma_t(th, t)
    L = L_t(th, t, lam, H)
    vecs = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(ne)] for _ in range(ne - 1)]
    dR = {P: [Fr(1) if e in [f for f, _ in th.paths[P]] else Fr(0) for e in range(ne)] for P in "ABC"}
    ok1 = True
    if H:
        for P in "ABC":
            es = [e for e, _ in th.paths[P]]
            for v in vecs[:3]:
                for i in range(H):
                    lhs = sum(s * y[e][i] * v[e] for e, s in th.paths[P])
                    dRv = sum(v[e] for e in es)
                    dm = Fr(0)
                    for k, u in enumerate(th.interior[P], start=1):
                        dsig = sum(v[e] for e in es[k:]) / R[P] - sig[P][k - 1] * dRv / R[P]
                        dm += dsig * Fr(lam[u][i])
                    rhs = L[i] / sinv * dRv / R[P] + R[P] * dm
                    ok1 &= (lhs == rhs)
    res["C1(1) phi_P"] = ok1
    if T == 0:
        Lr = [Fr(rng.randint(-9, 9), 7) for _ in range(H)]
        z = [0] * H
        cols = [col(Lr, z, H), col(z, Lr, H)] + [col(lam[u], z, H) for u in th.interior["A"]] + \
               [col([-x for x in lam[u]], lam[u], H) for u in th.interior["B"]] + \
               [col(z, [-x for x in lam[u]], H) for u in th.interior["C"]]
        E = det([[c[i] for c in cols] for i in range(2 * H)])
        res["C1(2) E = (-1)^(H+q) Q"] = (E == sg(H + q) * Qform(th, lam, H, Lr))
    # (3) iota_E vol = (-1)^{q-1} prod R^{n-1} Omega_R ^ dx
    ones = []
    for P in "ABC":
        es = [e for e, _ in th.paths[P]]
        for j in range(1, len(es)):
            xj = Fr(t[es[j]]) / R[P]
            v = [Fr(0)] * ne
            v[es[j]] += 1 / R[P]
            for e in es:
                v[e] -= xj / R[P]
            ones.append(v)
    OmR = [(R["A"], dR["B"], dR["C"]), (-R["B"], dR["A"], dR["C"]), (R["C"], dR["A"], dR["B"])]
    lhs3 = det([list(map(Fr, t))] + vecs)
    rhs3 = sg(q - 1) * form_eval_2_then_1(OmR, ones, vecs)
    for P in "ABC":
        rhs3 *= R[P] ** (th.n[P] - 1)
    res["C1(3) iota_E vol"] = (lhs3 == rhs3)
    # (4) A12 = Omega_R / (2U) on random pairs (rho_2 = U^{-1/2} A12)
    C = th.C3()
    M = [[sum(Fr(C[a][e] * C[b][e]) * Fr(t[e]) for e in range(ne)) for b in range(2)] for a in range(2)]
    U = det(M)
    Minv = [[M[1][1] / U, -M[0][1] / U], [-M[1][0] / U, M[0][0] / U]]
    dM = [[[Fr(C[a][e] * C[b][e]) for e in range(ne)] for b in range(2)] for a in range(2)]
    ok4 = True
    ex = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(ne)] for _ in range(3)]
    for (u, v) in [(ex[0], ex[1]), (ex[1], ex[2]), (ex[0], ex[2])]:
        A12 = -sum(Minv[a][b] * wedge2_on(dM[0][a], dM[b][1], u, v) for a in range(2) for b in range(2)) / 2
        om = sum(c * wedge2_on(a, b, u, v) for c, a, b in OmR)
        ok4 &= (A12 == om / (2 * U))
    res["C1(4) rho_2 = Omega_R/(2U^{3/2})"] = ok4
    allok = all(res.values())
    ok_all &= allok
    print("T=%d H=%d %s  " % (T, H, pqr) + "  ".join("%s:%s" % (k, v) for k, v in res.items()), flush=True)

# ---------------- C7 general sign formula on a wider family (not only tight graphs)
okw, cnt = True, 0
for p in range(1, 8):
    for q in range(1, 8):
        for r in range(1, 8):
            th = theta_std(p, q, r)
            Bf = incidence_full(th.nv, th.edges)
            rows = [[Bf[e][v] for e in range(th.ne)] for v in range(th.nv) if v != th.Np] + th.C3()
            d = det(rows)
            okw &= ((1 if d > 0 else -1) == sg(q + 1))
            cnt += 1
th = theta_std(1, 1, 1)
Bf = incidence_full(th.nv, th.edges)
d = det([[Bf[e][v] for e in range(th.ne)] for v in range(th.nv) if v != th.Np] + th.C3())
print("C7 sign (-1)^(q+1) on %d graphs theta(p,q,r), 1<=p,q,r<=7: %s ; theta(1,1,1): %s" % (cnt, okw, (1 if d > 0 else -1) == sg(2)))
print("ALL POINTWISE CHECKS PASS:", ok_all and okw, " time %.1f s" % (time.time() - t0))
