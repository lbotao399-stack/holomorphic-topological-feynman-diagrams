"""t1_rho_square.py -- rho_n ^ rho_n = 0, exact.

 (G)  every admissible graph, D = H+T <= 4, n <= 3 (2 random rational t each) and D <= 3, n = 4
      (1 point each): Pf ^ Pf = 0 as a form; odd n: Pf = 0 and the spanning-tree fiber integral
      (independent of the Pfaffian) = 0.
 (P)  universality on the cone P_n (coordinates m_ij, i <= j): Pf ^ Pf = 0 for n = 2, 4 at random rational M;
      tr (X^X)^j = 0 for X = M^{-1} dM, n = 2, 4, 6, j = 1..3.
 (C)  non-graphic subspaces: random integer C (n x N, rank n), M = C t C^T, n = 2, 4: Pf ^ Pf = 0.
 (H)  hafnian-involution identity (Balduf-Gaiotto, circuit form): for symmetric K of size 2L, L even,
      Q_E = sum_{E = S1 u S2} sgn(S1,S2) det K_{S1 S2} Haf K_{S1} Haf K_{S2} = 0:
      symbolic (sympy) L = 2, 4; random rational K, L = 6.
"""
import random, sys, time
from fractions import Fraction as Fr
from itertools import combinations
from multiprocessing import Pool
import tb
from tb import TopoData, wedge, fadd, fscale, one_form, pfaffian, inverse, mask_of, _wsign, hafnian, det
import graphs as gr


def graph_job(args):
    nv, edges, name, npts, seed = args
    G = tb.Graph(nv, edges, name)
    rng = random.Random(seed)
    out = []
    for _ in range(npts):
        t = tb.rand_t(G.ne, rng)
        D = TopoData(G, t)
        pf = D.rho_pfaffian()
        sq = wedge(pf, pf)
        if D.n % 2:
            gk = D.gkw_fiber()
            out.append((len(pf) == 0 and len(gk) == 0 and len(sq) == 0, len(pf)))
        else:
            out.append((len(sq) == 0, len(pf)))
    return name, out


def check_graphs(pool):
    print("(G) graphs")
    allok = True
    jobs = []
    for D in (1, 2, 3, 4):
        for n in (1, 2, 3):
            for G in gr.admissible(D, n):
                jobs.append(((D, n), (G.nv, G.edges, G.name, 2, hash((D, n, G.name)) % 10 ** 9)))
    for D in (1, 2, 3):
        for G in gr.admissible(D, 4):
            jobs.append(((D, 4), (G.nv, G.edges, G.name, 1, hash((D, 4, G.name)) % 10 ** 9)))
    res = pool.map(graph_job, [j[1] for j in jobs])
    summary = {}
    for (key, _), (name, out) in zip(jobs, res):
        s = summary.setdefault(key, [0, 0, 0])
        s[0] += 1
        s[1] += sum(1 for o in out if o[0])
        s[2] = max(s[2], max(o[1] for o in out))
        if not all(o[0] for o in out):
            allok = False
            print("   FAIL", key, name, out)
    for key in sorted(summary):
        s = summary[key]
        print("   D=%d n=%d: %3d graphs, %3d/%3d points with rho^rho = 0 (odd n: rho = 0 and fiber integral = 0); max #components of Pf: %d"
              % (key[0], key[1], s[0], s[1], s[1] if s[1] else 0, s[2]))
    return allok


def pn_forms(n, rng):
    idx = {}
    for i in range(n):
        for j in range(i, n):
            idx[(i, j)] = len(idx)
    M = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Fr(rng.randint(-5, 5), rng.randint(1, 4))
            M[i][j] = M[j][i] = v
        M[i][i] += 3 * n
    dM = [[{1 << idx[(min(i, j), max(i, j))]: Fr(1)} for j in range(n)] for i in range(n)]
    return M, dM, len(idx)


def check_pn(rng):
    print("(P) cone of positive matrices")
    ok = True
    for n in (2, 4):
        M, dM, dim = pn_forms(n, rng)
        Mi = inverse(M)
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                acc = {}
                for a in range(n):
                    for b in range(n):
                        acc = fadd(acc, wedge(dM[r][a], dM[s][b]), -Mi[a][b] / 2)
                A[r][s] = acc
                A[s][r] = fscale(acc, -1)
        pf = pfaffian(A)
        sq = wedge(pf, pf)
        print("   n=%d (dim P_n = %d): #components Pf = %d, Pf ^ Pf = 0: %s" % (n, dim, len(pf), len(sq) == 0))
        ok &= len(sq) == 0
    for n in (2, 4, 6):
        M, dM, dim = pn_forms(n, rng)
        Mi = inverse(M)
        X = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                acc = {}
                for k in range(n):
                    acc = fadd(acc, dM[k][j], Mi[i][k])
                X[i][j] = acc
        # Y = X ^ X
        Y = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                acc = {}
                for k in range(n):
                    acc = fadd(acc, wedge(X[i][k], X[k][j]))
                Y[i][j] = acc
        Yp = Y
        trs = []
        for jpow in range(1, 4):
            tr = {}
            for i in range(n):
                tr = fadd(tr, Yp[i][i])
            trs.append(len(tr) == 0)
            if jpow < 3:
                Yp = _matmul_forms(Yp, Y, n)
        print("   n=%d: tr (X^X)^j = 0, j=1,2,3: %s" % (n, trs))
        ok &= all(trs)
    return ok


def _matmul_forms(A, Bm, n):
    return [[_dot(A, Bm, i, j, n) for j in range(n)] for i in range(n)]


def _dot(A, Bm, i, j, n):
    acc = {}
    for k in range(n):
        acc = fadd(acc, wedge(A[i][k], Bm[k][j]))
    return acc


def _power(Y, p):
    out = Y
    n = len(Y)
    for _ in range(p - 1):
        out = _matmul_forms(out, Y, n)
    return out


def check_nongraphic(rng):
    print("(C) non-graphic cycle matrices (random integer C, entries in [-3,3])")
    ok = True
    for n, N in ((2, 6), (2, 9), (4, 9), (4, 11)):
        for trial in range(2):
            while True:
                C = [[rng.randint(-3, 3) for _ in range(N)] for _ in range(n)]
                if det([[sum(C[i][e] * C[j][e] for e in range(N)) for j in range(n)] for i in range(n)]) != 0:
                    break
            t = tb.rand_t(N, rng)
            M = [[sum(Fr(C[i][e] * C[j][e]) * t[e] for e in range(N)) for j in range(n)] for i in range(n)]
            Mi = inverse(M)
            dM = [[one_form([C[r][e] * C[a][e] for e in range(N)]) for a in range(n)] for r in range(n)]
            A = [[{} for _ in range(n)] for _ in range(n)]
            for r in range(n):
                for s in range(r + 1, n):
                    acc = {}
                    for a in range(n):
                        for b in range(n):
                            acc = fadd(acc, wedge(dM[r][a], dM[s][b]), -Mi[a][b] / 2)
                    A[r][s] = acc
                    A[s][r] = fscale(acc, -1)
            pf = pfaffian(A)
            sq = wedge(pf, pf)
            print("   n=%d N=%d trial %d: #components Pf = %d, Pf ^ Pf = 0: %s" % (n, N, trial, len(pf), len(sq) == 0))
            ok &= len(sq) == 0 and len(pf) > 0
    return ok


def QE(K, L2):
    """sum over ordered splits E = S1 u S2 (|S1| = |S2| = L2/2) of sgn * det K_{S1S2} Haf K_S1 Haf K_S2"""
    E = list(range(L2))
    h = L2 // 2
    tot = 0
    for S1 in combinations(E, h):
        S2 = tuple(e for e in E if e not in S1)
        sgn = _wsign(mask_of(S1), mask_of(S2))
        dK = det([[K[a][b] for b in S2] for a in S1]) if not hasattr(K[0][0], "free_symbols") else _sdet([[K[a][b] for b in S2] for a in S1])
        tot = tot + sgn * dK * _haf([[K[a][b] for b in S1] for a in S1]) * _haf([[K[a][b] for b in S2] for a in S2])
    return tot


def _haf(A):
    n = len(A)
    if n == 0:
        return 1
    tot = 0
    for j in range(1, n):
        idx = [k for k in range(1, n) if k != j]
        tot = tot + A[0][j] * _haf([[A[a][b] for b in idx] for a in idx])
    return tot


def _sdet(A):
    import sympy as sp
    return sp.Matrix(A).det(method="berkowitz")


def check_hafnian(rng):
    import sympy as sp
    print("(H) hafnian-involution identity")
    ok = True
    for L in (2, 4):
        L2 = 2 * L
        k = sp.symbols("k0:%d" % (L2 * L2))
        K = [[k[min(i, j) * L2 + max(i, j)] for j in range(L2)] for i in range(L2)]
        q = sp.expand(QE(K, L2))
        print("   L=%d (2L=%d), symbolic symmetric K: Q_E = %s" % (L, L2, q))
        ok &= (q == 0)
    for L in (6,):
        L2 = 2 * L
        for trial in range(2):
            K = [[Fr(0)] * L2 for _ in range(L2)]
            for i in range(L2):
                for j in range(i, L2):
                    K[i][j] = K[j][i] = Fr(rng.randint(-9, 9), rng.randint(1, 5))
            q = QE(K, L2)
            print("   L=%d, random rational symmetric K, trial %d: Q_E = %s" % (L, trial, q))
            ok &= (q == 0)
    # control: the same sum without the sign is nonzero
    L2 = 4
    K = [[Fr(rng.randint(1, 9)) for _ in range(L2)] for _ in range(L2)]
    K = [[K[min(i, j)][max(i, j)] for j in range(L2)] for i in range(L2)]
    tot = 0
    for S1 in combinations(range(L2), 2):
        S2 = tuple(e for e in range(L2) if e not in S1)
        tot += det([[K[a][b] for b in S2] for a in S1]) * _haf([[K[a][b] for b in S1] for a in S1]) * _haf([[K[a][b] for b in S2] for a in S2])
    print("   control (sign dropped), L=2: %s (nonzero expected)" % tot)
    return ok


def main():
    t0 = time.time()
    print(__doc__)
    rng = random.Random(7)
    with Pool(12) as pool:
        okG = check_graphs(pool)
    okP = check_pn(rng)
    okC = check_nongraphic(rng)
    okH = check_hafnian(rng)
    print("\nsummary: graphs %s, P_n %s, non-graphic %s, hafnian identity %s; time %.1fs" % (okG, okP, okC, okH, time.time() - t0))


if __name__ == "__main__":
    main()
