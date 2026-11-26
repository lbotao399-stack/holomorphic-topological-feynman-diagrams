"""V3 (claim a2): rho_n of the fan in the triangle basis; V4 (claim a3): Gaussian-source lemma
rho_n = U^{-1/2} E_kappa[prod_r psi_r[kappa]], kappa ~ N(0, Lambda/2), checked by brute-force Wick
expansion (Isserlis on kappa-monomials, no Pfaffian sign conventions) on fans and on other graphs."""
import itertools
import random
import sys
import time
from fractions import Fraction as Fr

from vlib import (Graph, fan_graph, fan_cycles, fundamental_cycles, spanning_tree_edges, normalize_cycles,
                  check_cycles, t_of_r, rand_r, A_forms, pfaffian_forms, wedge, fadd, one_form, matchings,
                  claim_gh, currents_mesh, laplacian, fdet, period, ia, ib)

rng = random.Random(11)
NS = [int(x) for x in sys.argv[1:]] or [2, 4, 6]


def triangle_basis(n, flip_if_negative):
    G = fan_graph(n)
    C = []
    for i in range(1, n + 1):
        row = [0] * G.N
        row[ia(i)] = 1
        row[ib(i + 1, n)] = 1
        row[ib(i, n)] = -1
        C.append(row)
    check_cycles(G, C)
    B = G.incidence()
    d = fdet([[B[e][v] for e in range(G.N)] for v in range(len(G.verts))] + C)
    if flip_if_negative and d < 0:
        C[0] = [-x for x in C[0]]
    return G, C, d


# ------------------------------------------------------------------ a2
print("a2: Pf(-1/2 dM M^-1 dM) [triangle basis] vs 2^{-n/2} sum_pi wedge_r C_r")
for n in NS:
    t0 = time.time()
    G, Craw, d = triangle_basis(n, False)
    print("  n=%d: det[B^T; C_triangles] = %s (unflipped)" % (n, d))
    for p in range(2):
        u, v = rand_r(n, rng)
        t = t_of_r(n, u, v)
        A, M = A_forms(G, Craw, t)
        Pf = pfaffian_forms(A, n)
        g, h, P, R, c = claim_gh(n, u, v)
        a = {i: t[ia(i)] for i in range(1, n + 1)}
        b = {i: t[ib(i, n)] for i in range(1, n + 2)}
        alpha = {r: fadd(one_form([Fr(1) / a[r] if e == ia(r) else 0 for e in range(G.N)]),
                         one_form([Fr(1) / b[r] if e == ib(r, n) else 0 for e in range(G.N)]), -1) for r in range(1, n + 1)}
        beta = {r: fadd(one_form([Fr(1) / a[r] if e == ia(r) else 0 for e in range(G.N)]),
                        one_form([Fr(1) / b[r + 1] if e == ib(r + 1, n) else 0 for e in range(G.N)]), -1) for r in range(1, n + 1)}
        rhs = {}
        for m in matchings(list(range(1, n + 1))):
            small = {i for (i, j) in m}
            w = {0: Fr(1)}
            for r in range(1, n + 1):
                f = g if r in small else h
                Cr = fadd({k: f[r] * x for k, x in alpha[r].items()}, {k: f[r + 1] * x for k, x in beta[r].items()}, -1)
                w = wedge(w, Cr)
            rhs = fadd(rhs, w)
        rhs = {k: x / 2 ** (n // 2) for k, x in rhs.items()}
        diff_plus = fadd(Pf, rhs, -1)
        diff_minus = fadd(Pf, rhs, +1)
        # ratio identities in r
        ok_rat = True
        for r in range(1, n + 1):
            vprev = v[r - 1] if r > 1 else Fr(1)
            if r < n:
                ok_rat &= a[r] / b[r] == (1 - u[r]) * vprev / u[r]
                ok_rat &= a[r] / b[r + 1] == (1 - u[r]) * v[r] / (1 - v[r])
            else:
                ok_rat &= a[n] / b[n] == (1 - u[n] - v[n]) * vprev / u[n]
                ok_rat &= a[n] / b[n + 1] == (1 - u[n] - v[n]) / v[n]
        print("    point %d: Pf - RHS = 0: %s;  Pf + RHS = 0: %s;  #terms(Pf) = %d;  ratio identities a_r/b_r, a_r/b_{r+1} in r: %s"
              % (p, diff_plus == {}, diff_minus == {}, len(Pf), ok_rat))
    print("    (%.1f s)" % (time.time() - t0))


# ------------------------------------------------------------------ a3
def wick_expectation(G, C, t):
    """E[prod_r psi_r[kappa]] for kappa ~ N(0, Lambda/2): psi_r[kappa] = sum_v kappa_v omega_{r,v};
    expand the wedge product with kappa-monomials kept symbolic, then E[kappa^m] by Isserlis."""
    nv = len(G.verts)
    n = len(C)
    Lap = laplacian(G, t)
    Cov = [[Lap[i][j] / 2 for j in range(nv)] for i in range(nv)]
    omega = []
    for r in range(n):
        om = []
        for v in range(nv):
            src = [Fr(int(k == v)) for k in range(nv)]
            y = currents_mesh(G, t, src)
            om.append(one_form([C[r][e] * y[e] for e in range(G.N)]))
        omega.append(om)
    # state: dict (monomial tuple sorted, form mask) -> coeff
    state = {((), 0): Fr(1)}
    for r in range(n):
        new = {}
        for (mono, mask), cf in state.items():
            for v in range(nv):
                for mk, x in omega[r][v].items():
                    if mask & mk:
                        continue
                    w = wedge({mask: cf}, {mk: x})
                    for k2, c2 in w.items():
                        key = (tuple(sorted(mono + (v,))), k2)
                        new[key] = new.get(key, 0) + c2
        state = {k: x for k, x in new.items() if x != 0}
    memo = {}

    def moment(mono):
        if mono in memo:
            return memo[mono]
        if not mono:
            return Fr(1)
        a = mono[0]
        rest = mono[1:]
        tot = Fr(0)
        for k in range(len(rest)):
            if Cov[a][rest[k]] != 0:
                tot += Cov[a][rest[k]] * moment(rest[:k] + rest[k + 1:])
        memo[mono] = tot
        return tot
    out = {}
    for (mono, mask), cf in state.items():
        mval = moment(mono)
        if mval:
            out[mask] = out.get(mask, 0) + cf * mval
    return {k: x for k, x in out.items() if x != 0}


def wheel(k):
    """hub = root k+1, rim 1..k; edges rim (i,i+1), (k,1), spokes (i, hub)."""
    hub = k + 1
    E = [(i, i % k + 1) for i in range(1, k + 1)] + [(i, hub) for i in range(1, k + 1)]
    return Graph(range(1, k + 1), hub, E)


def doubled_triangle():
    """vertices 1,2, root 3; edges (1,2), (1,2), (2,3), (1,3): n = 2."""
    return Graph([1, 2], 3, [(1, 2), (1, 2), (2, 3), (1, 3)])


def k4_plus():
    """K4 (root 4) with the edge (1,2) doubled: n = 4? 7 edges, 3 non-root vertices -> n = 4."""
    return Graph([1, 2, 3], 4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (2, 1)])


print("a3: Pf(-1/2 dM M^-1 dM) = E_kappa[prod_r psi_r[kappa]], kappa ~ N(0, Lambda/2) (Wick, exact)")
tests = [("fan2 pathtree", fan_cycles(2, "pathtree")), ("fan2 triangles", fan_cycles(2, "triangles")),
         ("fan4 pathtree", fan_cycles(4, "pathtree"))]
for name, G in [("doubled triangle n=2", doubled_triangle()), ("wheel W4 n=4", wheel(4)), ("K4 + doubled edge n=4", k4_plus())]:
    T = spanning_tree_edges(G, [])
    C = normalize_cycles(G, fundamental_cycles(G, T))
    tests.append((name, (G, C)))
for name, (G, C) in tests:
    t0 = time.time()
    for p in range(2):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(G.N)]
        A, M = A_forms(G, C, t)
        Pf = pfaffian_forms(A, len(C))
        E = wick_expectation(G, C, t)
        ok = fadd(Pf, E, -1) == {}
        # lemma: y[kappa] for kappa = Lambda tau / 2 equals t^{-1} B tau / 2
        nv = len(G.verts)
        tau = [Fr(rng.randint(-9, 9)) for _ in range(nv)]
        Lap = laplacian(G, t)
        kap = [sum(Lap[i][j] * tau[j] for j in range(nv)) / 2 for i in range(nv)]
        y = currents_mesh(G, t, kap)
        B = G.incidence()
        w = [sum(B[e][j] * tau[j] for j in range(nv)) / (2 * Fr(t[e])) for e in range(G.N)]
        print("  %-24s point %d: Pf == E_kappa[prod psi[kappa]]: %s; y[Lambda tau/2] == t^-1 B tau/2: %s (%d terms)"
              % (name, p, ok, y == w, len(Pf)))
    print("    (%.1f s)" % (time.time() - t0))
