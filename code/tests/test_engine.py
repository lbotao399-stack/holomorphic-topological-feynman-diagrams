"""Internal consistency checks of htk (exact arithmetic unless stated).

Run:  python3 code/tests/test_engine.py        (HTK_FAST=1 skips the slow cases)

Statements of htk/CONVENTIONS.md verified here:
  - sigma_T^{H+T} K_T is independent of the spanning tree T;  K_T = (-1)^{nH + C(n,2)C(H,2)} int rho^T ^ d^H y;
  - K -> (-1)^D K (edge reversal), sgn^{D+1} K (edge order), sgn^D K (vertex relabelling), D = H + T;
  - extraction sign formula;
  - series-parallel layer against the brute-force layer on every tight series-parallel graph
    of results/graph_counts.txt (PLAN.md E1), for T = 0, and on the two-loop theta graphs for T = 1;
  - independence of the reduction plan; local determinants = Jacobian of the polynomial currents;
  - n = 4, T = 1 (Pfaffian of a 4 x 4 matrix of two-forms).
"""
import itertools
import math
import os
import random
import re
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import htk
from htk import brute, catalog, oracle, sp_reduce, PiNumber, Graph
from htk import conventions as cv
from htk import spred
from htk.graph import polygon, theta, banana, fan, box, graph_A, graph_B, ladder
from htk.linalg import det, perm_sign, wedge_all

FAST = bool(os.environ.get("HTK_FAST"))
f = math.factorial
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")


def rsources(rng, G, H):
    return [[Fr(rng.randint(-5, 5)) for _ in range(H)] for _ in range(G.nv - 1)] if H else None


def rtimes(rng, G):
    return [Fr(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(G.ne)]


def all_trees(G):
    out = []
    B = G.incidence()
    for comb in itertools.combinations(range(G.ne), G.nv - 1):
        if det([B[e] for e in comb]) != 0:
            out.append(list(comb))
    return out


CASES = [(box(), 3, 0), (graph_A(), 3, 0), (graph_B(), 3, 0), (fan(1), 2, 0), (fan(2), 2, 0), (fan(3), 2, 0),
         (polygon(5), 4, 0), (theta(1, 2, 2)[0], 1, 1), (theta(1, 3, 3)[0], 2, 1), (theta(2, 2, 3)[0], 2, 1),
         (banana(3), 0, 1), (theta(1, 4, 4)[0], 3, 1), (theta(3, 3, 3)[0], 4, 0)]


def test_tree_independence():
    """sigma_T^{H+T} D_T is the same for all spanning trees; T = 0: D_T = sign_dy(n, H) D_dy."""
    rng = random.Random(201)
    for G, H, T in CASES:
        t, lam = rtimes(rng, G), rsources(rng, G, H)
        ref = brute.density(G, t, lam, T=T)
        trees = all_trees(G)
        if len(trees) > 40:
            trees = rng.sample(trees, 40)
        for tr in trees:
            DT = brute.density(G, t, lam, T=T, C=G.cycle_matrix(tr, positive=False))
            assert DT.scaled(G.sigma(tr) ** (H + T)).same_as(ref), (G.name, tr)
            if T == 0:
                assert DT.rational() == cv.sign_dy(G.loops, H) * brute.jacobian_density(G, t, lam, G.cotree(tr))
            # det[B^t ; C_T] = sigma_T * (number of spanning trees)
        B = G.incidence()
        W = [[B[e][v] for e in range(G.ne)] for v in range(G.nv - 1)] + G.cycle_matrix()
        assert det(W) == len(all_trees(G))
    print("  %d graphs" % len(CASES))


def test_exterior_algebra_vs_determinant():
    """T = 0: coefficient from the wedge product of the psi's (all e0) = determinant formula."""
    rng = random.Random(202)
    for G, H, T in CASES:
        if T:
            continue
        t, lam = rtimes(rng, G), rsources(rng, G, H)
        C = G.cycle_matrix()
        y = brute.currents(G, t, lam)
        assert brute.check_kirchhoff(G, t, lam, y)
        form = wedge_all(brute.psi_forms(G, t, y, C))
        U = brute.symanzik_U(G, t)
        from htk.linalg import basic_coefficient
        assert basic_coefficient(form, t) / U ** H == brute.density(G, t, lam, check_basic=True).rational()
    print("  ok")


def test_sign_rules():
    """Edge reversal: (-1)^D; edge order: sgn^{D+1}; vertex relabelling (root = vertex labelled last): sgn^D."""
    rng = random.Random(203)
    for G, H, T in CASES:
        D = H + T
        t, lam = rtimes(rng, G), rsources(rng, G, H)
        ref = brute.density(G, t, lam, T=T)
        if ref.coeff == 0:
            continue
        # edge reversal
        for k in (1, 2):
            which = rng.sample(range(G.ne), k)
            d2 = brute.density(G.reverse_edges(which), t, lam, T=T)
            assert d2.same_as(ref.scaled(cv.edge_reversal_sign(D, k))), (G.name, "reversal")
        # edge order
        for _ in range(3):
            order = list(range(G.ne))
            rng.shuffle(order)
            inv = [0] * G.ne
            for j, e in enumerate(order):
                inv[e] = j
            d2 = brute.density(G.permute_edges(order), [t[e] for e in order], lam, T=T)
            assert d2.same_as(ref.scaled(cv.edge_permutation_sign(D, perm_sign(inv)))), (G.name, "edge order")
        # vertex relabelling, including a change of root
        if H:
            full = lam + [[-sum(lam[v][i] for v in range(G.nv - 1)) for i in range(H)]]
        for _ in range(3):
            new = list(range(G.nv))
            rng.shuffle(new)
            G2 = G.relabel(new)
            lam2 = None
            if H:
                lam2 = [None] * G.nv
                for v in range(G.nv):
                    lam2[new[v]] = full[v]
            d2 = brute.density(G2, t, lam2, T=T)
            assert d2.same_as(ref.scaled(cv.vertex_relabelling_sign(D, perm_sign(new)))), (G.name, "relabel")
    print("  %d graphs" % len(CASES))


def test_extraction_sign_formula():
    """Position-space identity of CONVENTIONS.md Section 4 (S3 eq-extraction for D = 6)."""
    for G, H, T in CASES + [(polygon(7), 5, 1), (theta(1, 6, 6)[0], 5, 1), (theta(4, 4, 5)[0], 5, 1)]:
        D = H + T
        s = cv.extraction_sign(G, D)
        assert s == cv.extraction_sign_formula(G, D), G.name
        assert cv.extraction_sign(G, D, C=G.cycle_matrix(positive=False)) == s * G.sigma() ** D
    # S3 counting |I| = 6n+1, |V| = 5n+2: the formula gives (-1)^n
    for n in range(1, 9):
        N, V = 6 * n + 1, 5 * n + 2
        e = (N * (N - 1) // 2) * 15 + 6 * n + (V - 1) * n * 15
        assert (-1) ** e == (-1) ** n
        assert cv.varsigma_S3(n) == cv.sign_n(n)
    print("  ok")


def _count_line(d, n):
    path = os.path.join(ROOT, "results", "graph_counts.txt")
    if not os.path.exists(path):
        return None
    for line in open(path):
        m = re.match(r"d=(\d+) n=(\d+) \|V\|=\d+ \|I\|=\d+: total=(\d+) SP=(\d+) nonSP=\d+ CH=(\d+)", line)
        if m and int(m.group(1)) == d and int(m.group(2)) == n:
            return tuple(int(x) for x in m.groups()[2:])
    return None


def test_catalog_counts():
    """Enumeration, series-parallel test and chain reducibility against results/graph_counts.txt (S6)."""
    jobs = [(2, n) for n in range(1, 6)] + [(3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 1), (6, 1)]
    if not FAST:
        jobs.append((2, 6))
    for d, n in jobs:
        Gs = catalog.enumerate_tight(d, n)
        got = (len(Gs), sum(catalog.is_sp(G) for G in Gs), sum(catalog.chain_reducible(G, d) for G in Gs))
        want = _count_line(d, n)
        assert want is None or got == want, (d, n, got, want)
        if d == 2 and (n + 2) in catalog.LAMAN_COUNTS:
            assert got[0] == catalog.LAMAN_COUNTS[n + 2]
        print("  d=%d n=%d: total %d, SP %d, CH %d %s" % ((d, n) + got + ("" if want else "(no reference line)",)))
    assert [catalog.two_loop_thetas(d) for d in (2, 3)] == [[(1, 2, 2)], [(1, 3, 3), (2, 2, 3)]]
    assert catalog.two_loop_thetas(6) == [(1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5)]


def test_sp_all_tight_graphs_T0():
    """Every tight series-parallel graph of the table (PLAN.md E1): the planner finds tight chain
    moves, and the density (local determinants and generic Jacobian) equals the brute-force density
    at random rational points; the two integrals agree."""
    rng = random.Random(204)
    jobs = [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 1)]
    if not FAST:
        jobs.append((2, 5))
    count = 0
    for d, n in jobs:
        for G in catalog.tight_graphs(d, n):
            if not catalog.is_series_parallel(G):
                continue
            red = sp_reduce(G)
            assert red.tight_plan and red.D == d
            lam = rsources(rng, G, d)
            dl = red.density(lam, method="local", check=3, seed=count)
            Kl = dl.integrate()
            if G.ne <= (7 if FAST else 10):
                dg = red.density(lam, method="generic", check=1)
                assert dg.integrate() == Kl and dg.num == dl.num
            count += 1
    for pqr, G, interior in catalog.theta_graphs(5):
        red = sp_reduce(G)
        lam = rsources(rng, G, 5)
        Kl = red.kernel(lam, method="local", check=2)
        red3 = sp_reduce(G, prefer="large")              # one move with three chains
        assert len(red3.moves) == 1
        assert red3.kernel(lam, method="local", check=2) == Kl
        count += 1
    print("  %d graphs" % count)


def test_plan_independence():
    """Different reduction plans give the same kernel (T = 0 and T = 1)."""
    rng = random.Random(205)
    G = graph_A()
    lam = rsources(rng, G, 3)
    K1 = sp_reduce(G).kernel(lam)
    K2 = sp_reduce(G, prefer="large").kernel(lam)
    # the plan of S2: chain move (O,1,2,3) then the box through its vertices O, 3
    plan = [(5, 2, [[0, 1, 2], [3]]), (5, 2, [[7], [6, 5, 4]])]
    K3 = sp_reduce(G, plan=plan).kernel(lam)
    # generic plan without tightness
    red4 = spred.SPReduction(G, plan=spred.plan_generic(G))
    K4 = red4.kernel(lam, method="generic")
    assert K1 == K2 == K3 == K4
    for pqr, H in [((1, 2, 2), 1), ((2, 2, 3), 2), ((1, 3, 3), 2)]:
        G, _ = theta(*pqr)
        lam = rsources(rng, G, H)
        Ka = sp_reduce(G, prefer="large").kernel(lam, T=1)
        Kb = sp_reduce(G, prefer="small").kernel(lam, T=1)
        assert Ka == Kb, (pqr, Ka, Kb)
    print("  (A): K = %s for four plans; T = 1 thetas: one move = two moves" % K1)


def test_theta_T1_all_H():
    """Two-loop theta graphs at T = 1, H = 0..5 (D = 1..6):
    K = s (D_theta / 2) prod_P Gamma(H + 3/2 - n_P) / (n_P - 1)! / Gamma(H + 3/2), s = +-1 (printed);
    density = s (D_theta / 2) prod_P w_P^{H + 1/2 - n_P}."""
    rng = random.Random(206)
    table = []
    for D in range(1, 7):
        H = D - 1
        graphs = [((1, 1, 1),) + (banana(3), [[], [], []])] if D == 1 else catalog.theta_graphs(D)
        for pqr, G, interior in graphs:
            lam = rsources(rng, G, H)
            Dth = brute.theta_D(lam, interior) if H else Fr(1)
            red = sp_reduce(G, D=D, prefer="large")
            dens = red.density(lam, T=1, check=2)
            K = dens.integrate()
            c = sp.Rational(1, 2)
            for n_ in pqr:
                c *= sp.gamma(sp.Rational(2 * H + 3, 2) - n_) / f(n_ - 1)
            c /= sp.gamma(sp.Rational(2 * H + 3, 2))
            c = sp.nsimplify(sp.simplify(c / sp.pi))
            mag = PiNumber({1: Fr(int(c.p), int(c.q)) * Dth})
            s = 1 if K == mag else -1 if K == -mag else None
            assert s is not None, (pqr, H, K, mag)
            # sign observed for H = 0..5 with the labelling of graph.theta and the positive cycle basis
            assert s == (-1) ** ((H + 1) // 2 + (H + 1) * (pqr[1] + 1)), (pqr, H, s)
            assert dens.num.is_ground
            expo = dens.exponents()
            assert all(expo.get("w1_%d" % (i + 1), Fr(0)) == Fr(2 * H + 1, 2) - pqr[i] for i in range(3))
            table.append((H, pqr, s, c))
    print("  (H, pqr, s, c/pi):", "; ".join("%d %s %+d %s" % row for row in table))


def test_T1_four_loops():
    """n = 4, T = 1, H = 1 (|I| = 9, |V| = 6): Pfaffian of a 4 x 4 matrix of two-forms in the
    brute-force layer and in the ratio coordinates; all tight series-parallel graphs of (D, n) = (2, 4)."""
    rng = random.Random(207)
    vals = []
    for G in catalog.tight_graphs(2, 4):
        if not catalog.is_series_parallel(G):
            continue
        lam = rsources(rng, G, 1)
        red = sp_reduce(G, prefer="large")
        dens = red.density(lam, T=1, check=2)
        vals.append(dens.integrate())
        # oracle integrand (general-n branch) against the brute-force layer, pointwise
        import numpy as np
        t = rtimes(rng, G)
        ex = float(brute.density(G, t, lam, T=1))
        nu = oracle.density_batch(G, np.array([[float(x) for x in t]]), np.array(lam, dtype=float), 1)[0]
        assert abs(ex - nu) <= 1e-9 * max(abs(ex), 1e-300), (ex, nu)
    assert len(vals) == 5
    print("  K =", ", ".join(str(v) for v in vals))
    # n = 4, H = 2: the four-box ladder on R x C^2 (|I| = 13)
    G, cot, a, b = ladder(4)
    rng2 = random.Random(1)
    lam = rsources(rng2, G, 2)
    red = sp_reduce(G, prefer="small")
    dens = red.density(lam, T=1, check=2)
    K = dens.integrate()
    assert K == PiNumber({2: Fr(29127, 50)}), K
    if not FAST:
        assert sp_reduce(G, prefer="large").kernel(lam, T=1, check=1) == K
    print("  ladder(4), H = 2, T = 1: K = %s; exponents %s" % (K, {k: str(v) for k, v in dens.exponents().items()}))


def test_moments_and_symbolic():
    """Insertions: exact moments against the brute-force integrand and the oracle; symbolic sources at T = 1."""
    rng = random.Random(208)
    G = fan(2)
    lam = rsources(rng, G, 2)
    red = sp_reduce(G)
    powers = {(0, 0): 1, (3, 1): 2}
    exact = red.integrate(lam, powers=powers)
    val, err = oracle.kernel(G, lam, powers=powers)
    assert abs(val - float(exact)) < 1e-5 * max(1.0, abs(float(exact))), (val, exact)
    # pointwise: insertion = monomial of the Kirchhoff currents
    pt = red.random_point(rng)
    t = red.tmap(pt)
    y = brute.currents(G, t, lam)
    assert y == red.currents_at(lam, pt)
    z = [[Fr(rng.randint(-3, 3)) for _ in range(2)] for _ in range(G.ne)]
    X = red.displacement()
    vals = red._vals(pt)
    lhs = brute.exponent(y, z)
    rhs = sum(lam[v][i] * spred._evalp(X[v][e], vals) * z[e][i] for v in range(G.nv - 1) for e in range(G.ne) for i in range(2))
    assert lhs == rhs
    # symbolic sources, T = 1
    G, interior = theta(1, 2, 2)
    a, b, c = sp.symbols("a b c")
    red = sp_reduce(G, prefer="large")
    K = red.kernel([[a], [b], [c]], T=1)
    assert sp.expand(K - (-sp.pi / 3) * b * c) == 0, K
    print("  fan2 moment: exact %s, oracle %.9g; theta(1,2,2), H = 1: K = %s" % (exact, val, K))


def test_oracle_T1_slow():
    """|I| = 7, T = 1: theta(1,3,3) on R x C^2 (about 2 minutes); skipped with HTK_FAST."""
    if FAST:
        print("  skipped (HTK_FAST)")
        return
    rng = random.Random(209)
    G, interior = theta(1, 3, 3)
    while True:
        lam = rsources(rng, G, 2)
        exact = sp_reduce(G, prefer="large").kernel(lam, T=1)
        if float(exact) != 0:
            break
    val, err = oracle.kernel(G, lam, T=1)
    rel = abs(val - float(exact)) / abs(float(exact))
    assert rel < 2e-3, (val, float(exact))
    print("  theta(1,3,3), H = 2: exact %s = %.8f, oracle %.8f, rel.err %.1e (est. %.1e)" % (
        exact, float(exact), val, rel, err / abs(float(exact))))


TESTS = [test_tree_independence, test_exterior_algebra_vs_determinant, test_sign_rules,
         test_extraction_sign_formula, test_catalog_counts, test_sp_all_tight_graphs_T0,
         test_plan_independence, test_theta_T1_all_H, test_T1_four_loops, test_moments_and_symbolic,
         test_oracle_T1_slow]


def main(argv):
    names = set(argv[1:])
    failed, passed = [], []
    for t in TESTS:
        if names and t.__name__ not in names:
            continue
        t0 = time.time()
        print("%s" % t.__name__, flush=True)
        try:
            t()
            passed.append(t.__name__)
            print("  PASS (%.1fs)" % (time.time() - t0), flush=True)
        except Exception as exc:                                  # noqa: BLE001
            import traceback
            traceback.print_exc()
            failed.append(t.__name__)
            print("  FAIL: %r" % (exc,), flush=True)
    print("passed %d, failed %d" % (len(passed), len(failed)))
    if failed:
        print("failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
