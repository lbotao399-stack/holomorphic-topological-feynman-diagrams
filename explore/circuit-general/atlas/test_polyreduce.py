"""Exact checks of polyreduce against htk.brute and known kernels.

(1) series-parallel graphs: polynomial density = spred density value (fan(3) at D=2, theta(1,3,3) at D=3);
(2) C^2 three-loop K4 subdivision (nonSP/README.md): kernel polynomial with a source-free star;
(3) C^3 three-loop K4[1,1,1,2,2,3]a, four marks (nonSP/VERIFY.md section 3.2): -1/77760;
(4) pointwise identity g(xi) = D(t(xi)) |det[t; dt/dxi]| and Kirchhoff currents at random rational points,
    including a sourced star-triangle (rational transfer).
"""
import random
import sys
from fractions import Fraction as Fr

import networkx as nx

from polyreduce import FullReduction, NotPolynomial
from reduce_search import Network, search
from htk.graph import Graph, fan, theta
from htk import spred


def plan_from_search(G, S0, allow):
    net = Network(G.to_networkx() if hasattr(G, "to_networkx") else G, S0)
    res = search(net, allow)
    if not res["reducible"]:
        return None, res
    return [("S" if k == "S" else "Y", v) for v, k in res["sequence"]], res


def nxg(G):
    H = nx.Graph()
    H.add_nodes_from(range(G.nv))
    H.add_edges_from(G.edges)
    return H


def run_full(G, lam, S0, allow=("S", "Y0"), npts=3, label=""):
    plan, res = plan_from_search(nxg(G), S0, allow)
    if plan is None:
        print(label, "not reducible with", allow)
        return None
    red = FullReduction(G, plan)
    rng = random.Random(11)
    for _ in range(npts):
        xi = red.random_point(rng)
        g, b, Jt = red.check_point(lam, xi)
        assert g == b, (g, b)
    out = None
    try:
        dens = red.poly_density(lam)
        out = red.integrate(dens)
        # pointwise agreement of the polynomial density with the AD value
        for _ in range(npts):
            xi = red.random_point(rng)
            g, b, _ = red.check_point(lam, xi)
            vals = {red.el.coords[j]: xi[red.el.coords[j]] for j in range(red.el.k)}
            pv = dens(*[vals[c] for c in red.el.coords])
            assert Fr(int(pv.numerator), int(pv.denominator)) == g
        print(label, "plan", plan, "| pointwise ok (%d pts) | K =" % (2 * npts), out)
    except NotPolynomial as e:
        print(label, "plan", plan, "| pointwise ok (%d pts) | rational transfer:" % npts, e)
    return out


if __name__ == "__main__":
    print("(1) series-parallel")
    G = fan(3)
    lam = [[1, 2], [0, 0], [0, 0], [3, -1]]
    K = run_full(G, lam, [0, 3, 4], label="fan3 D=2")
    Ks = spred.sp_reduce(G).kernel(lam)
    print("   spred:", Ks, "equal:", K == Ks)
    G, interior = theta(1, 3, 3)
    lam = [[0, 0, 0]] * (G.nv - 1)
    lam = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [2, -1, 3], [5, 1, -2]][:G.nv - 1]
    K = run_full(G, lam, list(range(G.nv)), label="theta133 D=3")
    Ks = spred.sp_reduce(G).kernel(lam)
    print("   spred:", Ks, "equal:", K == Ks)

    print("(2) C^2 K4 subdivision: edges (01,02,03,14,42,13,23), root 4")
    G = Graph(5, [(0, 1), (0, 2), (0, 3), (1, 4), (4, 2), (1, 3), (2, 3)], "K4sub")
    def det2(a, b):
        return a[0] * b[1] - a[1] * b[0]
    for lam in ([[2, 1], [1, 3], [-1, 2], [0, 0]], [[1, -2], [3, 1], [2, 5], [0, 0]]):
        # lambda_3 = 0: star at 3 source-free
        l = [list(map(Fr, v)) for v in lam]
        l4 = [-sum(v[i] for v in l) for i in range(2)]
        A = [l[1][i] - l[2][i] for i in range(2)]
        B = [l[0][i] - l[3][i] for i in range(2)]
        C = [l[0][i] + l[3][i] for i in range(2)]
        Dv = l4
        ref = Fr(det2(A, Dv), 2304) * (-3 * det2(A, B) ** 2 + 6 * det2(B, C) ** 2 - det2(B, Dv) ** 2)
        K = run_full(G, lam, [0, 1, 2, 4], label="K4sub lam3=0")
        print("   nonSP/README formula:", ref, "equal:", K == ref)
    lam = [[2, 1], [1, 3], [-1, 2], [1, 1]]
    run_full(G, lam, [0, 1, 2, 3, 4], allow=("S", "Y0", "Y1"), label="K4sub all sourced")

    print("(3) C^3 K4[1,1,1,2,2,3]a, four marks lambda_0,3,6 = e1,e2,e3, root 7")
    E = [(0, 4), (0, 7), (1, 2), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4), (5, 6), (5, 7)]
    G = Graph(8, E, "K4_111223a")
    lam = {0: [1, 0, 0], 3: [0, 1, 0], 6: [0, 0, 1]}
    lamlist = [lam.get(v, [0, 0, 0]) for v in range(7)]
    K = run_full(G, lamlist, [0, 3, 6, 7], label="K4[1,1,1,2,2,3]a")
    print("   nonSP/VERIFY 3.2: -1/77760; equal:", K == Fr(-1, 77760))
    print("(4) the other five C^3 K4 subdivisions, four marks: polynomial reducibility")
    others = {
        "K4[1,1,1,1,3,3]": [(0, 4), (0, 7), (1, 2), (1, 5), (1, 7), (2, 3), (2, 4), (3, 6), (4, 5), (5, 6)],
        "K4[1,1,1,2,2,3]b": [(0, 4), (0, 7), (1, 2), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4), (5, 6), (6, 7)],
        "K4[1,1,1,2,2,3]c": [(0, 4), (0, 7), (1, 2), (1, 5), (2, 3), (2, 4), (3, 6), (3, 7), (4, 5), (5, 6)],
        "K4[1,1,2,2,2,2]a": [(0, 1), (0, 4), (0, 7), (1, 2), (2, 3), (2, 5), (3, 6), (3, 7), (4, 5), (5, 6)],
        "K4[1,1,2,2,2,2]b": [(0, 1), (0, 5), (0, 6), (1, 2), (1, 7), (2, 3), (3, 4), (3, 6), (4, 5), (4, 7)],
    }
    for name, E in others.items():
        G = Graph(8, E, name)
        biv = [v for v in range(8) if G.degree(v) == 2]
        lamlist = [[0, 0, 0] for _ in range(7)]
        for v, e in zip([b for b in biv if b != 7], ([1, 0, 0], [0, 1, 0], [0, 0, 1])):
            lamlist[v] = e
        plan, res = plan_from_search(nxg(G), biv, ("S", "Y0"))
        print("  ", name, "marks", biv, "polynomial:", plan is not None)
        if plan is None:
            run_full(G, lamlist, biv, allow=("S", "Y0", "Y1"), label="   " + name + " (with Y1)")
