"""Exact checks of htk.star (star-triangle moves, T in {0, 1}).

Run:  python3 code/tests/test_star.py        (HTK_FAST=1: fewer graphs and points)

Statements checked (explore/circuit-general/ht-d2/REPORT.md, Section 1), at random rational points:
  - forward(t(f, T)) = (f, T); currents of G = reverse current map of the reduced currents
    (source-free and sourced centres); U_Gamma = U_red * prod t / prod T * prod_v G_v;
  - J_aff = closed local Jacobian (single move), |J_aff| and U ratio in closed form (paired move);
  - D_Gamma(t) J_aff = Kron decomposition from reduced data (T = 1, H = 1; T = 0, H = 2), single and
    paired moves, sourced centres;
  - closed transformation formulas of psi (all components, H = 1, 2, 3), of the topological Pfaffian
    and of the period matrix (det M_gg, Schur complements), single move;
  - LLC local factors at T = 0, H = 2 (single: r1 r2 <y_ux, y_wx>; paired:
    -(r1 - r2)(1 - r3 - r4) <y_uA, y_AB><y_uB, y_AB>);
  - xi^0 identities at T = 1, H = 1 (single: (-1)^(n-1) (r1 r2 t_vx)^(1/2) D'_red[<y_ux, y_wx>_{tau,lambda}];
    paired: -(U/U_red)^(1/2)(r1 - r2)(1 - r3 - r4) D'_red[<y_uA, y_AB><y_uB, y_AB>_{tau,lambda}]);
  - T = 1 densities of odd-loop graphs vanish identically (n = 3 forests).
"""
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "explore", "circuit-general", "ht-d2"))

from htk import graph, brute, star  # noqa: E402
from forests import forests, forest_graph, spider_graph  # noqa: E402

FAST = bool(os.environ.get("HTK_FAST"))


def rl(rng, G, H, zero=()):
    lam = [[Fr(rng.randint(-5, 5)) for _ in range(H)] for _ in range(G.nv - 1)]
    for v in zero:
        lam[v] = [Fr(0)] * H
    return lam


def graphs():
    out = [graph.fan(2), graph.fan(4)]
    for n in ((3, 4) if FAST else (3, 4, 5)):
        out += [forest_graph(F)[0] for F in forests(n)]
    out += [spider_graph(a)[0] for a in ((1, 1, 2), (1, 2, 1, 2))]
    return out


def test_single():
    rng = random.Random(1)
    count = 0
    for G in graphs():
        for (v, u, w, x) in star.find_single(G)[:1 if FAST else 2]:
            mv = star.SingleMove(G, v)
            f, T = mv.random_point(rng)
            assert mv.jacobian(f, T) == mv.J_closed(f, T)
            mv.selfcheck(rl(rng, G, 1), 1, 1, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("kron", "closed"))
            mv.selfcheck(rl(rng, G, 2), 2, 0, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("kron", "closed"))
            mv.selfcheck(rl(rng, G, 2, (v,)), 2, 0, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("llc_T0",))
            for H in (1, 2, 3):
                f, T = mv.random_point(rng)
                assert mv.check_components(f, T, rl(rng, G, H))
            f, T = mv.random_point(rng)
            lhs, rhs = mv.xi0_T1_identity(f, T, rl(rng, G, 1, (v,)))
            assert lhs.same_as(rhs)
            if G.loops % 2:
                t = mv.tmap(f, T)
                assert brute.density(G, t, rl(rng, G, 1), T=1).coeff == 0
            count += 1
    print("  single moves checked on %d (graph, centre) pairs" % count)


def test_paired():
    rng = random.Random(2)
    count = 0
    for G in graphs():
        for (v1, v2, u, A, B) in star.find_paired(G)[:1]:
            mv = star.PairedMove(G, v1, v2)
            f, T = mv.random_point(rng)
            assert abs(mv.jacobian(f, T)) == mv.J_abs_closed(f, T)
            assert mv.U_ratio(f, T) == mv.U_ratio_closed(f, T)
            mv.selfcheck(rl(rng, G, 1), 1, 1, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("kron",))
            mv.selfcheck(rl(rng, G, 2), 2, 0, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("kron",))
            mv.selfcheck(rl(rng, G, 2, (v1, v2)), 2, 0, npoints=1, seed=rng.randint(0, 10 ** 6), methods=("llc_T0",))
            f, T = mv.random_point(rng)
            lhs, rhs = mv.xi0_T1_identity(f, T, rl(rng, G, 1, (v1, v2)))
            assert lhs.same_as(rhs)
            count += 1
    print("  paired moves checked on %d graphs" % count)


if __name__ == "__main__":
    for name in ("test_single", "test_paired"):
        print(name)
        globals()[name]()
    print("all passed")
