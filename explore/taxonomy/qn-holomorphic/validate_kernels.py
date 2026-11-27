"""Kernel validation of kern_h.GraphKernel (general H).  Log: validate_kernels.log
(a) H=2 fans: I = (-1)^{n(n+1)/2} K = det(lam_end1, lam_end2)^n / ((n!)^3 (n+1))   [src LLC eq-1-42; htk CONVENTIONS S1]
(b) H=3 rooted trees (labelled as rooted_trees.graphs): I(probe) = 6d-nilpotency REPORT table (n=4) and spred
(c) H=4 n=3: the 12 cone kernels of circuit-general/atlas/e1_cone_kernels.log (same graphs, same random sources)
(d) relabelling sign rule  K(pi Gamma) = sgn(pi_v)^H (-1)^{H k_rev} sgn(pi_e)^{H+1} K(Gamma)  on the rgraphs labelling
(e) H=1 banana n=1, H=5 trees n=1,2: values at random sources (new)
"""
import os
import sys
import random
import time
from fractions import Fraction as Fr
from itertools import cycle
from math import factorial

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'circuit-general', 'atlas'))

import kern_h
from kern_h import GraphKernel, modq, P
import rgraphs
from rooted_trees import graphs as graphs6
import htk
from htk.graph import from_networkx


def exact_I(gk, lam):
    """I = sign_n K by spred (Fractions), lam: v -> H-tuple."""
    K = gk.red.kernel({v: list(lam[v]) for v in lam}, check=0) if gk.red.tight_plan else \
        htk.kernel_exact(gk.G, {v: list(lam[v]) for v in lam}, check=0)
    return gk.sign_n * K


def numeric_I(gk, lam_list):
    L = len(lam_list)
    lam = {v: np.array([[lam_list[r][v][a] % P for r in range(L)] for a in range(gk.H)], dtype=np.int64)
           for v in gk.outputs}
    return gk.kernel_numeric(lam, L)


def part_a():
    print('(a) H=2 fans: I (rgraphs labelling) vs s_n det(lam_end1, lam_end2)^n/((n!)^3(n+1)), residues mod P', flush=True)
    for n in range(1, 7):
        g = rgraphs.rooted_trees(2, n)[0]
        gk = GraphKernel(2, g['nb'], g['edges'])
        lexp, coef = gk.kernel_symbolic(((),) * len(gk.spokes))
        u, w = gk.outputs
        syms = {nm: sp.Symbol(nm) for nm in gk.lam_names}
        det = syms['l%d_0' % u] * syms['l%d_1' % w] - syms['l%d_1' % u] * syms['l%d_0' % w]
        ref = sp.Poly(sp.expand(det ** n), *[syms[x] for x in gk.lam_names])
        den = factorial(n) ** 3 * (n + 1)
        refd = {m: modq(Fr(int(c), den)) for m, c in zip(ref.monoms(), ref.coeffs())}
        mine = {tuple(int(x) for x in lexp[t]): int(coef[t]) % P for t in range(lexp.shape[0]) if coef[t] % P}
        s_n = (-1) ** (n * (n - 1) // 2)
        ok = mine == {m: (s_n * v) % P for m, v in refd.items()}
        print('  n=%d ends %s: I = (-1)^{n(n-1)/2} det^n/((n!)^3(n+1)) [%d monomials]: %s' % (n, (u, w), len(mine), ok),
              flush=True)


def part_b(nmax=4):
    print('(b) H=3 rooted trees vs spred, REPORT table at n=4', flush=True)
    report4 = [Fr(-683, 447897600), Fr(7, 3359232), Fr(427, 44789760), Fr(-1441, 134369280),
               Fr(-3139, 2239488000), Fr(7, 3359232), Fr(-7, 3359232), Fr(-2743, 223948800),
               Fr(5843, 3359232000), Fr(12541, 6718464000), Fr(-37, 19906560), Fr(787, 895795200),
               Fr(1613, 4478976000), Fr(-77, 44789760), Fr(77, 44789760)]
    rng = random.Random(7)
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for n in range(1, nmax + 1):
        for k, edges in enumerate(graphs6(n), 1):
            t0 = time.time()
            nb = 2 * n + 1
            gk = GraphKernel(3, nb, edges)
            probe = dict(zip(gk.outputs, cycle(basis)))
            rnd = {v: tuple(rng.randint(-3, 3) for _ in range(3)) for v in gk.outputs}
            val = numeric_I(gk, [probe, rnd])
            ref0 = exact_I(gk, probe) if (n < 4 or gk.red.tight_plan or True) else None
            ok0 = int(val[0]) == modq(ref0)
            okr = None
            if gk.red.tight_plan:
                okr = int(val[1]) == modq(exact_I(gk, rnd))
            rep = (report4[k - 1] == ref0) if n == 4 else '-'
            print('  G%d.%d plan %s: I(probe) = %s  XP==spred %s, random %s, =REPORT %s  [%.1fs]'
                  % (n, k, 'tight' if gk.red.tight_plan else 'non-tight', ref0, ok0, okr, rep, time.time() - t0),
                  flush=True)


def part_c():
    print('(c) H=4 n=3 cone kernels vs e1_cone_kernels.log', flush=True)
    from atlas import set_cone
    logv = [Fr(34835017, 41472000), Fr(8447, 2073600), Fr(-18845633, 3888000), Fr(862993709, 18662400),
            Fr(2201, 6912), Fr(-143507, 248832), Fr(7143371, 37324800), Fr(47569949, 7464960),
            Fr(-2937451, 4976640), Fr(-146710723, 7776000), Fr(-11557, 259200), Fr(7914137, 11197440)]
    rng = random.Random(20261123)
    items = set_cone(4, 3)
    mine = rgraphs.rooted_trees(4, 3)
    print('  set_cone(4,3): %d graphs; rgraphs.rooted_trees(4,3): %d' % (len(items), len(mine)), flush=True)
    for idx, (G, O, S0) in enumerate(items):
        Gh, lab = from_networkx(G, root=O)
        lam = {lab[v]: [rng.randint(-3, 3) for _ in range(4)] for v in S0 if v != O}
        t0 = time.time()
        gk = GraphKernel(4, Gh.nv - 1, Gh.edges) if all(e[1] == Gh.nv - 1 or Gh.nv - 1 not in e for e in Gh.edges) \
            else None
        if gk is None:
            # orient spokes towards the root (edge reversal sign (-1)^{H k} = +1 at H = 4)
            edges = [(a, b) if b == Gh.nv - 1 or a != Gh.nv - 1 else (b, a) for a, b in Gh.edges]
            gk = GraphKernel(4, Gh.nv - 1, edges)
        lam_o = {v: tuple(lam.get(v, [0] * 4)) for v in gk.outputs}
        assert all(v in gk.outputs for v in lam), 'sources outside output vertices'
        val = numeric_I(gk, [lam_o])
        K = Fr(int(kern_h.ratrec(int(val[0]) * gk.sign_n % P)))
        print('  #%d K = %s  log %s  equal %s  [%.1fs]' % (idx, K, logv[idx], K == logv[idx], time.time() - t0),
              flush=True)


def part_d():
    print('(d) relabelling sign rule on rgraphs labelling (H=3 n<=3 against rooted_trees labelling; H=2,4 random)',
          flush=True)
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    rng = random.Random(3)

    def perm_sign(p):
        p = list(p)
        s = 1
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]:
                    s = -s
        return s
    for H, n in ((2, 2), (2, 3), (3, 2), (3, 3), (4, 2)):
        for g in rgraphs.rooted_trees(H, n):
            nb, edges = g['nb'], g['edges']
            gk = GraphKernel(H, nb, edges)
            # random relabelling of bulk vertices, edge order, reversal of non-spoke edges
            perm = list(range(nb))
            rng.shuffle(perm)
            perm.append(nb)
            order = list(range(len(edges)))
            rng.shuffle(order)
            new = []
            krev = 0
            for e in order:
                a, b = perm[edges[e][0]], perm[edges[e][1]]
                if b != nb and rng.randrange(2):
                    a, b = b, a
                    krev += 1
                new.append((a, b))
            gk2 = GraphKernel(H, nb, new)
            lam = {v: tuple(rng.randint(-3, 3) for _ in range(H)) for v in gk.outputs}
            lam2 = {perm[v]: lam[v] for v in lam}
            I1, I2 = exact_I(gk, lam), exact_I(gk2, lam2)
            pred = perm_sign(perm[:nb]) ** H * (-1) ** (H * krev) * perm_sign(order) ** (H + 1)
            print('  H=%d n=%d %s: K2/K1 = %s predicted %s  %s' % (H, n, edges, I2 / I1 if I1 else None, pred,
                                                                   I2 == pred * I1), flush=True)


def part_e():
    print('(e) H=1 banana, H=5 trees (random integer sources)', flush=True)
    rng = random.Random(11)
    for H, n in ((1, 1), (5, 1), (5, 2)):
        for g in rgraphs.rooted_trees(H, n):
            t0 = time.time()
            gk = GraphKernel(H, g['nb'], g['edges'])
            lam = {v: tuple(rng.randint(-3, 3) for _ in range(H)) for v in gk.outputs}
            I = exact_I(gk, lam)
            val = numeric_I(gk, [lam])
            print('  H=%d n=%d %s |Aut|=%d plan %s: I = %s  XP==spred %s  [%.1fs]' % (
                H, n, g['edges'], g['aut'], 'tight' if gk.red.tight_plan else 'non-tight', I,
                int(val[0]) == modq(I), time.time() - t0), flush=True)


if __name__ == '__main__':
    parts = sys.argv[1:] or ['a', 'b', 'c', 'd', 'e']
    for p in parts:
        {'a': part_a, 'b': part_b, 'c': part_c, 'd': part_d, 'e': part_e}[p]()
