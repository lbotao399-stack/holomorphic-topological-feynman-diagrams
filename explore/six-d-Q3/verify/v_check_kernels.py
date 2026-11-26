"""C3, C4: exact pointwise comparison of the explorer's densities and displacement maps with the
brute-force density of S2 eq-kernel-def (own nodal Kirchhoff solver, own chart t(r), vk.py)."""
import random
import sys
import time
from fractions import Fraction as Fr

import vk
from vk import GRAPHS, Network, density, reduction_factors, J_formula, X_formula, exponent, det3, random_point

NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
rng = random.Random(20261120)

MOVE_COTREES = {'G1': [0, 4, 7], 'G2': [0, 6, 4], 'G3': [0, 4, 7], 'G4': [0, 6, 4], 'G5': [0, 4, 7]}
EXPL_EPS = {'G1': -1, 'G2': -1, 'G3': 1, 'G4': 1, 'G5': -1}


def rvec(lo=-6, hi=6):
    return tuple(Fr(rng.randint(lo, hi)) for _ in range(3))


def s2_J(name, r, l):
    if name == 'box':
        return det3(l[1], l[2], l[3])
    r1, r2, r3, r4, r5, r6 = r
    if name == 'A':
        return (r1 ** 2 * (r5 * det3(l[1], l[2], l[4]) + r6 * det3(l[1], l[2], l[5]))
                * ((1 - r2) * det3(l[1], l[4], l[5]) + (1 - r3) * det3(l[2], l[4], l[5])))
    return (r1 * (1 - r1) * (r5 * det3(l[1], l[2], l[4]) + r6 * det3(l[1], l[2], l[5]))
            * ((1 - r2) * det3(l[1], l[4], l[5]) + (1 - r3) * det3(l[2], l[4], l[5])))


def s2_X(name, r, z):
    if name == 'box':
        z1, z2 = z
        return {j: vk.vadd((r[j - 1], z1), (1 - r[j - 1], z2)) for j in (1, 2, 3)}
    r1, r2, r3, r4, r5, r6 = r
    z1, z2, z3 = z
    W = vk.vadd((1 - r1, z1), (r1, z2), (-1, z3))
    Y = vk.vadd((1, z3), (r4, W))
    X = {1: vk.vadd((r2, z1), (1 - r2, Y)), 4: vk.vadd((1, z3), (r5, W)), 5: vk.vadd((1, z3), (r6, W))}
    X[2] = vk.vadd((r3, z1), (1 - r3, Y)) if name == 'A' else vk.vadd((r3, z2), (1 - r3, Y))
    return X


S2_SIGN = {'box': 1, 'A': 1, 'B': -1}

t00 = time.time()
for name in ['box', 'A', 'B', 'G1', 'G2', 'G3', 'G4', 'G5']:
    t0 = time.time()
    g = GRAPHS[name]
    net = Network(g['edges'], g['verts'])
    trees = []
    while len(trees) < 4:
        T = net.random_tree(rng)
        if T not in trees:
            trees.append(T)
    if name in MOVE_COTREES:
        mt = sorted(set(range(len(g['edges']))) - set(MOVE_COTREES[name]))
        trees = [mt] + trees
    elif name == 'box':
        trees = [[1, 2, 3]] + trees
    outs = g['outputs']
    ok_J = ok_tree = ok_red = ok_ratio = ok_cur = ok_X = True
    nres = 0
    for ip in range(NPTS):
        r = random_point(name, rng)
        # (a) sources at the outputs only: explorer's compact J (or S2's J) with the I-sign
        lam = {v: rvec() for v in outs}
        dens = [density(name, r, lam, T)[0] for T in trees]
        if len(set(dens)) != 1:
            ok_tree = False
        ref = (vk.I_SIGN[name] * J_formula(name, r, lam)) if name.startswith('G') else S2_SIGN[name] * s2_J(name, r, lam)
        if dens[0] != ref:
            ok_J = False
            print('   J mismatch', name, ip, dens[0], ref)
        if dens[0] != 0:
            nres += 1
        # displacement map
        nz = 4 if name.startswith('G') else len(g['spokes'])
        z = [rvec(-9, 9) for _ in range(nz)]
        _, y, _ = density(name, r, lam, trees[0])
        lhs = exponent(name, y, z)
        X = X_formula(name, r, z) if name.startswith('G') else s2_X(name, r, z)
        rhs = sum(sum(lam[v][a] * X[v][a] for a in range(3)) for v in outs)
        if lhs != rhs:
            ok_X = False
            print('   X mismatch', name, ip, lhs, rhs)
        # (b) sources at all vertices: s * prod(local factors) through reduced networks
        lam2 = {v: rvec() for v in g['verts']}
        d2 = [density(name, r, lam2, T)[0] for T in trees[:3]]
        if len(set(d2)) != 1:
            ok_tree = False
        rr, prod, okc = reduction_factors(name, r, lam2)
        if rr != r:
            ok_ratio = False
        if not okc:
            ok_cur = False
        sgn = vk.S_SIGN.get(name, {'box': 1, 'A': -1, 'B': 1}.get(name))
        if d2[0] != sgn * prod:
            ok_red = False
            print('   product mismatch', name, ip, d2[0], prod)
    sT = {}
    for T in trees[:2]:
        sT[tuple(sorted(set(range(len(g['edges']))) - set(T)))] = net.tree_sign(T)[0]
    print('%s: %d points (%d nonzero), trees (cotrees) %s; density tree-independent %s; density == I-sign*J %s; '
          'X map %s; general sources: chart ratios recovered %s, reduced currents == full %s, density == s*prod(local) %s; '
          'eps_T det B_T %s  [%.1fs]'
          % (name, NPTS, nres, [sorted(set(range(len(g['edges']))) - set(T)) for T in trees], ok_tree, ok_J, ok_X,
             ok_ratio, ok_cur, ok_red, sT, time.time() - t0), flush=True)
    if name in EXPL_EPS:
        mt = trees[0]
        alt = vk.perm_sign(mt + MOVE_COTREES[name]) * vk.det([[net.B[e][i] for e in mt] for i in range(len(g['verts']))])
        print('   move cotree %s: eps_T det B_T = %d with cotree in increasing order (S2 eq-kernel-def); '
              '%d with cotree in the order %s; explorer %d'
              % (MOVE_COTREES[name], net.tree_sign(mt)[0], alt, MOVE_COTREES[name], EXPL_EPS[name]))
    if name in ('A', 'B'):
        mt = sorted(set(range(len(g['edges']))) - {0, 4})
        print('   cotree [0, 4]: eps_T det B_T = %d' % net.tree_sign(mt)[0])
print('total %.1fs' % (time.time() - t00))
