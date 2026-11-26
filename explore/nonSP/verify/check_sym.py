"""Exact pointwise automorphism identities with the fresh solver (vk.density_exact).
For a vertex permutation pi that is a graph automorphism, (P t)_{pi(e)} = t_e.  Lebesgue measure on the simplex is
P-invariant, so D(P t, lam') = sigma D(t, lam) for all t implies K(lam') = sigma K(lam)."""
import sys, random
from fractions import Fraction as Fr
from vk import G, density_exact
sys.path.insert(0, '..')
from graphs import K4SUB, PRISM, K33
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
rng = random.Random(7)


def Pt(E, pi, t):
    idx = {frozenset(e): k for k, e in enumerate(E)}
    out = [None] * len(E)
    for k, (a, b) in enumerate(E):
        out[idx[frozenset((pi.get(a, a), pi.get(b, b)))]] = t[k]
    return out


def ratio(nv, E, pi, lam, lam2, npts=4):
    g = G(nv, E)
    rs = []
    for _ in range(npts):
        t = [Fr(rng.randint(1, 50), rng.randint(1, 50)) for _ in E]
        a = density_exact(g, t, lam)
        b = density_exact(g, Pt(E, pi, t), lam2)
        rs.append(b / a if a != 0 else ('D=0', b))
    return rs


def lamlist(nv, d, H):
    return [[Fr(x) for x in d.get(v, [0] * H)] for v in range(nv - 1)]


if __name__ == "__main__":
    # a7: K4[1,1,2,2,2,2]a, pi = (3 5)(4 7), sources e1,e2,e3,-(sum) at 1,4,6,7 (same sources on both sides)
    name, E = K4SUB[4]
    lam = lamlist(8, {1: [1, 0, 0], 4: [0, 1, 0], 6: [0, 0, 1]}, 3)
    print("a7", name, "pi=(3 5)(4 7): D(Pt,lam)/D(t,lam) =", [str(x) for x in ratio(8, E, {3: 5, 5: 3, 4: 7, 7: 4}, lam, lam)])
    # same, with generic sources at the degree-two vertices: D(Pt, lam o pi^-1)/D(t, lam)
    l = {1: [2, -1, 3], 4: [1, 4, -2], 6: [-3, 1, 1]}; l[7] = [-sum(l[v][i] for v in l) for i in range(3)]
    l2 = {1: l[1], 6: l[6], 5: [0, 0, 0], 3: [0, 0, 0], 7: l[4], 4: l[7]}
    print("a7 generic sources, lam' = lam o pi^-1: D(Pt,lam')/D(t,lam) =", [str(x) for x in ratio(8, E, {3: 5, 5: 3, 4: 7, 7: 4}, lamlist(8, l, 3), lamlist(8, l2, 3))])
    # b3: K33, sources at 0,1,3,5 generic, pi = (2 4)
    l = {0: [1, 0], 1: [0, 1], 3: [Fr(1, 3), Fr(2, 5)]}
    L = lamlist(6, l, 2)
    print("b3 K33 src {0,1,3,5}, pi=(2 4):", [str(x) for x in ratio(6, K33, {2: 4, 4: 2}, L, L)])
    l = {0: [3, -1], 1: [Fr(2, 7), 5], 3: [-2, Fr(1, 3)]}
    L = lamlist(6, l, 2)
    print("b3 K33 src {0,1,3,5} (other lam), pi=(2 4):", [str(x) for x in ratio(6, K33, {2: 4, 4: 2}, L, L)])
    l = {0: [1, 0], 1: [0, 1], 2: [-1, -1]}
    L = lamlist(6, l, 2)
    print("b3 K33 src {0,1,2}, pi=(3 5):", [str(x) for x in ratio(6, K33, {3: 5, 5: 3}, L, L)])
    l = {1: [1, 2], 2: [-3, 1], 3: [Fr(1, 2), -1], 5: [2, 2]}
    L = lamlist(6, l, 2)
    print("b3 K33 src {1,2,3,5} (+root 5 implied), pi=(0 4):", [str(x) for x in ratio(6, K33, {0: 4, 4: 0}, L, L)])
    # b5: C^2 n=3 graph, sources at 3,2,4 (lam0 = lam1 = 0), pi = (0 1)
    l = {3: [1, 0], 2: [Fr(1, 3), Fr(2, 5)]}
    L = lamlist(5, l, 2)
    print("b5 src {3,2,4}, pi=(0 1):", [str(x) for x in ratio(5, C2N3, {0: 1, 1: 0}, L, L)])
    # b2: prism sources {0,1,2} (lam2 = -(e1+e2)) vs sources {0,1,3} (lam3 = -(e1+e2)), pi = (2 3)(4 5)
    La = lamlist(6, {0: [1, 0], 1: [0, 1], 2: [-1, -1]}, 2)
    Lb = lamlist(6, {0: [1, 0], 1: [0, 1], 3: [-1, -1]}, 2)
    print("b2 prism: D(Pt, src{0,1,3})/D(t, src{0,1,2}), pi=(2 3)(4 5):", [str(x) for x in ratio(6, PRISM, {2: 3, 3: 2, 4: 5, 5: 4}, La, Lb)])
