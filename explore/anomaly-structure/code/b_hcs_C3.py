"""Part (b): holomorphic Chern-Simons on C^3 (H = 3, T = 0).

1. Kernels with four legs (engine, symbolic sources at the external vertices, lambda = 0 at the
   trivalent vertices), reduced by momentum conservation sum_i p_i = 0:
     box (n = 1), theta(1,3,3), theta(2,2,3) (n = 2).
   Independent pointwise brute-force check of the densities is in b_pointwise.py.
2. S_4 behaviour of Delta = det(p_1, p_2, p_3) modulo conservation.
3. Colour weights (orientation rule of cochain.py), their totally (anti)symmetric parts.
4. dim (Lambda^4 g^*)^g and dim (S^4 g^*)^g, exact, for sl2, gl2, sl3, gl3, sl4.
5. The obstruction cochains omega_n = sum_Gamma |Aut|^{-1} K_Gamma W_Gamma (alternated, modulo
   conservation): n = 2 vanishes identically; n = 1 is nonzero and closed; exactness rank of n = 1.
"""
import sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import theta, polygon
import lie, cochain as cc

H = 3
mom = cc.Momenta(6, H)
P = mom.p


def det3(a, b, c):
    return sp.Matrix([a, b, c]).T.det()


Delta = det3(P[0], P[1], P[2])


def kernel_4legs(G, ext):
    """K with symbolic sources p_i at ext[i], zero elsewhere, as a polynomial in p_0..p_3."""
    lam = {}
    for v in range(G.nv - 1):
        lam[v] = list(P[ext.index(v)]) if v in ext else [0] * H
    K = htk.kernel_exact(G, lam, T=0)
    # root source = -sum; if the root is external, its momentum is p_{ext.index(root)} = -sum others
    return sp.expand(K)


def reduce_mod_conservation(expr, k=4):
    return sum(sp.Rational(c) * mom.mono_expr(m, k, [P[i] for i in range(k - 1)])
               for m, c in mom.canon(expr, k).items())


print("=== 1. kernels with four legs, H = 3, T = 0 ===")
box = polygon(4)
graphs = {}
# box: all four vertices external (root 3 carries -sum)
Kb = kernel_4legs(box, [0, 1, 2]) if False else None
lam = {0: list(P[0]), 1: list(P[1]), 2: list(P[2])}
Kb = sp.expand(htk.kernel_exact(box, lam, T=0))
Kb_red = sp.factor(reduce_mod_conservation(Kb))
print("box: edges", box.edges, " K(p0,p1,p2; root = -sum) =", sp.factor(Kb), "; mod conservation:", Kb_red)
graphs['box'] = dict(G=box, ext=[0, 1, 2, 3], K=Kb, aut=4 * 2 // 1)
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, interior = theta(*pqr)
    ext = [v for Pp in interior for v in Pp]
    K = kernel_4legs(G, ext)
    Kred = sp.factor(reduce_mod_conservation(K))
    print("theta%s: edges %s, external vertices %s (path order A,B,C)" % (pqr, G.edges, ext))
    print("   K =", sp.factor(K))
    print("   mod conservation: K =", Kred, ";  K / Delta^2 =", sp.simplify(Kred / Delta ** 2))
    graphs['theta%d%d%d' % pqr] = dict(G=G, ext=ext, K=K)
graphs['theta133']['aut'] = 4
graphs['theta223']['aut'] = 4
graphs['box']['aut'] = 8
print("|Aut|: box 8 (dihedral), theta(1,3,3) 4, theta(2,2,3) 4")

print("\n=== 2. Delta under S_4 modulo conservation ===")
for sigma in itertools.permutations(range(4)):
    D_s = det3(P[sigma[0]], P[sigma[1]], P[sigma[2]])
    diff = reduce_mod_conservation(sp.expand(D_s - cc.perm_sign(sigma) * Delta))
    assert diff == 0
print("Delta(p_s1, p_s2, p_s3) = sgn(s) Delta(p_1, p_2, p_3) mod conservation for all 24 s: True")

algs = [lie.sl(2), lie.gl(2), lie.sl(3), lie.gl(3)]


def sym_parts(W):
    k = W.ndim
    A = sum(cc.perm_sign(s) * np.transpose(W, cc.transpose_for(s)) for s in itertools.permutations(range(k)))
    S = sum(np.transpose(W, cc.transpose_for(s)) for s in itertools.permutations(range(k)))
    return A, S


print("\n=== 3. colour weights: 24 x (anti)symmetric parts (integer, kinv scaled) ===")
colour = {}
for g in algs:
    for name, d in graphs.items():
        W, scale = cc.colour_weight(g, d['G'].nv, d['G'].edges, d['ext'])
        A, S = sym_parts(W)
        colour[(g.name, name)] = W
        print("%s %-9s |W|_max %6d  antisym part zero: %s   sym part zero: %s   (kinv scale^|I| = %d)" % (
            g.name, name, np.abs(W).max(), not A.any(), not S.any(), scale))


print("\n=== 4. invariant alternating / symmetric 4-tensors ===")


def rank_mod_p(A, p=2147483629):
    """rank of an integer matrix over GF(p); rank over Q >= rank over GF(p)."""
    M = np.array(A, dtype=np.int64) % p
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if piv.size == 0:
            continue
        i = r + piv[0]
        M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for j in nz:
            if j != r:
                M[j] = (M[j] - M[j, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r


def invariant_alt_dim(g):
    """dim (Lambda^4 g*)^g: nullity of the invariance equations on the basis e^a^e^b^e^c^e^d,
    a<b<c<d.  Returns (nullity over GF(p), exact nullity or None)."""
    d = g.dim
    basis = list(itertools.combinations(range(d), 4))
    if not basis:
        return 0, 0
    index = {b: i for i, b in enumerate(basis)}
    # (delta_e w)(X_a,X_b,X_c,X_d) = -sum_slots w(.., [X_e, X_slot], ..);  linear in the coordinates
    A = np.zeros((d * len(basis), len(basis)), dtype=np.int64)
    for e in range(d):
        ad = g.fu[e]
        for ri, (a, b, c, dd) in enumerate(basis):
            row = e * len(basis) + ri
            tup = [a, b, c, dd]
            for slot in range(4):
                for x in np.nonzero(ad[tup[slot]])[0]:
                    new = tup[:]
                    new[slot] = int(x)
                    if len(set(new)) < 4:
                        continue
                    sgn = cc.perm_sign(new)
                    A[row, index[tuple(sorted(new))]] += sgn * int(ad[tup[slot], x])
    A = A[np.any(A != 0, axis=1)]
    rp = rank_mod_p(A)
    null_p = len(basis) - rp
    exact = None
    if null_p > 0:
        from sympy.polys.matrices import DomainMatrix
        from sympy import ZZ
        DM = DomainMatrix([[ZZ(int(x)) for x in row] for row in A.tolist()], A.shape, ZZ)
        exact = len(basis) - DM.rank()
    else:
        exact = 0
    return null_p, exact


for g in [lie.sl(2), lie.gl(2), lie.sl(3), lie.gl(3), lie.sl(4), lie.gl(4)]:
    t0 = time.time()
    null_p, exact = invariant_alt_dim(g)
    print("%s: dim(Lambda^4 g*)^g = %s  (nullity mod p %d; centre dim %d; %.1fs)" % (
        g.name, exact, null_p, g.center_dim(), time.time() - t0), flush=True)


print("\n=== 5. obstruction cochains (canonical form modulo conservation) ===")
# integer scaling of the kernels: box K = c1 Delta-type, thetas c2 ...
def int_scale(expr):
    expr = sp.expand(expr)
    den = sp.ilcm(*[sp.Rational(c).q for c in sp.Poly(expr, *[x for v in P[:4] for x in v]).coeffs()])
    return den


for g in algs:
    # n = 2: sum_Gamma K_Gamma W_Gamma / |Aut|; common factor: 1/4 and lcm of denominators
    s133 = int_scale(graphs['theta133']['K'])
    s223 = int_scale(graphs['theta223']['K'])
    s = sp.ilcm(s133, s223)
    w133 = cc.alternate(mom, graphs['theta133']['K'] * s, colour[(g.name, 'theta133')], 4)
    w223 = cc.alternate(mom, graphs['theta223']['K'] * s, colour[(g.name, 'theta223')], 4)
    print("%s n=2: omega(theta133) == 0: %s; omega(theta223) == 0: %s" % (g.name, cc.is_zero(w133), cc.is_zero(w223)))
    # n = 1 box
    sb = int_scale(graphs['box']['K'])
    wb = cc.alternate(mom, graphs['box']['K'] * sb, colour[(g.name, 'box')], 4)
    dwb = cc.d_CE(mom, g, wb, 4)
    print("%s n=1: omega(box) == 0: %s; d omega(box) == 0: %s; monomials %d" % (
        g.name, cc.is_zero(wb), cc.is_zero(dwb), len(wb)), flush=True)
