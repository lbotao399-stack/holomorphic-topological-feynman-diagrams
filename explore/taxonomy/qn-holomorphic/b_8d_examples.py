"""Fork B: 8d (C^{4|1}, sl2) explicit Q_1 = (1/2) hatQ_pentagon on the lowest lifted class [tr gamma gamma]
(sector #gamma = 2, d = 0, k = 2; run_lifts b_lifts_8d_sl2_k4_d1.log) and its image class; structure check of
hatQ_1(gamma^A gamma^B) against  sum (pentagon kernel) x (path contraction V G V G V G V).  Log: b_8d_examples.log"""
import itertools
from fractions import Fraction as Fr
import numpy as np
import polyops as po
import theories
from wick import Algebra
from explicit import graph_ops
import lifts as lf

H, sl = theories.build('8d', 'sl2')
alg = Algebra(sl, H)
g = theories.gauge('sl2')
names = ['c[%s]' % n for n in ('E12', 'E21', 'H')] + ['gam[%s]' % n for n in ('E12', 'E21', 'H')]
cartan = [2]
pos = [0]
th = lf.LiftTheory(alg, cartan=[2], posroots=[0])
G1 = graph_ops(alg, 1)[0]
k1 = po.modq(Fr(1, 2))
monos = lf.sector_monomials(th, (2,), (0, 0, 0, 0), 2)
inv = lf.invariants(th, monos)
O = lf.vec_to_poly(inv[0], monos)
print('source class (sector #gamma=2, d=0, k=2):')
print(po.show(O, names, H))
img = lf.drop_c0(th, po.pscale(G1.apply(O), k1))
print('Q_1 O  (%d monomials), first terms:' % len(img))
print(po.show(img, names, H, maxterms=12))
# derivative structure: every monomial has one d_a on each of four letters with a permutation of (1,2,3,4)
perm_ok = all(sorted(tuple(al) for _, al in m) == sorted([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]) for m in img)
print('all monomials: four letters with d_1, d_2, d_3, d_4 once each:', perm_ok)
# antisymmetry in the derivative labels: permuting the derivative directions by sigma multiplies by sgn(sigma)
def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s
ok = True
for sg in itertools.permutations(range(4)):
    new = {}
    for m, v in img.items():
        word = [(A, tuple(al[sg.index(a)] for a in range(4))) for A, al in m]
        s, mm = po.canon(word, alg.pi)
        po.padd(new, mm, s * v * perm_sign(sg))
    if po.psum(new, po.pscale(img, -1)):
        ok = False
print('Q_1 O is totally antisymmetric in the four derivative directions (eps^{abcd} structure):', ok)
cnt = {}
for m in img:
    key = tuple(sorted(names[A].split('[')[0] for A, _ in m))
    cnt[key] = cnt.get(key, 0) + 1
print('letter content of Q_1 O:', cnt)

# closed form: compare with T = eps^{abcd} tr( d_d gamma [d_a c, [d_b c, d_c c]] ) (fundamental trace of sl2 matrices)
mats = g.mats


import sympy as _sp
_kinv = _sp.Matrix(g.kappa.tolist()).inv()


def mat_letter(block_off, a_dir):
    """c: sum_A E_A c^A;  gamma (dual block): sum_{A,B} (kappa^{-1})^{AB} E_B gamma_A  (g^v ~ g by kappa)."""
    M = [[{} for _ in range(2)] for _ in range(2)]
    al = tuple(1 if i == a_dir else 0 for i in range(4))
    for A in range(3):
        coefs = {A: 1} if block_off == 0 else {B: _kinv[A, B] for B in range(3) if _kinv[A, B] != 0}
        for B, cf in coefs.items():
            cf = po.modq(Fr(int(_sp.fraction(cf)[0]), int(_sp.fraction(cf)[1])))
            for i in range(2):
                for j in range(2):
                    if mats[B][i, j]:
                        po.padd(M[i][j], ((block_off + A, al),), cf * int(mats[B][i, j]))
    return M


def pmul(p, q):
    out = {}
    for m1, v1 in p.items():
        for m2, v2 in q.items():
            s, m = po.canon(m1 + m2, alg.pi)
            if s:
                po.padd(out, m, s * v1 * v2)
    return out


def mmul(X, Y):
    return [[po.psum(*[pmul(X[i][k], Y[k][j]) for k in range(2)]) for j in range(2)] for i in range(2)]


def msub(X, Y):
    return [[po.psum(X[i][j], po.pscale(Y[i][j], -1)) for j in range(2)] for i in range(2)]


def comm(X, Y):          # graded: entries of odd-letter matrices anticommute; [X,Y] = XY + YX for odd X, Y
    return [[po.psum(a, b) for a, b in zip(r1, r2)] for r1, r2 in zip(mmul(X, Y), mmul(Y, X))]


T = {}
for sg in itertools.permutations(range(4)):
    a, b, c_, d = sg
    Ca, Cb, Cc = mat_letter(0, a), mat_letter(0, b), mat_letter(0, c_)
    Gd = mat_letter(3, d)
    outer = mmul(Ca, mmul(Cb, Cc))                      # plain matrix product (the [ , [ , ]] form vanishes under eps)
    tr = po.psum(*[pmul(Gd[i][k], outer[k][i]) for i in range(2) for k in range(2)])
    for m, v in tr.items():
        po.padd(T, m, perm_sign(sg) * v)
ratios = set()
for m in set(T) | set(img):
    a_, b_ = img.get(m, 0), T.get(m, 0)
    if b_ == 0:
        ratios.add('img-only' if a_ else None)
    else:
        ratios.add(po.ratrec(a_ * pow(b_, po.P - 2, po.P) % po.P))
ratios.discard(None)
print('T terms %d; Q_1 tr(gamma gamma) / [eps^{abcd} tr(d_d gamma d_a c d_b c d_c c)] :' % len(T), ratios)


def eps_sum(builder):
    out = {}
    for sg in itertools.permutations(range(4)):
        a, b, c_, d = sg
        X, Y, Z = mat_letter(0, a), mat_letter(0, b), mat_letter(0, c_)
        Gm = mat_letter(3, d)
        for m, v in builder(X, Y, Z, Gm).items():
            po.padd(out, m, perm_sign(sg) * v)
    return out


def tr(M):
    return po.psum(M[0][0], M[1][1])


cands = {
    'tr(g X Y Z)': lambda X, Y, Z, Gm: tr(mmul(Gm, mmul(X, mmul(Y, Z)))),
    'tr(X g Y Z)': lambda X, Y, Z, Gm: tr(mmul(X, mmul(Gm, mmul(Y, Z)))),
    'tr(X Y g Z)': lambda X, Y, Z, Gm: tr(mmul(X, mmul(Y, mmul(Gm, Z)))),
    'tr(g X) tr(Y Z)': lambda X, Y, Z, Gm: pmul(tr(mmul(Gm, X)), tr(mmul(Y, Z))),
}
cv = {k: eps_sum(f) for k, f in cands.items()}
keys = list(cv)
monos = sorted(set(img).union(*[set(v) for v in cv.values()]))
rows = [[cv[k].get(m, 0) for k in keys] for m in monos]
rhs = [img.get(m, 0) for m in monos]
from nilp_explicit import rank_solve
rk, aug, sol = rank_solve(rows, rhs, po.P)
print('fit Q_1 tr(gamma gamma) = sum_i x_i eps^{abcd} T_i (X=d_a c, Y=d_b c, Z=d_c c, g=d_d gamma, fundamental traces):')
print('  rank %d, augmented %d (candidates %s)' % (rk, aug, keys))
if sol is None:
    # least: solve on independent subset
    import sympy
    M = sympy.Matrix([[po.ratrec(x) if x else 0 for x in r] for r in rows])
    print('  candidate rank over Q:', M.rank())
else:
    print('  solution:', {k: str(po.ratrec(s)) for k, s in zip(keys, sol)})
