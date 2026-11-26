"""K_Gamma(lam o pi; z o pi) = sgn(pi_E)^{D+1} sgn(pi_V)^D K_Gamma(lam; z), D = 4 (CONVENTIONS 2.2), for the
automorphisms of theta(1,4,4) (P2 <-> P3: +1) and theta(3,3,3) (P_i <-> P_j: -1), at z = 0 and at first order in z,
from the stored Cauchy-Binet weights; and the closed forms of K_144(lam;0), K_234(lam;0)."""
import sys, glob, pickle, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
res = {}
for f in glob.glob((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/kernels_c4_*.pkl'):
    res.update(pickle.load(open(f, 'rb')))
H = 4
rng = random.Random(2)


def ev(W, outputs, lam):
    tot = Fr(0)
    for m, c in W.items():
        t = c
        for i, v in enumerate(outputs):
            for a in range(H):
                if m[i][a]:
                    t *= Fr(lam[v][a]) ** m[i][a]
        tot += t
    return tot


def K(pqr, lam, z):
    r = res[pqr]
    out = ev(r['weights'][((0,) * H,) * 3], r['outputs'], lam)
    for j in range(3):
        for a in range(H):
            ks = tuple(tuple(1 if (jj == j and aa == a) else 0 for aa in range(H)) for jj in range(3))
            out += Fr(z[j + 1][a]) * ev(r['weights'][ks], r['outputs'], lam)
    return out


def br(lam, S):
    import sympy as sp
    return int(sp.Matrix([lam[v] for v in S]).det())


for pqr, vperm, zperm, sign in [((1, 4, 4), {1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}, {1: 1, 2: 3, 3: 2}, 1),
                                ((3, 3, 3), {1: 3, 2: 4, 3: 1, 4: 2, 5: 5, 6: 6}, {1: 2, 2: 1, 3: 3}, -1),
                                ((3, 3, 3), {1: 1, 2: 2, 3: 5, 4: 6, 5: 3, 6: 4}, {1: 1, 2: 3, 3: 2}, -1)]:
    ok = True
    for _ in range(3):
        lam = {v: [rng.randint(-6, 6) for _ in range(H)] for v in range(1, 7)}
        z = {j: [rng.randint(-4, 4) for _ in range(H)] for j in (1, 2, 3)}
        lam2 = {v: lam[vperm[v]] for v in lam}
        z2 = {j: z[zperm[j]] for j in z}
        ok &= K(pqr, lam2, z2) == sign * K(pqr, lam, z)
    print('theta%s, outputs %s, inputs %s: K(lam o pi; z o pi) = %+d K(lam; z) (z = 0 and first order, 3 points): %s'
          % (pqr, vperm, zperm, sign, ok))
ok1 = ok2 = True
for _ in range(4):
    lam = {v: [rng.randint(-6, 6) for _ in range(H)] for v in range(1, 7)}
    z0 = {j: [0] * H for j in (1, 2, 3)}
    A = br(lam, (1, 2, 3, 4)) + 2 * br(lam, (1, 2, 3, 5)) + 3 * br(lam, (1, 2, 3, 6))
    ok1 &= K((1, 4, 4), lam, z0) == Fr(-1, 11520) * A * (br(lam, (1, 4, 5, 6)) + 2 * br(lam, (2, 4, 5, 6)) + 3 * br(lam, (3, 4, 5, 6)))
    ok2 &= K((2, 3, 4), lam, z0) == Fr(-1, 17280) * A * (3 * br(lam, (1, 4, 5, 6)) + 2 * br(lam, (2, 4, 5, 6)) + 4 * br(lam, (3, 4, 5, 6)))
print('K_144(lam;0) = -(1/11520)(<1234>+2<1235>+3<1236>)(<1456>+2<2456>+3<3456>):', ok1)
print('K_234(lam;0) = -(1/17280)(<1234>+2<1235>+3<1236>)(3<1456>+2<2456>+4<3456>):', ok2)
