"""d_2 = [Q_2] on H(Q_0) for W = kappa = 0 theories, sectors (J0, f, af) -> (J0, f, af-2).
Usage: python3 lift.py <theory> <J0max>.  Log: lift_<theory>.log"""
import sys, time
from fractions import Fraction as Fr
from ht3d import Engine, load_kernels, pstr
from cohom import sector_basis, Echelon, kernel_basis, poly_to_vec, vec_to_poly
from theories import *

TH = {"su2adj": lambda: su2_adj(1), "su2pure": lambda: su2_adj(0),
      "u1": lambda: u1([1, -1]), "u1q": lambda: u1([1, 2, -3])}
name, J0max = sys.argv[1], int(sys.argv[2])
th = TH[name]()
E = Engine(th, load_kernels(4))
print("theory", th.name, "W = 0, kappa = 0")

cache = {}
def sector(J0, f, af):
    key = (J0, f, af)
    if key in cache:
        return cache[key]
    bas = sector_basis(th, J0, f, af)
    idx = {g: {m: i for i, m in enumerate(b)} for g, b in bas.items()}
    # Q0 images and kernels
    Z, Bim = {}, {}
    for g, b in bas.items():
        imgs = []
        for m in b:
            P = E.Q0({m: Fr(1)})
            imgs.append(poly_to_vec(P, idx.get(g + 1, {})) if P else {})
        Z[g] = kernel_basis(None, imgs, None)
        ech = Echelon()
        for v in imgs:
            if v:
                ech.add(v)
        Bim[g + 1] = ech  # image in degree g+1
    H = {g: len(Z[g]) - (Bim[g].rank() if g in Bim else 0) for g in bas}
    cache[key] = (bas, idx, Z, Bim, H)
    return cache[key]

t0 = time.time()
total = {}
for J0 in range(1, J0max + 1):
    for af in range(0, J0 + 1):
        for f in range(-af, J0 + af + 1):
            bas, idx, Z, Bim, H = sector(J0, f, af)
            nz = {g: h for g, h in H.items() if h}
            if not nz:
                continue
            line = "J0=%d f=%+d af=%d: dim C = %s ; dim H(Q0) = %s" % (J0, f, af, {g: len(b) for g, b in sorted(bas.items())}, dict(sorted(nz.items())))
            if af >= 2:
                tb, tidx, tZ, tBim, tH = sector(J0, f, af - 2)
                ranks = {}
                for g, h in nz.items():
                    if not tH.get(g + 1):
                        continue
                    base = tBim.get(g + 1, Echelon())
                    ech = Echelon(); ech.rows = dict(base.rows)
                    r0 = ech.rank()
                    reps = []
                    for zv in Z[g]:
                        P = vec_to_poly(zv, bas[g])
                        Q2P = E.Q2(P)
                        if Q2P:
                            v = poly_to_vec(Q2P, tidx[g + 1])
                            if ech.add(v):
                                reps.append((P, Q2P))
                    if ech.rank() - r0:
                        ranks[g] = ech.rank() - r0
                        for P, Q2P in reps[:2]:
                            print("    lifted: O =", pstr(P, 8), "\n       Q2 O =", pstr(Q2P, 8))
                if ranks:
                    line += " ; rank d2 (gh -> gh+1, af -> af-2) = %s" % ranks
            print(line); sys.stdout.flush()
print("time %.0fs" % (time.time() - t0))
