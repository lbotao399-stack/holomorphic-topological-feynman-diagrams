"""Independent check (single_particle.py functions only) of symmetry/REPORT.txt 2.1:
  rank( phi = del(F . -) : Z(A0)(W) -> A2(W - 1^5) ) = B^1(gamma(W)),  gamma(W) = (W_z, -W_w),
  B^1(gamma) = [gamma_z >= 0, gamma_theta in {0,1}^3, gamma != 0];
  dim ker phi(beta - 1^5) - coker(Wn(beta)) = [beta in M1] + [beta = 1^5],  Wn(beta) = (beta_z - 1, -beta_theta - 1).
usage: python3 c5_symmetry_spotcheck.py zmax tmax"""
import sys
from itertools import product
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twisted-holography-N4/code")
import single_particle as sp

def phi_data(W):
    A0 = sp.bulk_basis_A0(W)
    k0 = sp.kernel_basis([sp.div_A0(e) for e in A0], None)
    vecs = []
    for comb in k0:
        v = {}
        for i, c in comb.items():
            for k, cc in sp.F_times(A0[i]).items():
                v[k] = v.get(k, 0) + c * cc
        dv = {}
        for k, c in v.items():
            if not c: continue
            for kk, cc in sp.div_A2(k).items():
                dv[kk] = dv.get(kk, 0) + c * cc
        vecs.append({k: c for k, c in dv.items() if c})
    r = sp.rank(vecs)
    return len(k0), r

def B1(g):
    return int(min(g[:2]) >= 0 and all(x in (0, 1) for x in g[2:]) and any(g))

zmax, tmax = int(sys.argv[1]), int(sys.argv[2])
bad = 0; n = 0
for bz in product(range(zmax + 1), repeat=2):
    if bz[0] > bz[1]: continue
    for bt in product(range(tmax + 1), repeat=3):
        if not (bt[0] <= bt[1] <= bt[2]): continue
        beta = bz + bt
        Wr = tuple(b - 1 for b in beta)
        gam = (Wr[0], Wr[1], -Wr[2], -Wr[3], -Wr[4])
        z0, r = phi_data(Wr)
        Wn = (bz[0] - 1, bz[1] - 1) + tuple(-b - 1 for b in bt)
        cok = sp.bulk_coker(Wn)[0]
        M1 = int(any(beta) and all(x in (0, 1) for x in bt) and min(bz) == 0)
        one = int(beta == (1, 1, 1, 1, 1))
        ok = (r == B1(gam)) and (z0 - r - cok == M1 + one)
        n += 1
        if not ok:
            bad += 1
            print("MISMATCH", beta, "dimZ0", z0, "rank phi", r, "B1(gamma)", B1(gam), "ker", z0 - r, "coker", cok, "M1", M1, "1^5", one)
print(f"weights {n} (beta_z in [0,{zmax}]^2, beta_theta in [0,{tmax}]^3 mod S2xS3), mismatches {bad}")
