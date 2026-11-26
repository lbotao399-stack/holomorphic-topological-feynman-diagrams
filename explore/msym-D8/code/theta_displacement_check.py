"""Displacement map of theta(p,q,r) with O at a pole, by direct Kirchhoff solution at random rational t:
    - sum_e y_e . z_e = sum_v lam_v . X_v,   X_v = sigma_v z_{j(v)} + (1 - sigma_v) Z_N,
    Z_N = sum_j c_j z_j,  c_j = (1/R_j) / sum_i (1/R_i),  R_j = sum_{e in P_j} t_e,
    sigma_v = (sum of t_e on P_j between v and N) / R_j   (z_j on the spoke of P_j, oriented away from O).
Pentagon (S2 Fig. fig-box analogue): X_j = r_j z_1 + (1 - r_j) z_2, r_j = (1/T) sum_{k>j} t_k."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
from graphs_c4 import theta_O, pentagon, THETAS
import brute_c4 as bf
H = 4
rng = random.Random(17)
for pqr in THETAS:
    G, spokes, outputs, paths = theta_O(*pqr)
    ok = True
    for _ in range(5):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(G.ne)]
        lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in outputs}
        z = {j: [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(H)] for j in (1, 2, 3)}
        y, _ = bf.kirchhoff(G.nv, G.edges, t, lam)
        lhs = -sum(s * sum(y[e][a] * z[j][a] for a in range(H)) for e, (j, s) in spokes.items())
        R = [sum(t[e] for e in P) for P, _ in paths]
        cj = [(1 / Rj) / sum(1 / Ri for Ri in R) for Rj in R]
        ZN = [sum(cj[j] * z[j + 1][a] for j in range(3)) for a in range(H)]
        rhs = Fr(0)
        for j, (P, ints) in enumerate(paths):
            for k, v in enumerate(ints):
                sig = sum(t[e] for e in P[k + 1:]) / R[j]
                X = [sig * z[j + 1][a] + (1 - sig) * ZN[a] for a in range(H)]
                rhs += sum(lam[v][a] * X[a] for a in range(H))
        ok &= (lhs == rhs)
    print('theta%s: -sum y.z == sum lam_v.(sigma_v z_j + (1-sigma_v) Z_N) at 5 random rational points:' % (pqr,), ok)
G, spokes, outputs = pentagon()
ok = True
for _ in range(5):
    t = [Fr(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(5)]
    lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in outputs}
    z = {j: [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(H)] for j in (1, 2)}
    y, _ = bf.kirchhoff(G.nv, G.edges, t, lam)
    lhs = -sum(s * sum(y[e][a] * z[j][a] for a in range(H)) for e, (j, s) in spokes.items())
    T = sum(t)
    rhs = sum(lam[v][a] * ((sum(t[v + 1:]) / T) * z[1][a] + (1 - sum(t[v + 1:]) / T) * z[2][a]) for v in outputs for a in range(H))
    ok &= lhs == rhs
print('pentagon: -sum y.z == sum lam_j.(r_j z_1 + (1-r_j) z_2) at 5 random rational points:', ok)
