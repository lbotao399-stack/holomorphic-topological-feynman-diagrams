"""(H+1)-gon kernel at (H,0): K = c_H det(lam_0..lam_{H-1}), c_H = (-1)^H / H!.

(a) own Kirchhoff layer, density D(t) from every omitted edge, exact; (b) htk.kernel_exact (engine, own
self-check against its brute layer); (c) alternation identity sgn(s) det(p_s(0..H-1)) = det(p_0..p_{H-1});
(d) cross-checks with existing project pieces (fan(1)/LLC, box/S2, msym-D8 pentagon edge order,
type-i hexagon, banana(2)).
"""
import sys
from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial
import random

from common import Log, det, perm_sign, P

sys.path.insert(0, str(P / "code"))
import htk  # noqa: E402
from htk.graph import polygon, Graph, fan, box  # noqa: E402

log = Log("polygon.log")
rng = random.Random(20261127)


def own_density(H, t, lam):
    """polygon m = H+1: edges e_j = (j, j+1 mod m), root = H; returns list of D(t) from each omitted edge"""
    m = H + 1
    U = sum(t)
    Q = []
    acc = [0] * H
    for e in range(H):
        acc = [a + b for a, b in zip(acc, lam[e])]
        Q.append(list(acc))
    Q.append([0] * H)                                   # Q_H = sum of all sources = 0 at the root edge
    ell = [-sum(F(t[e]) * Q[e][i] for e in range(m)) / U for i in range(H)]
    y = [[ell[i] + Q[e][i] for i in range(H)] for e in range(m)]
    # Kirchhoff checks: B^T y = lam (non-root vertices), sum t y = 0 (cycle C = (1,..,1))
    for v in range(H):
        out_e, in_e = v, (v - 1) % m
        assert [y[out_e][i] - y[in_e][i] for i in range(H)] == [F(x) for x in lam[v]]
    assert all(sum(F(t[e]) * y[e][i] for e in range(m)) == 0 for i in range(H))
    Ds = []
    for k in range(m):
        cols = [e for e in range(m) if e != k]
        coef = det([[y[e][i] for e in cols] for i in range(H)]) / U ** H   # [Omega]_{dt_cols}
        Ds.append(F((-1) ** k) / t[k] * coef)
    return Ds, U


log("(a) own layer: D(t) U^{H+1} = (-1)^H det(lam_0..lam_{H-1}) for every omitted edge")
for H in range(1, 9):
    ok = True
    for trial in range(4):
        t = [F(rng.randint(1, 50), rng.randint(1, 9)) for _ in range(H + 1)]
        lam = [[rng.randint(-7, 7) for _ in range(H)] for _ in range(H)]
        Ds, U = own_density(H, t, lam)
        Delta = det([[lam[v][i] for v in range(H)] for i in range(H)])
        target = F((-1) ** H) * Delta / U ** (H + 1)
        ok &= all(D == target for D in Ds)
    log("  H=%d: %s  =>  K = (-1)^H Delta * vol(Delta^H) = %s * Delta" % (H, ok, F((-1) ** H, factorial(H))))

log("(b) htk.kernel_exact(polygon(H+1)) / det(lam_0..lam_{H-1})")
for H in range(1, 7):
    vals = []
    for trial in range(2):
        lam = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(H)]
        Delta = det([[lam[v][i] for v in range(H)] for i in range(H)])
        if Delta == 0:
            continue
        K = htk.kernel_exact(polygon(H + 1), lam, T=0, check=1)
        vals.append(F(K) / Delta)
    log("  H=%d: c_H = %s   (expected %s): %s" % (H, vals, F((-1) ** H, factorial(H)),
                                                  all(v == F((-1) ** H, factorial(H)) for v in vals)))

log("(c) sum p = 0  =>  sgn(s) det(p_s(0),..,p_s(H-1)) = det(p_0..p_{H-1}) for all s in S_{H+1}")
for H in range(1, 6):
    p = [[rng.randint(-6, 6) for _ in range(H)] for _ in range(H)]
    p.append([-sum(c) for c in zip(*p)])
    D0 = det([[p[v][i] for v in range(H)] for i in range(H)])
    ok = all(perm_sign(s) * det([[p[s[v]][i] for v in range(H)] for i in range(H)]) == D0
             for s in permutations(range(H + 1)))
    log("  H=%d: %s (|S_%d| = %d)" % (H, ok, H + 1, factorial(H + 1)))

log("(d) existing pieces")
lam2 = [[3, -2], [5, 7]]
d2 = det([[3, 5], [-2, 7]])
Kf = htk.kernel_exact(fan(1), lam2, check=1)
Kp = htk.kernel_exact(polygon(3), lam2, check=1)
log("  H=2 fan(1) (LLC triangle) K = %s = %s det;  polygon(3) K = %s = %s det;  I^S1 = -K_fan = %s det"
    " (LLC eq-1-42: det/2)" % (Kf, F(Kf) / d2, Kp, F(Kp) / d2, -F(Kf) / d2))
lam3 = [[1, 2, 0], [0, 1, 3], [2, -1, 1]]
d3 = det([[lam3[v][i] for v in range(3)] for i in range(3)])
Kb = htk.kernel_exact(box(), lam3, check=1)
log("  H=3 box() (S2, edges (3,0),(0,1),(1,2),(2,3)) K = %s det  (anomalies REPORT c_box = -1/6)" % (F(Kb) / d3))
pent = Graph(5, [(4, 0), (0, 1), (1, 2), (2, 3), (3, 4)], "pentagon_msymD8")
lam4 = [[1, 0, 2, 0], [0, 1, 1, 3], [2, 1, 0, 1], [1, -1, 1, 1]]
d4 = det([[lam4[v][i] for v in range(4)] for i in range(4)])
Kq = htk.kernel_exact(pent, lam4, check=1)
log("  H=4 msym-D8 pentagon edge order K = %s det  (msym-D8 VERIFY C5: <0123>/24)" % (F(Kq) / d4))
lam5 = [[int(i == j) + (i == 4 and j == 0) for j in range(5)] for i in range(5)]
d5 = det([[lam5[v][i] for v in range(5)] for i in range(5)])
Kh = htk.kernel_exact(polygon(6), lam5, check=1)
log("  H=5 polygon(6) K = %s det  (type-i README: -det/120)" % (F(Kh) / d5))
Kb2 = htk.kernel_exact(Graph(2, [(0, 1), (0, 1)], "banana2"), [[4]], check=1)
Kp2 = htk.kernel_exact(polygon(2), [[4]], check=1)
log("  H=1 polygon(2) K = %s * lam_0;  banana(2) (one edge reversed, (-1)^D) K = %s * lam_0" % (F(Kp2) / 4, F(Kb2) / 4))
log("c_H = (-1)^H/H!:", [str(F((-1) ** H, factorial(H))) for H in range(1, 9)])
log.close()
