"""Comparison of the closed two-loop formulas with S1 (two-loop fan), S2 (graphs A, B), S3 (five
theta graphs, D_pqr sums, H = 0 remark), S4 (ladder n = 2).  Exact arithmetic, random rational
sources; every kernel also compared with htk.kernel_exact.  Writes sources_results.json."""
import json
import random
import sys

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import *  # noqa
import htk
from htk import conventions as cv
from htk.simplex import PiNumber

rng = random.Random(7)


def rv(H):
    return [Fr(rng.randint(-9, 9), rng.randint(1, 3)) for _ in range(H)]


def d2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def d3(a, b, c):
    return det([[a[i], b[i], c[i]] for i in range(3)])


def dn(*vs):
    return det([[v[i] for v in vs] for i in range(len(vs))])


out = {}

# ------------------------------------------------------------------ S1: two-loop fan (H = 2)
G = htk.graph.fan(2)                         # path 0 -> 1 -> 2, root 3 (S1 vertex i = htk i-1)
th = Theta(G.nv, G.edges)
assert (th.N, th.Np, th.pqr) == (1, 3, (1, 2, 2))
res = []
for case in ("lambda2=0", "general"):
    for _ in range(3):
        l1, l2, l3 = rv(2), rv(2), rv(2)
        if case == "lambda2=0":
            l2 = [Fr(0), Fr(0)]
        lam = full_sources(4, 3, {0: l1, 1: l2, 2: l3})
        K = kernel_htk_closed(th, lam, 2, 0)
        Ke = Fr(htk.kernel_exact(G, [l1, l2, l3], 0))
        I1 = cv.S1_kernel(K, 2)
        f42 = d2(l1, l3) ** 2 / 24                                                   # S1 eq-1-42
        spider = Fr(1, 24) * d2(l1, [a + 2 * b for a, b in zip(l3, l2)]) * d2(l3, [a + 2 * b for a, b in zip(l1, l2)])
        res.append({"case": case, "K_closed_eq_htk_exact": K == Ke, "I_S1": str(I1),
                    "eq-1-42": str(f42), "I_S1 == eq-1-42": I1 == f42,
                    "I_S1 == -spider_two_arm_form": I1 == -spider})
out["S1_fan2"] = res
print("S1 fan n=2:")
for r_ in res:
    print("  ", r_)

# ------------------------------------------------------------------ S2: graphs A, B (H = 3)
res = {}
for name, G, zero, f in (
        ("A", htk.graph.graph_A(), 2,
         lambda L: Fr(1, 432) * (2 * d3(L[0], L[1], L[3]) + d3(L[0], L[1], L[4])) * (d3(L[0], L[3], L[4]) + 2 * d3(L[1], L[3], L[4]))),
        ("B", htk.graph.graph_B(), 2,
         lambda L: -Fr(1, 288) * (2 * d3(L[0], L[1], L[3]) + d3(L[0], L[1], L[4])) * (d3(L[0], L[3], L[4]) + d3(L[1], L[3], L[4])))):
    th = Theta(G.nv, G.edges)
    rows = []
    for case in ("S2 sources", "general"):
        for _ in range(3):
            L = [rv(3) for _ in range(5)]
            if case == "S2 sources":
                L[zero] = [Fr(0)] * 3
            lam = full_sources(6, 5, {v: L[v] for v in range(5)})
            K = kernel_htk_closed(th, lam, 3, 0)
            Ke = Fr(htk.kernel_exact(G, L, 0))
            I2 = cv.S2_kernel(K, 2)
            row = {"case": case, "K_closed_eq_htk_exact": K == Ke, "I_S2": str(I2)}
            if case == "S2 sources":
                row["S2 value"] = str(f(L))
                row["equal"] = I2 == f(L)
            rows.append(row)
    res[name] = {"N": th.N, "Np": th.Np, "pqr": th.pqr, "interior": th.interior, "rows": rows}
    print("S2 graph", name, "N =", th.N, "N' =", th.Np, th.pqr, th.interior)
    for r_ in rows:
        print("  ", r_)
out["S2"] = res

# ------------------------------------------------------------------ S3: H = 5 thetas, D_pqr sums, H = 0
c_S3 = {(1, 6, 6): Fr(1, 158400), (2, 5, 6): Fr(1, 285120), (3, 4, 6): Fr(1, 332640),
        (3, 5, 5): Fr(1, 798336), (4, 4, 5): Fr(1, 997920)}


def Dsum(pqr, I, lam):
    A = [lam[u] for u in I["A"]]
    B = [lam[u] for u in I["B"]]
    C = [lam[u] for u in I["C"]]
    if pqr == (1, 6, 6):
        return dn(*B) * dn(*C)
    if pqr == (2, 5, 6):
        return dn(A[0], *B) * dn(*C)
    if pqr == (3, 4, 6):
        return dn(A[0], A[1], *B) * dn(*C)
    if pqr == (3, 5, 5):
        return sum((-1) ** d * dn(A[0], A[1], *[B[k] for k in range(4) if k != d - 1]) * dn(B[d - 1], *C)
                   for d in range(1, 5))
    if pqr == (4, 4, 5):
        return sum((-1) ** (d + 1) * dn(A[0], A[1], A[2], *[B[k] for k in range(3) if k != d - 1]) * dn(B[d - 1], *C)
                   for d in range(1, 4))


rows = []
for pqr, cS in c_S3.items():
    nv, edges = theta_graph(*pqr)
    th = Theta(nv, edges)
    G = htk.Graph(nv, edges)
    L = [rv(5) for _ in range(nv - 1)]
    lam = full_sources(nv, nv - 1, {v: L[v] for v in range(nv - 1)})
    K = kernel_htk_closed(th, lam, 5, 1)
    KS3 = cv.S3_kernel(K, 2)
    Dt = D_theta(th, lam)
    ok_c = gamma_half_ratio(pqr, 5) == cS
    ok_K = KS3 == cS * Dt
    ok_D = Dsum(pqr, th.interior, lam) == Dt
    ok_e = htk.kernel_exact(G, L, 1) == PiNumber({1: K})
    rows.append({"pqr": pqr, "c_closed": str(gamma_half_ratio(pqr, 5)), "c_S3": str(cS), "c equal": ok_c,
                 "K_S3 == c pi D_theta": ok_K, "S3 explicit D_pqr == D_theta": ok_D, "htk exact": ok_e})
    print("S3", rows[-1])
out["S3"] = rows
K0 = gamma_half_ratio((1, 1, 1), 0)
out["S3_H0"] = {"K/pi": str(K0), "(1/4pi) K": str(K0 / 4)}
print("S3 H=0 remark: K/pi =", K0, " (1/4pi)K =", K0 / 4)

# ------------------------------------------------------------------ S4: ladder n = 2 (root at a corner)
G, cot, a, b = htk.graph.ladder(2)
th = Theta(G.nv, G.edges, root=G.root)
print("S4 ladder(2): N =", th.N, "N' =", th.Np, "pqr", th.pqr, "interior", th.interior, "root", G.root)
rows = []
for _ in range(3):
    A_, B_, C_ = rv(3), rv(3), rv(3)
    lam = {a(0): A_, b(0): B_, a(2): [-(x + y + z) for x, y, z in zip(A_, B_, C_)]}
    lamfull = full_sources(G.nv, G.root, lam)
    assert lamfull[G.root] == C_
    K = kernel_htk_closed(th, lamfull, 3, 0)
    Ke = Fr(htk.kernel_exact(G, lamfull[:-1], 0))
    I4 = cv.S4_integral(K, G, cot, 3)
    Delta = d3(A_, B_, C_)
    rows.append({"K_closed_eq_htk_exact": K == Ke, "I_S4/Delta^2": str(I4 / Delta ** 2),
                 "S4 closed form n=2": "-1/432", "equal": I4 / Delta ** 2 == Fr(-1, 432)})
    print("  ", rows[-1])
out["S4_ladder2"] = rows

with open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/sources_results.json", "w") as f_:
    json.dump(out, f_, indent=1, default=str)
