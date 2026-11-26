"""Source normalizations (C7 tail, C10, C11, C8 S4 ladder, C5 S3 items), from the definitions in S1, S2, S3, S4.

S1: eq-operatope integrand  wedge_{e not in T} (-dy^1_e ^ dy^2_e)  (density w.r.t. iota_E vol), and the explicit
    fan theorem eq-1-40/eq-1-38 at n = 2 integrated with sympy.
S2: eq-kernel-def integrand  eps_T det(B_T) wedge_{e not in T} d^3 y_e, for every spanning tree T;
    general-source values re-derived from S2's moves eq-chain/eq-bubble + box (eq-A-currents, eq-B-currents).
S4: check_ladder convention  wedge_{c in cotree} d^3 y_c, ladder(2), corner sources.
S3: explicit D_166..D_445 (thm-theta), c_pqr (R-c), varsigma_theta = -1 (proof of thm-theta).
"""
import sys, time, random, math, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.dont_write_bytecode = True
from vcore import *
import sympy as sp

rng = random.Random(4711)
t0 = time.time()


def spanning_trees(nv, edges):
    ne = len(edges)
    out = []
    for T in itertools.combinations(range(ne), nv - 1):
        par = list(range(nv))

        def f(x):
            while par[x] != x:
                x = par[x]
            return x
        ok = True
        for e in T:
            a, b = f(edges[e][0]), f(edges[e][1])
            if a == b:
                ok = False
                break
            par[a] = b
        if ok:
            out.append(list(T))
    return out


def dyform_density(nv, edges, t, lam, H, cot, sign_per_edge=1, root=None):
    """[wedge_{c in cot} (sign_per_edge * dy^1_c ^ ... ^ dy^H_c)](v) / det[t; v]  with two random vector sets."""
    y = kirchhoff(nv, edges, t, lam, root=root)
    dy = kirchhoff_dt(nv, edges, t, y, root=root)
    ne = len(edges)
    vals = []
    for _ in range(2):
        while True:
            vecs = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(ne)] for _ in range(ne - 1)]
            dd = det([list(map(Fr, t))] + vecs)
            if dd != 0:
                break
        rows = [[matvec(dy[c][i], v) for v in vecs] for c in cot for i in range(H)]
        vals.append(sign_per_edge ** len(cot) * det(rows) / dd)
    assert vals[0] == vals[1]
    return vals[0]


def K_C3_T0(th, lam, H):
    fac = Fr(1)
    for n_ in th.pqr:
        fac /= math.factorial(n_ - 1)
    c0 = fac * Fr(math.prod(math.factorial(H - n_) for n_ in th.pqr), math.factorial(H + 1))
    mu = {th.N: Fr(1), th.Np: Fr(0)}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            mu[u] = Fr(th.n[P] - k, th.n[P])
    Lbar = [sum(mu[v] * Fr(lam[v][i]) for v in mu) for i in range(H)]
    return th.order_sign * sg(H * (H + 1) // 2 + 1) * c0 * Qform(th, lam, H, Lbar)


def sgn_htk(th, root):
    Bf = incidence_full(th.nv, th.edges)
    rows = [[Bf[e][v] for e in range(th.ne)] for v in range(th.nv) if v != root] + th.C3()
    return 1 if det(rows) > 0 else -1


def br(lam, vs):
    return det([[Fr(lam[v][i]) for v in vs] for i in range(len(vs))])


# ============================================================================ S1
print("=== S1 fan n=2 ===")
l = sp.symbols("l1:4")  # stand-ins: we integrate with symbolic brackets
r1, r2, r3, r4 = sp.symbols("r1:5")
b12, b13, b23 = sp.symbols("b12 b13 b23")
# antisymmetric pairing <x, y> on span(lambda_1, lambda_2, lambda_3), coefficients in the basis
def pair(x, y):
    B = {(0, 1): b12, (0, 2): b13, (1, 2): b23}
    s = 0
    for i in range(3):
        for j in range(3):
            if i < j:
                s += x[i] * y[j] * B[(i, j)]
            elif i > j:
                s -= x[i] * y[j] * B[(j, i)]
    return sp.expand(s)
L1, L2, L3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
y23 = [r3 * r1, r3, -r4]                               # eq-1-38, n = 2: r3 (r1 l1 + l2) - r4 l3
integrand = r2 * pair(L1, [L2[k] - y23[k] for k in range(3)]) * pair([r1, 1, 0], L3)   # eq-1-40, n = 2
I = sp.integrate(integrand, (r4, 0, 1 - r3))
I = sp.integrate(I, (r3, 0, 1))
I = sp.integrate(I, (r1, 0, 1), (r2, 0, 1))
I = sp.factor(sp.expand(I))
claim = -sp.Rational(1, 24) * pair(L1, [0, 2, 1]) * pair(L3, [1, 2, 0])
spider_rhs = sp.Rational(1, 24) * pair(L1, [0, 2, 1]) * pair(L3, [1, 2, 0])  # (v11, v21, A) = (1, 3, 2)
print("eq-1-40 (n=2) integrated:", I)
print("  equals claim -(1/24)[l1, l3+2 l2][l3, l1+2 l2]:", sp.expand(I - claim) == 0)
print("  at l2 = 0 equals eq-1-42 det(l1,l3)^2/24:", sp.expand(I.subs({b12: 0, b23: 0}) - b13 ** 2 / 24) == 0)
print("  equals minus RHS of eq-spider-two-arm-value with (v11,v21,A)=(1,3,2):", sp.expand(I + spider_rhs) == 0)
I_S1_poly = I

# pointwise: eq-operatope integrand on fan(2) (S1 vertex i = label i-1, root 3)
nv, edges = 4, [(0, 1), (0, 3), (1, 2), (1, 3), (2, 3)]
th = theta_from_graph(nv, edges, N=1, Np=3)
print("fan(2) as theta:", th.pqr, "paths", th.paths, "order_sign", th.order_sign)
trees = spanning_trees(nv, edges)
ok_tree, ratios = True, set()
for _ in range(3):
    t = rand_t(5, rng)
    lam = rand_sources(nv, 3, 2, rng)
    vals = [dyform_density(nv, edges, t, lam, 2, [e for e in range(5) if e not in T], sign_per_edge=-1) for T in trees]
    ok_tree &= all(v == vals[0] for v in vals)
    cl, _ = density_closed(th, t, lam, 2, 0)
    ratios.add(vals[0] / cl)
print("S1 integrand: tree-independent (%d trees): %s ; D_S1 / D_closed[C3] = %s" % (len(trees), ok_tree, ratios))
rS1 = ratios.pop()
okv = True
for _ in range(3):
    lam = rand_sources(nv, 3, 2, rng)
    IS1 = rS1 * K_C3_T0(th, lam, 2)
    sub = {b12: br(lam, [0, 1]), b13: br(lam, [0, 2]), b23: br(lam, [1, 2])}
    okv &= (sp.Rational(IS1.numerator, IS1.denominator) == I_S1_poly.subs(sub))
    Kh = sgn_htk(th, 3) ** 2 * K_C3_T0(th, lam, 2)
    okv &= (IS1 == -Kh)
print("S1: I^S1 (eq-operatope pointwise x own closed value) = eq-1-40 value, and I^S1 = -K_htk:", okv)

# ============================================================================ S2
print("=== S2 ===")
b = {k: sp.Symbol("b" + k) for k in ["123", "124", "125", "145", "245", "345", "uwq", "uwa", "uwb", "qab", "uab", "wab"]}
r = sp.symbols("s1:7")
JA = r[0] ** 2 * (r[3] * b["123"] + r[4] * b["124"] + r[5] * b["125"]) * \
     (b["345"] + (1 - r[1]) * b["145"] + (1 - r[2]) * b["245"])
IA = sp.integrate(JA, (r[5], 0, r[4]), (r[4], 0, r[3]), (r[3], 0, 1), (r[2], 0, r[1]), (r[1], 0, 1), (r[0], 0, 1))
claimA = sp.Rational(1, 432) * (3 * b["123"] + 2 * b["124"] + b["125"]) * (b["145"] + 2 * b["245"] + 3 * b["345"])
print("(A) from eq-chain + box with lambda_3 != 0: ", sp.factor(IA), " equals claim:", sp.expand(IA - claimA) == 0)
print("    at lambda_3 = 0 equals eq-IA:", sp.expand(IA.subs({b["123"]: 0, b["345"]: 0}) - sp.Rational(1, 3) * (b["124"] / 12 + b["125"] / 24) * (b["145"] / 6 + b["245"] / 3)) == 0)
JB = r[0] * (1 - r[0]) * (r[3] * b["uwq"] + r[4] * b["uwa"] + r[5] * b["uwb"]) * \
     (b["qab"] + (1 - r[1]) * b["uab"] + (1 - r[2]) * b["wab"])
IB = -sp.integrate(JB, (r[5], 0, r[4]), (r[4], 0, r[3]), (r[3], 0, 1), (r[2], 0, 1), (r[1], 0, 1), (r[0], 0, 1))
claimB = -sp.Rational(1, 288) * (3 * b["uwq"] + 2 * b["uwa"] + b["uwb"]) * (2 * b["qab"] + b["uab"] + b["wab"])
print("(B) from eq-bubble + box with lambda_q != 0:", sp.factor(IB), " equals claim:", sp.expand(IB - claimB) == 0)
print("    at lambda_q = 0 equals eq-IB:", sp.expand(IB.subs({b["uwq"]: 0, b["qab"]: 0}) + sp.Rational(1, 6) * (b["uwa"] / 12 + b["uwb"] / 24) * (b["uab"] / 2 + b["wab"] / 2)) == 0)


def S2_pointwise(nv, edges, N, Np, name, value_fn):
    th = theta_from_graph(nv, edges, N=N, Np=Np)
    Bf = incidence_full(nv, edges)
    trees = spanning_trees(nv, edges)
    ok_tree, ratios = True, set()
    for _ in range(2):
        t = rand_t(len(edges), rng)
        lam = rand_sources(nv, nv - 1, 3, rng)
        vals = []
        for T in trees:
            cot = [e for e in range(len(edges)) if e not in T]
            epsT = perm_sign(T + cot)
            BT = [[Bf[e][v] for e in T] for v in range(nv - 1)]
            vals.append(epsT * det(BT) * dyform_density(nv, edges, t, lam, 3, cot))
        ok_tree &= all(v == vals[0] for v in vals)
        cl, _ = density_closed(th, t, lam, 3, 0)
        ratios.add(vals[0] / cl)
    rr = ratios.pop() if len(ratios) == 1 else None
    okv = True
    for _ in range(3):
        lam = rand_sources(nv, nv - 1, 3, rng)
        IS2 = rr * K_C3_T0(th, lam, 3)
        okv &= (IS2 == value_fn(lam))
        okv &= (IS2 == -sgn_htk(th, nv - 1) ** 3 * K_C3_T0(th, lam, 3))
    print("S2 %s: theta%s, eq-kernel-def integrand tree-independent (%d trees): %s ; ratio to D_closed[C3]: %s ;"
          " I^S2 = general-source formula and = -K_htk: %s" % (name, th.pqr, len(trees), ok_tree, rr, okv))


# (A): vertices 1..5 -> 0..4, O = 5; edges (O1),(12),(23),(O3),(34),(45),(5O); N = S2 vertex 3 (label 2), N' = O
edgesA = [(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)]
fA = lambda lam: Fr(1, 432) * (3 * br(lam, [0, 1, 2]) + 2 * br(lam, [0, 1, 3]) + br(lam, [0, 1, 4])) * \
    (br(lam, [0, 3, 4]) + 2 * br(lam, [1, 3, 4]) + 3 * br(lam, [2, 3, 4]))
S2_pointwise(6, edgesA, 2, 5, "(A)", fA)
# (B): u, w, q, a, b -> 0..4, O = 5; N = q (label 2), N' = O
edgesB = [(5, 0), (0, 2), (5, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
fB = lambda lam: -Fr(1, 288) * (3 * br(lam, [0, 1, 2]) + 2 * br(lam, [0, 1, 3]) + br(lam, [0, 1, 4])) * \
    (2 * br(lam, [2, 3, 4]) + br(lam, [0, 3, 4]) + br(lam, [1, 3, 4]))
S2_pointwise(6, edgesB, 2, 5, "(B)", fB)

# ============================================================================ S4 ladder(2)
print("=== S4 ladder(2) ===")
a_ = lambda i: 2 * i
b_ = lambda i: 2 * i + 1
n = 2
edgesL, cot = [], []
for i in range(n - 1):
    cot.append(len(edgesL))
    edgesL += [(a_(i + 1), a_(i)), (a_(i), b_(i)), (b_(i), b_(i + 1))]
edgesL += [(a_(n - 1), b_(n - 1)), (b_(n - 1), b_(n)), (b_(n), a_(n)), (a_(n), a_(n - 1))]
cot.append(len(edgesL) - 3)
nvL = 2 * n + 2
thL = theta_from_graph(nvL, edgesL, N=2, Np=3)
print("edges", edgesL, "cotree", cot, "theta", thL.pqr, "interior", thL.interior)
ratios = set()
for _ in range(3):
    t = rand_t(len(edgesL), rng)
    lam = rand_sources(nvL, 5, 3, rng)
    ratios.add(dyform_density(nvL, edgesL, t, lam, 3, cot, root=5) / density_closed(thL, t, lam, 3, 0)[0])
print("D_S4 / D_closed[C3] over 3 points:", ratios)
rS4 = ratios.pop()
vals = []
for _ in range(3):
    while True:
        A, B, C = ([rng.randint(-4, 4) for _ in range(3)] for _ in range(3))
        Dl = det([A, B, C])
        if Dl != 0:
            break
    lam = [[Fr(0)] * 3 for _ in range(nvL)]
    lam[a_(0)], lam[b_(0)], lam[b_(n)] = [list(map(Fr, A)), list(map(Fr, B)), list(map(Fr, C))]
    lam[a_(n)] = [-(Fr(A[c]) + B[c] + C[c]) for c in range(3)]
    vals.append(rS4 * K_C3_T0(thL, lam, 3) / Dl ** 2)
print("I^S4 / Delta^2 at 3 corner-source sets:", [str(v) for v in vals], " (claim -1/432)")

# ============================================================================ S3
print("=== S3 ===")
okD = True
for pqr in [(1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5)]:
    th = theta_std(*pqr)
    lam = rand_sources(th.nv, th.Np, 5, rng)
    A_, B_, C_ = th.interior["A"], th.interior["B"], th.interior["C"]
    if pqr == (1, 6, 6):
        Dx = br(lam, B_) * br(lam, C_)
    elif pqr in ((2, 5, 6), (3, 4, 6)):
        Dx = br(lam, A_ + B_) * br(lam, C_)
    elif pqr == (3, 5, 5):
        Dx = sum(sg(d) * br(lam, A_ + [B_[k] for k in range(4) if k != d - 1]) * br(lam, [B_[d - 1]] + C_) for d in range(1, 5))
    else:
        Dx = sum(sg(d + 1) * br(lam, A_ + [B_[k] for k in range(3) if k != d - 1]) * br(lam, [B_[d - 1]] + C_) for d in range(1, 4))
    okD &= (Dx == D_theta(th, lam, 5))
print("S3 explicit D_166, D_256, D_346, D_355, D_445 equal D_theta:", okD)
okc = True
Rc = {(1, 6, 6): Fr(1, 158400), (2, 5, 6): Fr(1, 285120), (3, 4, 6): Fr(1, 332640), (3, 5, 5): Fr(1, 798336), (4, 4, 5): Fr(1, 997920)}
for pqr, c in Rc.items():
    fac = Fr(1)
    for n_ in pqr:
        fac /= math.factorial(n_ - 1)
    cc = fac * ghalf(6 - pqr[0]) * ghalf(6 - pqr[1]) * ghalf(6 - pqr[2]) / ghalf(6) / 2
    okc &= (cc == c)
    # K[C3] = (-1)^{15} c pi D_theta ; K^S3 = varsigma_theta K[C3] with varsigma_theta = -1 => + c pi D_theta (eq-Ktheta0)
    okc &= (-1 * sg(15) * cc == c)
print("S3 R-c values = Gamma formula, and varsigma_theta = -1 gives K^S3 = +c pi D_theta = -K[C3]:", okc)
print("time %.1f s" % (time.time() - t0))
