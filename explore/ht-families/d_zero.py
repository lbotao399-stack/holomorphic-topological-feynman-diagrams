"""(d) Graph 2 = triangle with an ear on each side, degrees [2,2,2,4,4,4],
edges [(0,1),(0,2),(0,4),(0,5),(1,2),(1,3),(1,4),(3,4),(4,5)], (H,T) = (1,1), n = 4.

 D1  density g (small plan: three ears, then the central move) is not identically zero
 D2  integral of g over the central blocks only (w4, x4_2), with the ear ratios kept as variables,
     vanishes identically (symbolic sources); integral over the ear blocks only does not
 D3  g = c(ears) * d(sqrt Pi, l) / d(x4_2_2, w4_2) with c independent of the central coordinates,
     Pi = t_9 t_10 t_11 / T^3, l = current of the reduced edge 9 (exact, rational points)
 D4  first order in z: K_1 = int g (-sum_e y_e z_e) at integer sources and z  (nonzero?)
"""
import sys, random, time
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
from htk.simplex import dirichlet
import hf_brute as hb

edges = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)]
G = htk.Graph(6, edges)
red = spred.sp_reduce(G, prefer="small")
print(red.describe())
lsyms = sp.symbols("l0:5")
lam = [[s] for s in lsyms]
t0 = time.time()
dens = red.density(lam, T=1, check=3)
print("D1  exps:", dens.exponents())
print("    #terms of the polynomial part:", len(dens.num.terms()), " (%.1fs)" % (time.time() - t0))

# ---- D2: partial Dirichlet integration
ring = dens.ring
coords = red.coords
nc = len(coords)
blocks = red.blocks
central = [blk for blk in blocks if blk.move == 3]
ears = [blk for blk in blocks if blk.move != 3]
print("    central blocks:", [blk.bary for blk in central], " ear blocks:", [blk.bary for blk in ears])


def partial(dens, blks):
    """integrate over the given blocks; returns dict {monomial in remaining coords + params (tuple): (q, k)}
    as sympy expression"""
    cidx = {c: i for i, c in enumerate(coords)}
    bidx = {b: i for i, b in enumerate(red.bary)}
    intset = set(c for blk in blks for c in blk.indep)
    expr = sp.Integer(0)
    gens = [sp.Symbol(str(g)) for g in ring.gens]
    for mon, c in dens.num.terms():
        q, kpi = Fr(int(c.numerator), int(c.denominator)), 0
        for blk in blks:
            e0 = dens.exps[bidx[blk.bary[0]]]
            ex = (e0,) + tuple(dens.exps[bidx[b]] + mon[cidx[b]] for b in blk.indep)
            dq, dk = dirichlet(ex)
            q *= dq
            kpi += dk
        rest = sp.Integer(1)
        for i, g in enumerate(gens):
            if str(g) in intset:
                continue
            rest *= g ** mon[i]
        expr += sp.Rational(q.numerator, q.denominator) * sp.pi ** kpi * rest
    # monomial prefactor of the remaining blocks
    pref = sp.Integer(1)
    for b, e in zip(red.bary, dens.exps):
        if e != 0 and not any(b in blk.bary for blk in blks):
            pref *= sp.Symbol(b) ** sp.Rational(e.numerator, e.denominator)
    return sp.expand(expr), pref


t0 = time.time()
Fc, prefc = partial(dens, central)
print("D2  int over central blocks of g  = (%s) * [%s]   -> identically zero: %s  (%.1fs)" % (
    prefc, "polynomial with %d terms" % len(sp.Add.make_args(Fc)) if Fc != 0 else "0", Fc == 0, time.time() - t0))
Fe, prefe = partial(dens, ears)
print("    int over ear blocks of g      = (%s) * [polynomial with %d terms]  -> identically zero: %s" % (
    prefe, len(sp.Add.make_args(Fe)), Fe == 0))

# ---- D3: Jacobian structure, exact at rational points (numeric sources)
rng = random.Random(4)
lamn = [[Fr(rng.choice([-5, -3, -2, 2, 3, 7]))] for _ in range(5)]
densn = red.density(lamn, T=1, check=2)
mv4 = red.moves[3]
mv1 = red.moves[0]
# l = current of reduced edge 9 (= y_red of move 1) as polynomial in the coordinates
lpoly = sum((mv1.yred[v] * int(lamn[v][0]) for v in range(5)), red.ring.zero)
wb, xb = mv4.wblock, mv4.xblocks[1]
print("D3  central move:", mv4, " w-block", wb.bary, " x-block", xb.bary)
cid = {c: i for i, c in enumerate(coords)}


def evalp(p, pt):
    vals = tuple(Fr(pt[c]) for c in coords)
    return spred._evalp(p, vals)


def ratio_at(pt):
    pt = red.complete_point(pt)
    rat, rad = densn.evaluate(tuple(Fr(pt[c]) for c in coords))
    # sqrt Pi = w4_1 * sqrt(x4_2_1 x4_2_2 w4_2)  (t_9 = 1/w4_1, t_10 = x4_2_1/w4_2, t_11 = x4_2_2/w4_2)
    w1, w2 = pt[wb.bary[0]], pt[wb.bary[1]]
    x1, x2 = pt[xb.bary[0]], pt[xb.bary[1]]
    t9, t10, t11 = 1 / w1, x1 / w2, x2 / w2
    T = t9 + t10 + t11
    Pi = t9 * t10 * t11 / T ** 3
    # independent central coords: w2 (= w4_2) and x2 (= x4_2_2), with w1 = 1 - w2, x1 = 1 - x2
    # sqrt Pi = s * sqrt(r),  r = x1 x2 w2, s = w1:   d sqrt Pi = sqrt(r) [ds + s dr/(2r)]
    r = x1 * x2 * w2
    assert Pi == w1 ** 2 * r
    ds = {"w": Fr(-1), "x": Fr(0)}
    dr = {"w": x1 * x2, "x": (x1 - x2) * w2}
    dsq = {a: ds[a] + w1 * dr[a] / (2 * r) for a in ("w", "x")}       # times sqrt(r)
    lp = lpoly
    dl_w = evalp(lp.diff(red.ring.gens[cid[wb.bary[1]]]), pt)
    dl_x = evalp(lp.diff(red.ring.gens[cid[xb.bary[1]]]), pt)
    jac = dsq["x"] * dl_w - dsq["w"] * dl_x                              # d(sqrtPi, l)/d(x2, w2) / sqrt(r)
    # g = rat sqrt(rad);  ratio^2 = g^2 / (jac^2 r)
    return rat, rad, jac, r


for trial in range(3):
    base = red.random_point(rng)
    ear_pt = {c: base[c] for c in coords if not (c in wb.bary or c in xb.bary)}
    vals = []
    for _ in range(3):
        pc = dict(ear_pt)
        a, b = rng.randint(1, 20), rng.randint(1, 20)
        pc[wb.bary[1]] = Fr(a, a + rng.randint(1, 20))
        pc[xb.bary[1]] = Fr(b, b + rng.randint(1, 20))
        rat, rad, jac, r = ratio_at(pc)
        # c^2 = rat^2 rad / (jac^2 r);  sign(c) = sign(rat)*sign(jac)
        c2 = rat ** 2 * rad / (jac ** 2 * r)
        sg = (1 if rat > 0 else -1) * (1 if jac > 0 else -1)
        vals.append((sg, c2))
    print("    ear point %d: (sign, c^2) at three central points: %s  -> equal: %s" % (trial, vals, len(set(vals)) == 1))

# ---- D4: first order in z
z = [[Fr(rng.randint(-4, 4), rng.randint(1, 3))] for _ in range(G.ne)]
t0 = time.time()
Kz = red.taylor(lamn, z, 1, T=1, dens=densn)
print("D4  lambda =", [int(l[0]) for l in lamn], " z =", [str(x[0]) for x in z])
print("    K_0 = %s,  K_1 = %s  (%.1fs)" % (Kz[0], Kz[1], time.time() - t0))

# ---- D5: second order in z, and the moment with insertion y_0^2 y_3 (powers), at integer sources
t0 = time.time()
Kz2 = red.taylor(lamn, z, 2, T=1, dens=densn)
print("D5  K_2 (second order in z) = %s ;  moment <(y_0)^2 y_3> = %s  (%.1fs)" % (
    Kz2[2], red.integrate(lamn, T=1, powers={(0, 0): 2, (3, 0): 1}, dens=densn), time.time() - t0))

# ---- D6: every current is affine in l (current of reduced edge 9) at fixed ear coordinates:
#      d(y_e, l)/d(x4_2_2, w4_2) == 0 identically (symbolic sources)
cw = red.ring.gens[cid[wb.bary[1]]]
cx = red.ring.gens[cid[xb.bary[1]]]
ysym = red.currents(lam)
ring2 = ysym[0][0].ring
lsym = red._setup(red._flat(lam))
gens2 = {str(g): g for g in ring2.gens}
gw, gx = gens2[wb.bary[1]], gens2[xb.bary[1]]
l9 = None
# l from the response polynomials of move 1 (y_red of edge 9) with symbolic sources
lsy = ring2.zero
for v in range(5):
    lsy += mv1.yred[v].set_ring(ring2) * gens2["l%d" % v]
bad = 0
for e in range(G.ne):
    ye = ysym[e][0]
    J = ye.diff(gx) * lsy.diff(gw) - ye.diff(gw) * lsy.diff(gx)
    if J != 0:
        bad += 1
print("D6  #edges with d(y_e, l)/d(x4_2_2, w4_2) != 0 (symbolic):", bad, "of", G.ne)
print("    (nontriviality: dl/dw4_2 != 0: %s, dl/dx4_2_2 != 0: %s, #edges with dy_e/dw4_2 != 0: %d)" % (
    lsy.diff(gw) != 0, lsy.diff(gx) != 0, sum(1 for e in range(G.ne) if ysym[e][0].diff(gw) != 0)))
