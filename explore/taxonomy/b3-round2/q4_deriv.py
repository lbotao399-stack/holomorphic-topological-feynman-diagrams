"""(Q^2)_{hbar^4} = 0 at H = 1 on monomials with derivative letters (kappa = 0, cubic theories):
  Q2^ Q2^ + {Q0, Q4^fan(4)} + {Q0, T^_{5,0}} = 0     (u = v = 1 fixed by the z = 0 test, qn-ht §1)
fan(4): exact htk Taylor moments (fan4_taylor.pkl, edge order (-v,w)).
tadpole(5,0) (not series-parallel): z = 0 moment (1/630) l0^4 [src 3dN2 §3.3]; order-m moments
M_p = t_p l0^{4+m} (|p| = m <= dmax) are unknowns, solved exactly from the identity.
Usage: python3 q4_deriv.py <theory> <dmax> <nsample>   theory in {u1, gl2f, su2adj}   Log: q4_deriv_<theory>_<dmax>.log"""
import sys, os, itertools, random, time, pickle
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
QN = os.path.join(HERE, "..", "qn-ht")
sys.path.insert(0, QN)
import sympy as sp
import hteng, theories as T
from hteng import Shape, perm_sign, psum

name, dmax, nsample = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
ctrl = len(sys.argv) > 4 and sys.argv[4] == "ctrl"  # negative control: order >= 1 fan(4) moments doubled
log = open(os.path.join(HERE, "q4_deriv_%s_%d%s.log" % (name, dmax, "_ctrl" if ctrl else "")), "w")
def out(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()
if name == "gl2f":
    import ref_theories as RT
    thr = RT.gl_theory(2, 1, 0)
    th = hteng.Theory(1, thr.G, thr.f, bf=[(m[0], m[1]) for m in thr.multiplets], name="GL(2)+1F")
elif name == "u1":
    th = T.u1_3d([1])
elif name == "su2adj":
    th = hteng.Theory(1, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))], name="SU(2)+adj")
t0 = time.time()
# ---------------- fan(4) moments from htk Taylor
R = 5
fe = []
for i in range(4):
    fe += [(i, i + 1), (i, R)]
fe.append((4, R))
mine = sorted(fe, key=lambda e: (-e[0], e[1]))
dat = pickle.load(open(os.path.join(HERE, "fan4_taylor.pkl"), "rb"))
assert dat["edges"] == mine
Mtay = len(dat["K"]) - 1
assert perm_sign(fe, mine) ** 3 == 1  # h1_q4_shapes sign = +1, consistent with K_0 = +a^2 b^2/180 here
a, b = sp.symbols("a b")
zsym = [sp.Symbol("z%d" % v) for v in range(5)]
# root edges processed in hteng.F order: v = 4, 3, 2, 1, 0
rootorder_fan = [4, 3, 2, 1, 0]
fanmom = {}
for m in range(Mtay + 1):
    P = sp.Poly(sp.expand(dat["K"][m] / sp.pi ** 2), *(zsym + [a, b]))
    for mon, c in P.terms():
        zp, (ia, ib) = mon[:5], mon[5:]
        sgn = (-1) ** sum(zp)
        fac = 1
        for p in zp:
            fac *= sp.factorial(p)
        key = tuple((zp[v],) for v in rootorder_fan)
        val = Fr(int(sp.numer(c * sgn * fac)), int(sp.denom(c * sgn * fac)))
        if ctrl and m >= 1:
            val *= 2
        fanmom.setdefault(key, {})[((ia,), (ib,))] = val
def mfan(key):
    if sum(p[0] for p in key) > Mtay:
        raise ValueError("fan(4) Taylor order too low")
    return fanmom.get(tuple(key), {})
sh_fan = Shape(6, mine, 2, [0, 4], mfan, "fan(4)")
# ---------------- tadpole(5,0)
te = sorted([(min(x, y), max(x, y)) for x, y in [(i, (i + 1) % 5) for i in range(5)]] + [(v, R) for v in range(1, 5)])
mt = sorted(te, key=lambda e: (-e[0], e[1]))
sg2 = perm_sign(te, mt) ** 3
out("fan(4) Taylor order %d available%s" % (Mtay, "; NEGATIVE CONTROL: order>=1 fan moments x2" if ctrl else ""))
out("tadpole(5,0) edges (order -v,w):", mt, "; z=0 moment sign", sg2)
rootorder_tad = [v for v in range(4, -1, -1) for (x, w) in mt if x == v and w == R]
unk = []
for m in range(1, dmax + 1):
    for p in itertools.product(range(m + 1), repeat=4):
        if sum(p) == m:
            unk.append(tuple((x,) for x in p))
out("tadpole root edges in key order: (%s,5); unknown moments t_p: %d" % (rootorder_tad, len(unk)))
def mtad0(key):
    if any(p[0] for p in key):
        return {}
    return {((4,),): Fr(1, 630) * sg2}
def mtad_unk(j):
    def f(key, j=j):
        if tuple(key) == unk[j]:
            return {((4 + sum(p[0] for p in key),),): Fr(1)}
        return {}
    return f
E2 = hteng.Engine(th)
Ef = hteng.Engine(th, shapes=[sh_fan])
Et0 = hteng.Engine(th, shapes=[Shape(6, mt, 2, [0], mtad0, "tad0")])
Etu = [hteng.Engine(th, shapes=[Shape(6, mt, 2, [0], mtad_unk(j), "tad%d" % j)]) for j in range(len(unk))]
def anti(Eng, O):
    return psum(E2.Q0(Eng.Q2(O)), Eng.Q2(E2.Q0(O)))
# ---------------- monomials: degree 4..5 (cubic: Q4 needs >= 4 letters), total derivative 1..dmax
pool0 = th.letter_pool()
monos = []
for deg in (4, 5):
    for combo in itertools.combinations_with_replacement(pool0, deg):
        af = sum(1 for L in combo if L[0] == 'B' or L[0][0] == 's')
        if af < 4:
            continue
        for dist in itertools.product(range(dmax + 1), repeat=deg):
            if not 1 <= sum(dist) <= dmax:
                continue
            lets = [(L[0], L[1], (d,)) for L, d in zip(combo, dist)]
            O = th.mono(*lets)
            if O and O not in monos:
                monos.append(O)
rng = random.Random(7)
sample = monos if len(monos) <= nsample else rng.sample(monos, nsample)
out("theory %s: %d monomials (deg 4-5, N_B + N_psi >= 4, 1 <= total derivative <= %d), sample %d" % (th.name, len(monos), dmax, len(sample)))
rows = []
cnt = dict(X=0, F=0, T0=0, TU=0)
for O in sample:
    X = E2.Q2(E2.Q2(O))
    Yf = anti(Ef, O)
    Y0 = anti(Et0, O)
    Yu = [anti(E, O) for E in Etu]
    cnt['X'] += bool(X); cnt['F'] += bool(Yf); cnt['T0'] += bool(Y0); cnt['TU'] += any(bool(y) for y in Yu)
    keys = set(X) | set(Yf) | set(Y0)
    for y in Yu:
        keys |= set(y)
    for m in keys:
        rows.append([y.get(m, Fr(0)) for y in Yu] + [-(X.get(m, Fr(0)) + Yf.get(m, Fr(0)) + Y0.get(m, Fr(0)))])
out("nonzero on: Q2^2 %d, {Q0,Q4fan} %d, {Q0,T50 z=0} %d, {Q0,T50 unknown moments} %d ; equations %d  [%.0fs]"
    % (cnt['X'], cnt['F'], cnt['T0'], cnt['TU'], len(rows), time.time() - t0))
# exact elimination
nu = len(unk)
Rm = [r[:] for r in rows if any(r)]
piv = []
r = 0
for col in range(nu):
    p = next((i for i in range(r, len(Rm)) if Rm[i][col] != 0), None)
    if p is None:
        continue
    Rm[r], Rm[p] = Rm[p], Rm[r]
    pv = Rm[r][col]
    Rm[r] = [x / pv for x in Rm[r]]
    for i in range(len(Rm)):
        if i != r and Rm[i][col] != 0:
            f = Rm[i][col]
            Rm[i] = [x - f * y for x, y in zip(Rm[i], Rm[r])]
    piv.append(col)
    r += 1
incons = sum(1 for i in range(r, len(Rm)) if Rm[i][nu] != 0)
out("rank %d of %d unknowns; inconsistent equations: %d" % (r, nu, incons))
for i, col in enumerate(piv):
    lhs = " + ".join("%s*t%s" % (Rm[i][c], tuple(p[0] for p in unk[c])) for c in range(nu) if Rm[i][c] != 0)
    out("  %s = %s" % (lhs, Rm[i][nu]))
und = [tuple(p[0] for p in unk[c]) for c in range(nu) if c not in piv]
out("undetermined: %s" % und)
out("[%.0fs]" % (time.time() - t0))
