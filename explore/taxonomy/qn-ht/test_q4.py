"""(Q^2)_{hbar^4} = 0 at H = 1 on undifferentiated monomials (z = 0, kappa = 0, cubic W):
  Q2^hat Q2^hat + u {Q0, Q4^fan(4)} + v {Q0, Q4^tadpole(5,0)} = 0,
fan(4): K = (pi^2/180)(l0 l4)^2, tadpole(5,0): K = (pi^2/630) l0^4 (3dN2 REPORT Sec. 3.2-3.3, converted to
the (-v,w) edge order by sgn(pi)^{D+1}); expected u = v = 1 (N_4 / N_2^2 = +1).  Q2^hat Q2^hat uses the
derivative moments of the thetas.  Usage: python3 test_q4.py <theory> <maxdeg> <nsample> [B4X]  Log: test_q4_<theory>.log"""
import sys, itertools, random, time
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng, theories as T, specs

TH = {
    "u1": lambda: T.u1_3d([1]),
    "su2adj": lambda: hteng.Theory(1, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))]),
    "gl2f": None,
    "sqed1": lambda: T.sqed_3d(1),
    "sqed2": lambda: T.sqed_3d(2),
    "su2N4_1": lambda: T.su2_N4_3d(1),
    "su2N8": lambda: specs.su2_N8_sl2(),
    "su2N8so3": lambda: T.su2_N8_3d(),
}

name, maxdeg, nsample = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
if name == "gl2f":
    import ref_theories as RT
    thr = RT.gl_theory(2, 1, 0)
    th = hteng.Theory(1, thr.G, thr.f, bf=[(m[0], m[1]) for m in thr.multiplets], name="GL(2)+1F")
else:
    th = TH[name]()
log = open("test_q4_%s%s.log" % (name, "_" + sys.argv[4] if len(sys.argv) > 4 else ""), "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()
t0 = time.time()
E2 = hteng.Engine(th)
E4 = hteng.Engine(th, shapes=hteng.h1_q4_shapes())
pool = th.letter_pool()
monos = []
mode = sys.argv[4] if len(sys.argv) > 4 else "all"
if mode == "B4X":
    Bs = [L for L in pool if L[0] == 'B']
    for b4 in itertools.combinations_with_replacement(Bs, 4):
        for X in pool:
            O = th.mono(*(list(b4) + [X]))
            if O and O not in monos:
                monos.append(O)
else:
    for deg in range(4, maxdeg + 1):
        for combo in itertools.combinations_with_replacement(pool, deg):
            O = th.mono(*combo)
            if O:
                monos.append(O)
rng = random.Random(11)
sample = monos if len(monos) <= nsample else rng.sample(monos, nsample)
eqs = []
nX = nYf = nYt = 0
for O in sample:
    X = E2.Q2(E2.Q2(O))
    Yf = hteng.psum(E2.Q0(E4.Q2(O, [0])), E4.Q2(E2.Q0(O), [0]))
    Yt = hteng.psum(E2.Q0(E4.Q2(O, [1])), E4.Q2(E2.Q0(O), [1]))
    nX += bool(X); nYf += bool(Yf); nYt += bool(Yt)
    for m in set(X) | set(Yf) | set(Yt):
        eqs.append((X.get(m, Fr(0)), Yf.get(m, Fr(0)), Yt.get(m, Fr(0))))
out("theory %s: %d undifferentiated monomials of degree 4..%d (sample of %d); Q2^2 != 0 on %d, "
    "{Q0,Q4fan} != 0 on %d, {Q0,T50} != 0 on %d; %d coefficient equations  [%.0fs]"
    % (th.name, len(sample), maxdeg, len(monos), nX, nYf, nYt, len(eqs), time.time() - t0))
# exact solve of X + u Yf + v Yt = 0 (least: Gaussian elimination on [Yf Yt | -X])
rows = [(yf, yt, -x) for x, yf, yt in eqs if (x, yf, yt) != (0, 0, 0)]
# rank analysis
def solve(rows):
    piv = []
    R = [list(r) for r in rows]
    used = []
    r = 0
    for col in range(2):
        p = next((i for i in range(r, len(R)) if R[i][col] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        pv = R[r][col]
        R[r] = [x / pv for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][col] != 0:
                f = R[i][col]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        used.append(col)
        r += 1
    incons = sum(1 for i in range(r, len(R)) if R[i][2] != 0)
    sol = {used[i]: R[i][2] for i in range(r)}
    return sol, incons, used
sol, incons, used = solve(rows)
out("solution: u = %s, v = %s; inconsistent equations: %d (unknowns determined: %s)"
    % (sol.get(0, "undetermined"), sol.get(1, "undetermined"), incons, ["u", "v"] if used == [0, 1] else used))
out("[%.0fs]" % (time.time() - t0))
