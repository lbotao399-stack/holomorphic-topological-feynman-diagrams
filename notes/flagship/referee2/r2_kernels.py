# referee2: own computation of the (H,T)=(2,1) theta kernels and of two J_1-vertex kernels.
# No project imports.  Exact Fractions except where marked float.
import sys, random, itertools
from fractions import Fraction as Fr
from math import factorial
import mpmath as mp
import numpy as np
from scipy import integrate

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 20261121)

# ---------------- graphs, engine labelling (flagship Sec. 4.2) ----------------
G = {
 '133': dict(edges=[(0,5),(0,1),(1,2),(2,5),(0,3),(3,4),(4,5)],
             C=[[-1,1,1,1,0,0,0],[-1,0,0,0,1,1,1]],
             paths={'A':[0],'B':[1,2,3],'C':[4,5,6]}),
 '223': dict(edges=[(0,1),(1,5),(0,2),(2,5),(0,3),(3,4),(4,5)],
             C=[[1,1,-1,-1,0,0,0],[-1,-1,0,0,1,1,1]],
             paths={'A':[0,1],'B':[2,3],'C':[4,5,6]}),
}
NV, ROOT, NE = 6, 5, 7

def det(Mx):
    M = [[Fr(x) for x in r] for r in Mx]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c+1, n):
            f = M[r][c]/M[c][c]
            for k in range(c, n): M[r][k] -= f*M[c][k]
    return d

def solve(A, b):
    n = len(A); M = [[Fr(x) for x in A[i]] + [Fr(b[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]/M[c][c]
                for k in range(c, n+1): M[r][k] -= f*M[c][k]
    return [M[i][n]/M[i][i] for i in range(n)]

def Bmat(edges):  # B_{ev} = delta_{v,e0} - delta_{v,e1}, v != root
    return [[(1 if v == e[0] else 0) - (1 if v == e[1] else 0) for v in range(NV-1)] for e in edges]

def currents(g, t, lam):
    """lam: list of 2-vectors for v=0..4 (root gets -sum).  Kirchhoff: B^T y = lam, sum_e C_re t_e y_e = 0."""
    B = Bmat(g['edges']); C = g['C']
    A = [[B[e][v] for e in range(NE)] for v in range(NV-1)] + [[C[r][e]*t[e] for e in range(NE)] for r in range(2)]
    ys = []
    for i in range(2):
        ys.append(solve(A, [lam[v][i] for v in range(NV-1)] + [0, 0]))
    return [(ys[0][e], ys[1][e]) for e in range(NE)]   # y_e = (y^1_e, y^2_e)

# ---------------- exterior algebra on dt_0..dt_6 ----------------
def wedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            if set(ka) & set(kb): continue
            idx = list(ka) + list(kb)
            # sign of sorting permutation
            s = 1
            arr = idx[:]
            for i in range(len(arr)):
                for j in range(len(arr)-1-i):
                    if arr[j] > arr[j+1]: arr[j], arr[j+1] = arr[j+1], arr[j]; s = -s
            k = tuple(arr); out[k] = out.get(k, 0) + s*va*vb
    return {k: v for k, v in out.items() if v != 0}

def one(coeffs): return {(e,): c for e, c in enumerate(coeffs) if c != 0}
def scal(c, a): return {k: c*v for k, v in a.items()}
def add(*fs):
    out = {}
    for f in fs:
        for k, v in f.items(): out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}

def density_times_U72(g, t, lam):
    """returns (D(t) U^{7/2}, U, all components consistent?) with Omega = U^{-2} rho_2 psi^1_1 psi^1_2 psi^2_1 psi^2_2 = D iota_E vol."""
    C = g['C']
    M = [[sum(C[r][e]*C[s][e]*t[e] for e in range(NE)) for s in range(2)] for r in range(2)]
    U = M[0][0]*M[1][1] - M[0][1]*M[1][0]
    dM = [[one([C[r][e]*C[s][e] for e in range(NE)]) for s in range(2)] for r in range(2)]
    # rho_2 * U^{3/2} = -(M11 dM12 dM22 - M12 dM11 dM22 + M22 dM11 dM12)/2   (flagship eq. (1.3), checked below against Pf def)
    num = scal(Fr(-1, 2), add(scal(M[0][0], wedge(dM[0][1], dM[1][1])),
                               scal(-M[0][1], wedge(dM[0][0], dM[1][1])),
                               scal(M[1][1], wedge(dM[0][0], dM[0][1]))))
    # independent: Pf(-1/2 dM adj(M) dM)_{12}
    adj = [[M[1][1], -M[0][1]], [-M[1][0], M[0][0]]]
    A12 = {}
    for k in range(2):
        for l in range(2):
            A12 = add(A12, scal(Fr(-1, 2)*adj[k][l], wedge(dM[0][k], dM[l][1])))
    assert A12 == num, "rho_2 formula vs Pfaffian definition"
    y = currents(g, t, lam)
    psi = {(i, r): one([C[r][e]*y[e][i] for e in range(NE)]) for i in range(2) for r in range(2)}
    form = num
    for i in range(2):
        for r in range(2):
            form = wedge(form, psi[(i, r)])
    # form = U^{3/2} * U^2 * Omega ; Omega = D iota_E vol, iota_E vol = sum_j (-1)^j t_j dt_{0..^j..6}
    vals = []
    for j in range(NE):
        key = tuple(k for k in range(NE) if k != j)
        vals.append(form.get(key, 0) / ((-1)**j * t[j]))
    ok = all(v == vals[0] for v in vals)
    return vals[0], U, ok

def br(a, b): return a[0]*b[1] - a[1]*b[0]
def rnd(): return Fr(random.randint(-9, 9), random.randint(1, 5))
def rvec(): return (rnd(), rnd())

out = []
def P(*s):
    line = ' '.join(str(x) for x in s); print(line); out.append(line)

P("== 1. positive cycle basis: det[B^T; C] (engine C of flagship Sec. 4.2)")
for name, g in G.items():
    B = Bmat(g['edges'])
    Mx = [[B[e][v] for e in range(NE)] for v in range(NV-1)] + g['C']
    P("theta", name, "det =", det(Mx))

P("== 2. density at random t, random sources lambda_0..lambda_4 (lambda_N = lambda_0 arbitrary)")
for name, g in G.items():
    for trial in range(4):
        lam = [rvec() for _ in range(5)]
        t = [Fr(random.randint(1, 30), random.randint(1, 7)) for _ in range(NE)]
        DU, U, ok = density_times_U72(g, t, lam)
        p = lam[1:5]   # p_v = lambda_v, v = 1..4
        ratio = DU / (br(p[0], p[1]) * br(p[2], p[3]))
        RA = sum(t[e] for e in g['paths']['A']); RB = sum(t[e] for e in g['paths']['B']); RC = sum(t[e] for e in g['paths']['C'])
        P("theta", name, "trial", trial, "all 7 components equal:", ok, " U == RA RB + RB RC + RC RA:", U == RA*RB + RB*RC + RC*RA,
          " D U^{7/2}/([p1p2][p3p4]) =", ratio)

# ---------------- integrals ----------------
from math import gamma, pi
def dirichlet(ex):  # int_{Delta^2} prod w_P^{ex_P - 1} dw = prod Gamma(ex_P)/Gamma(sum ex_P), exact via mpmath
    q = lambda e: mp.mpf(Fr(e).numerator)/Fr(e).denominator
    return mp.fprod([mp.gamma(q(e)) for e in ex]) / mp.gamma(q(sum(ex)))

mp.mp.dps = 40
P("== 3. base integral  I = int_{Delta^6} U^{-7/2} dt  (path reduction + inversion R_P = 1/w_P -> Dirichlet)")
for name, g in G.items():
    n = [len(g['paths'][P_]) for P_ in 'ABC']
    I_exact = dirichlet([Fr(7, 2) - k for k in n]) / mp.fprod([factorial(k-1) for k in n])
    # independent float check in R-space: int_{Delta^2} prod R_P^{n_P-1}/(n_P-1)! U^{-7/2} dR, map R_A=u, R_B=(1-u)v, R_C=(1-u)(1-v)
    def f(u, v):
        RA = u; RB = (1-u)*v; RC = (1-u)*(1-v)
        Ux = RA*RB + RB*RC + RC*RA
        return (1-u) * RA**(n[0]-1) * RB**(n[1]-1) * RC**(n[2]-1) / Ux**mp.mpf(3.5) / mp.fprod([factorial(k-1) for k in n])
    mp.mp.dps = 25
    I_num = mp.quad(f, [0, 1], [0, 1])
    mp.mp.dps = 40
    c = Fr(-1, 2) if name == '133' else Fr(1, 2)
    P("theta", name, "n_P =", n, " I_exact =", mp.nstr(I_exact, 25), " = pi *", mp.nstr(I_exact/mp.pi, 20),
      " I_num(R-space, tanh-sinh) =", mp.nstr(I_num, 20), " rel.diff =", mp.nstr(abs(I_num/I_exact - 1), 3))
    P("   K/([p1p2][p3p4]) = (D U^{7/2} ratio) * I =", c, "* I =", mp.nstr(c*I_exact, 25), " ; claimed", "-pi/20" if name == '133' else "+pi/30",
      "=", mp.nstr(-mp.pi/20 if name == '133' else mp.pi/30, 25))

# ---------------- J_1-vertex kernels ----------------
# Physical sources: lambda_0 = 0, lambda_v = p_v (v=1..4), sum p = 0 => lambda_root = 0.
# Insertion [q_h1 q_h2], q at an external leg = p_v, at an edge tail (e,0) = -y_e, at an edge head (e,1) = +y_e (flagship Prop. 6.2).
def tpoint(g, w, xs):
    """t from w (R_P = 1/w_P) and x-fractions per path."""
    t = [None]*NE
    for P_, wv in zip('ABC', w):
        R = 1/wv
        for e, x in zip(g['paths'][P_], xs[P_]):
            t[e] = R*x
    return t

def qval(h, y, p):
    if h[0] == 'x': return p[h[1]-1]
    e, end = h
    return (-y[e][0], -y[e][1]) if end == 0 else (y[e][0], y[e][1])

# exact quadrature rules on simplices (barycentric nodes, weights summing to 1)
RULE2 = {1: [((Fr(1),), Fr(1))],
         2: [((Fr(1), Fr(0)), Fr(1, 6)), ((Fr(1, 2), Fr(1, 2)), Fr(2, 3)), ((Fr(0), Fr(1)), Fr(1, 6))],            # Simpson, deg 3
         3: [((Fr(1, 2), Fr(1, 2), Fr(0)), Fr(1, 3)), ((Fr(1, 2), Fr(0), Fr(1, 2)), Fr(1, 3)), ((Fr(0), Fr(1, 2), Fr(1, 2)), Fr(1, 3))]}  # edge midpoints, deg 2
RULE3 = {1: RULE2[1], 2: RULE2[2],
         3: [((Fr(1), Fr(0), Fr(0)), Fr(3, 60)), ((Fr(0), Fr(1), Fr(0)), Fr(3, 60)), ((Fr(0), Fr(0), Fr(1)), Fr(3, 60)),
             ((Fr(1, 2), Fr(1, 2), Fr(0)), Fr(8, 60)), ((Fr(1, 2), Fr(0), Fr(1, 2)), Fr(8, 60)), ((Fr(0), Fr(1, 2), Fr(1, 2)), Fr(8, 60)),
             ((Fr(1, 3), Fr(1, 3), Fr(1, 3)), Fr(27, 60))]}  # deg 3

def xavg(g, w, ins, rule, p):
    lam = [(Fr(0), Fr(0))] + [tuple(Fr(c) for c in pv) for pv in p]
    n = {P_: len(g['paths'][P_]) for P_ in 'ABC'}
    tot = Fr(0)
    for (xa, wa) in rule[n['A']]:
        for (xb, wb) in rule[n['B']]:
            for (xc, wc) in rule[n['C']]:
                t = tpoint(g, w, {'A': xa, 'B': xb, 'C': xc})
                y = currents(g, t, lam)
                tot += wa*wb*wc*br(qval(ins[0], y, p), qval(ins[1], y, p))
    return tot

def affine_test(g, w, p):
    lam = [(Fr(0), Fr(0))] + [tuple(Fr(c) for c in pv) for pv in p]
    def rx(n):
        a = [Fr(random.randint(1, 9)) for _ in range(n)]; s = sum(a); return [x/s for x in a]
    def yat(xs): return currents(g, tpoint(g, w, xs), lam)
    x0 = {P_: rx(len(g['paths'][P_])) for P_ in 'ABC'}
    x1 = {P_: rx(len(g['paths'][P_])) for P_ in 'ABC'}
    x2 = {P_: rx(len(g['paths'][P_])) for P_ in 'ABC'}
    # affine in x  <=>  y(x0 + s(x1-x0) + u(x2-x0)) - ... second differences vanish: check y(m) = (y(x1)+y(x2))/2 at midpoint, and 3-point collinearity
    xm = {P_: [(a+b)/2 for a, b in zip(x1[P_], x2[P_])] for P_ in 'ABC'}
    xt = {P_: [(a+2*b)/3 for a, b in zip(x0[P_], x1[P_])] for P_ in 'ABC'}
    y0, y1, y2, ym, yt = yat(x0), yat(x1), yat(x2), yat(xm), yat(xt)
    ok1 = all(ym[e][i] == (y1[e][i]+y2[e][i])/2 for e in range(NE) for i in range(2))
    ok2 = all(yt[e][i] == (y0[e][i]+2*y1[e][i])/3 for e in range(NE) for i in range(2))
    return ok1 and ok2

CASES = [  # (graph, vertex v, (h1, h2), claimed K^(v)/pi as function of b_ij)
 ('223', 3, (('x', 3), (4, 1)), lambda b: -b[1][2]*(b[1][3]+b[2][3])**2/84),
 ('133', 0, ((0, 0), (1, 0)),   lambda b:  b[1][2]*(b[1][3]+b[2][3])*(b[1][2]-2*b[1][3]-b[2][3])/252),
 ('133', 5, ((0, 1), (3, 1)),   lambda b: -b[1][2]*(b[1][3]+b[2][3])*(2*b[1][2]-b[1][3]-2*b[2][3])/252),
]

P("== 4. J_1-vertex kernels K^(v)(lambda;0) = int D(t) [q_h1 q_h2] dt, own Kirchhoff + exact quadrature + Dirichlet")
for (name, v, ins, claim) in CASES:
    g = G[name]
    n = [len(g['paths'][P_]) for P_ in 'ABC']
    p3 = [(random.randint(-7, 7), random.randint(-7, 7)) for _ in range(3)]
    p = p3 + [(-sum(q[0] for q in p3), -sum(q[1] for q in p3))]
    b = [[None]*5 for _ in range(5)]
    for i in range(1, 5):
        for j in range(1, 5): b[i][j] = Fr(br(p[i-1], p[j-1]))
    # density at the physical source configuration
    lam = [(Fr(0), Fr(0))] + [tuple(Fr(c) for c in pv) for pv in p]
    t = [Fr(random.randint(1, 30), random.randint(1, 7)) for _ in range(NE)]
    DU, U, ok = density_times_U72(g, t, lam)
    c = DU / (b[1][2]*b[3][4])
    # affine x-dependence of the currents at fixed R
    wr = [Fr(random.randint(1, 9)) for _ in range(3)]; s = sum(wr); wr = [x/s for x in wr]
    aff = affine_test(g, wr, p)
    # x-average: degree-2 rule vs degree-3 rule at random w
    agree = all(xavg(g, ww, ins, RULE2, p) == xavg(g, ww, ins, RULE3, p) for ww in
                [[x/sum(a) for x in a] for a in [[Fr(random.randint(1, 9)) for _ in range(3)] for _ in range(2)]])
    # fit g(w) = E_x[ins] as polynomial of degree <= 2 in (wA, wB), wC = 1 - wA - wB, verify at 4 extra points
    mons = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
    pts = []
    while len(pts) < 10:
        a = [Fr(random.randint(1, 13)) for _ in range(3)]; s = sum(a); pts.append([x/s for x in a])
    rows = [[ww[0]**i * ww[1]**j for (i, j) in mons] for ww in pts[:6]]
    rhs = [xavg(g, ww, ins, RULE2, p) for ww in pts[:6]]
    coef = solve(rows, rhs)
    fit_ok = all(sum(cf*ww[0]**i*ww[1]**j for cf, (i, j) in zip(coef, mons)) == xavg(g, ww, ins, RULE2, p) for ww in pts[6:])
    # exact integral: int_{Delta^2(w)} prod w_P^{5/2-n_P} wA^i wB^j dw ; wC exponent from expansion is already absorbed (basis in wA, wB only)
    val = mp.mpf(0)
    for cf, (i, j) in zip(coef, mons):
        ex = [Fr(7, 2) - n[0] + i, Fr(7, 2) - n[1] + j, Fr(7, 2) - n[2]]
        val += mp.mpf(cf.numerator)/cf.denominator * dirichlet(ex)
    val /= mp.fprod([factorial(k-1) for k in n])
    K = mp.mpf(c.numerator)/c.denominator * b[1][2]*b[3][4] * val
    Kc = claim(b) * mp.pi
    # independent float check: w-space quadrature with w_P = s_P^2 substitution and E_x from the raw rule (no fit), scipy dblquad
    def gw(wa, wb):
        wc = 1 - wa - wb
        lamf = [np.zeros(2)] + [np.array(pv, float) for pv in p]
        B = Bmat(g['edges']); Cm = g['C']
        tot = 0.0
        for (xa, wta) in RULE2[n[0]]:
            for (xb, wtb) in RULE2[n[1]]:
                for (xc, wtc) in RULE2[n[2]]:
                    tt = [0.0]*NE
                    for P_, wv, xx in zip('ABC', (wa, wb, wc), (xa, xb, xc)):
                        for e, x in zip(g['paths'][P_], xx): tt[e] = float(x)/wv
                    A = np.array([[B[e][vv] for e in range(NE)] for vv in range(NV-1)] + [[Cm[r][e]*tt[e] for e in range(NE)] for r in range(2)], float)
                    Y = np.linalg.solve(A, np.array([[lamf[vv][0], lamf[vv][1]] for vv in range(NV-1)] + [[0, 0], [0, 0]], float))
                    def q(h):
                        if h[0] == 'x': return np.array(p[h[1]-1], float)
                        return -Y[h[0]] if h[1] == 0 else Y[h[0]]
                    a1, a2 = q(ins[0]), q(ins[1])
                    tot += float(wta*wtb*wtc)*(a1[0]*a2[1]-a1[1]*a2[0])
        return tot
    def integrand(sb, sa):   # wA = sa^2 * (1 - ...) : use wA = sa^2, wB = (1-sa^2) sb^2 ... keep simple: triangle map
        wa = sa*sa; wb = (1-wa)*sb*sb; wc = 1 - wa - wb
        jac = 2*sa * (1-wa) * 2*sb
        return jac * wa**(2.5-n[0]) * wb**(2.5-n[1]) * wc**(2.5-n[2]) * gw(wa, wb)
    num, err = integrate.dblquad(integrand, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-12)
    num *= float(c) * float(b[1][2]*b[3][4]) / float(np.prod([factorial(k-1) for k in n]))
    P(f"theta({name[0]},{name[1]},{name[2]}) v={v} insertion [q_{ins[0]} q_{ins[1]}]  p =", p)
    P("   D U^{7/2}/([p1p2][p3p4]) at physical sources:", c, "(components consistent:", ok, ")")
    P("   currents affine in x at fixed R:", aff, "; degree-2 rule == degree-3 rule:", agree, "; degree-2 fit in w verified at 4 extra points:", fit_ok)
    P("   K^(v) exact (Dirichlet) =", mp.nstr(K, 30))
    P("   claimed (table, L')     =", mp.nstr(Kc, 30), "  -> ratio", mp.nstr(K/Kc, 25) if Kc != 0 else 'n/a')
    P("   float (scipy dblquad, w = s^2 maps, raw rule, no fit) =", repr(num), " rel.diff to exact:", float(abs(num/float(K) - 1)))

with open(__file__.replace('.py', '.log') if len(sys.argv) < 3 else sys.argv[2], 'w') as fh:
    fh.write('\n'.join(out) + '\n')
