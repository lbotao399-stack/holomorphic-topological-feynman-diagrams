"""VERIFY (b3), (c1): four-leg kernels, independent of htk and of the explorer's code.

1. H = 3, T = 0, exact (GM rules in own SP coordinates, own Kirchhoff solver):
   theta(1,3,3), theta(2,2,3) in engine labelling, sources p1..p4 at vertices 1..4, lambda_0 = 0,
   compared with the claimed polynomials at random integer momenta (no conservation imposed),
   and reduced mod conservation; box polygon(4) with sources p1,p2,p3 at vertices 0,1,2.
2. S2 comparison: graphs A and B in S2 labelling (sec7, eq-IA, eq-IB), I^{S2} = -K (n = 2).
3. H = 2, T = 1: theta(1,3,3), theta(2,2,3), sources at all non-root vertices (lambda_N included):
   pointwise exact identity  q |det J| = s (D_theta/2) prod_P w_P^{H - n_P}  (g = q U^{-1/2}|det J|,
   U = 1/(w_A w_B w_C)), then Dirichlet; closed form with mpmath; float numerical integral of the
   density computed from the definition at quadrature nodes (Gauss-Legendre after smoothing maps).
"""
import sys, time, random, itertools, math
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import vk, sp_coords as spc

rng = random.Random(4242)


def det3(a, b, c):
    return vk.det([[a[i], b[i], c[i]] for i in range(3)])


def br(p, i, j, k):  # 1-based brackets
    return det3(p[i - 1], p[j - 1], p[k - 1])


def rand_mom(k, H, lo=-7, hi=7):
    return [[rng.randint(lo, hi) for _ in range(H)] for _ in range(k)]


print("=== 1. H = 3, T = 0 four-leg kernels ===", flush=True)
H = 3
claims = {
    (1, 3, 3): lambda p: -Fr(1, 432) * (2 * br(p, 1, 2, 3) + br(p, 1, 2, 4)) * (2 * br(p, 1, 3, 4) + br(p, 2, 3, 4)),
    (2, 2, 3): lambda p: Fr(1, 288) * (2 * br(p, 1, 2, 3) + br(p, 1, 2, 4)) * (br(p, 1, 3, 4) + br(p, 2, 3, 4)),
}
for ns, f in claims.items():
    V, edges, interior = spc.paths_graph(ns)
    ext = [v for P in interior for v in P]
    assert ext == [1, 2, 3, 4]
    ok = True
    for trial in range(4):
        p = rand_mom(4, H)
        lam = {v: (p[ext.index(v)] if v in ext else [0, 0, 0]) for v in range(V - 1)}
        K = spc.integrate_T0(ns, lam, H, 1)
        K2 = spc.integrate_T0(ns, lam, H, 2) if trial == 0 else K
        ok &= (K == f(p)) and (K == K2)
        print("  theta%s p = %s: K = %s, claim %s, GM s=1 vs s=2 equal %s" % (ns, p, K, f(p), K == K2), flush=True)
    # conservation
    okc = True
    for trial in range(3):
        p = rand_mom(3, H)
        p.append([-(p[0][i] + p[1][i] + p[2][i]) for i in range(H)])
        lam = {v: (p[ext.index(v)] if v in ext else [0, 0, 0]) for v in range(V - 1)}
        K = spc.integrate_T0(ns, lam, H, 1)
        Dl = det3(p[0], p[1], p[2])
        target = -Fr(1, 432) * Dl ** 2 if ns == (1, 3, 3) else Fr(0)
        okc &= (K == target)
        print("  theta%s with p4 = -p1-p2-p3, Delta = %s: K = %s, claim %s" % (ns, Dl, K, target), flush=True)
    print("theta%s: K equals the claimed polynomial at 4 random points: %s; mod conservation: %s" % (ns, ok, okc))
    pD = [[3, 4, 3], [3, 4, -1], [-2, 4, 3], [-4, -12, -5]]
    lam = {v: (pD[ext.index(v)] if v in ext else [0, 0, 0]) for v in range(V - 1)}
    print("  explorer's numeric point (Delta = %s): K = %s" % (det3(*pD[:3]), spc.integrate_T0(ns, lam, H, 1)))

# box: polygon(4), sources p1, p2, p3 at vertices 0, 1, 2
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
from v_a_polynomiality import polygon_integral
okb = True
for trial in range(3):
    p = rand_mom(3, H)
    lam = {0: p[0], 1: p[1], 2: p[2]}
    K = polygon_integral(4, lam, 3, 2)
    okb &= (K == -det3(*p) / 6)
    print("  box p = %s: K = %s, -Delta/6 = %s" % (p, K, -det3(*p) / 6))
print("box (polygon(4), edges (k,k+1 mod 4), root 3): K = -Delta/6: %s" % okb)
pD = [[3, 4, 3], [3, 4, -1], [-2, 4, 3]]
print("  explorer's numeric point: K_box = %s" % polygon_integral(4, {0: pD[0], 1: pD[1], 2: pD[2]}, 3, 2))

print("\n=== 2. S2 graphs A, B (S2 labelling), I^S2 = -K ===", flush=True)
# S2 vertex v -> v-1, O = 5
GA = (6, [(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)], [[0, 1, 2], [3], [6, 5, 4]])
GB = (6, [(5, 0), (0, 2), (5, 1), (1, 2), (2, 3), (3, 4), (4, 5)], [[0, 1], [2, 3], [6, 5, 4]])
for name, (nv, edges, paths) in (("A", GA), ("B", GB)):
    ok = True
    for trial in range(3):
        if name == "A":
            l1, l2, l4, l5 = rand_mom(4, 3)
            lam = {0: l1, 1: l2, 2: [0, 0, 0], 3: l4, 4: l5}
            L = {1: l1, 2: l2, 4: l4, 5: l5}
            b = lambda i, j, k: det3(L[i], L[j], L[k])
            I_claim = Fr(1, 3) * (b(1, 2, 4) / 12 + b(1, 2, 5) / Fr(24)) * (b(1, 4, 5) / Fr(6) + b(2, 4, 5) / Fr(3))
        else:
            lu, lw, la, lb = rand_mom(4, 3)
            lam = {0: lu, 1: lw, 2: [0, 0, 0], 3: la, 4: lb}
            L = {'u': lu, 'w': lw, 'a': la, 'b': lb}
            b = lambda i, j, k: det3(L[i], L[j], L[k])
            I_claim = -Fr(1, 6) * (b('u', 'w', 'a') / Fr(12) + b('u', 'w', 'b') / Fr(24)) * (b('u', 'a', 'b') / Fr(2) + b('w', 'a', 'b') / Fr(2))
        K = spc.integrate_paths(nv, edges, paths, lam, 3, 1)
        ok &= (-K == I_claim)
        print("  graph %s: -K = %s, S2 eq-I%s = %s" % (name, -K, name, I_claim))
    print("graph %s: -K = S2 eq-I%s(lambda;0) at 3 random source sets: %s" % (name, name, ok), flush=True)

print("\n=== 3. H = 2, T = 1 thetas ===", flush=True)
H = 2


def Dtheta(ns, lam, interior):
    """det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]] (2H x 2H)."""
    LA = [lam[v] for v in interior[0]]
    LB = [lam[v] for v in interior[1]]
    LC = [lam[v] for v in interior[2]]
    cols = []
    for l in LA:
        cols.append(list(l) + [0] * H)
    for l in LB:
        cols.append(list(l) + list(l))
    for l in LC:
        cols.append([0] * H + list(l))
    return vk.det([[cols[j][i] for j in range(len(cols))] for i in range(2 * H)])


def gamma_closed(ns, H, s):
    val = mp.mpf(1)
    for n_ in ns:
        val *= mp.gamma(H + mp.mpf(3) / 2 - n_) / mp.factorial(n_ - 1)
    return s * val / (2 * mp.gamma(H + mp.mpf(3) / 2))


def dirichlet_half(a):
    """int_{Delta^2} prod w_P^{a_P} = prod Gamma(a_P+1) / Gamma(3 + sum a)."""
    return mp.fprod([mp.gamma(x + 1) for x in a]) / mp.gamma(3 + sum(a))


for ns, s_claim, Kclaim in (((1, 3, 3), -1, -mp.pi / 20), ((2, 2, 3), 1, mp.pi / 30)):
    V, edges, interior = spc.paths_graph(ns)
    C = vk.positive_cycles(V, edges, vk.dfs_tree(V, edges))
    ok = True
    for trial in range(5):
        Dt = 0
        while Dt == 0:
            lam = {v: [rng.randint(-6, 6) for _ in range(H)] for v in range(V - 1)}
            Dt = Dtheta(ns, lam, interior)
        w = [Fr(rng.randint(1, 20)) for _ in range(3)]
        sw = sum(w)
        w = [x / sw for x in w]
        xs = []
        for n_ in ns:
            xx = [Fr(rng.randint(1, 20)) for _ in range(n_)]
            sx = sum(xx)
            xs.append([x / sx for x in xx])
        t, J = spc.jac(list(ns), w, xs)
        q, U = vk.density(V, edges, t, lam, H, 1, C)
        dJ = vk.det(J)
        Uw = 1 / (w[0] * w[1] * w[2])
        target = s_claim * Dt / 2
        for P in range(3):
            target *= w[P] ** (H - ns[P])
        ok &= (U == Uw) and (q * abs(dJ) == target) and Dt != 0
    print("theta%s T=1: U = 1/(w_A w_B w_C) and q|det J| = s D_theta/2 prod w^(H-n_P) with s = %d at 5 random points/sources: %s" % (
        ns, s_claim, ok), flush=True)
    a = [H + mp.mpf(1) / 2 - n_ for n_ in ns]
    Kd = s_claim * dirichlet_half(a) / mp.fprod([mp.factorial(n_ - 1) for n_ in ns]) / 2
    print("  K/D_theta from Dirichlet = %s ; closed form = %s ; claim %s ; |diff| = %.1e" % (
        mp.nstr(Kd, 15), mp.nstr(gamma_closed(ns, H, s_claim), 15), mp.nstr(Kclaim, 15), float(abs(Kd - Kclaim))))


# ---------------- float density from the definition, vectorised, T = 1, n = 2
def batch_density_T1(nv, edges, C, t, lam, H):
    """t: (P, N) array. Returns D(t) (float) for each row."""
    P, N = t.shape
    g = nv - 1
    B = np.array(vk.incidence(nv, edges), dtype=float)  # N x g
    Cm = np.array(C, dtype=float)  # 2 x N
    L = np.einsum('ev,pe,ew->pvw', B, 1.0 / t, B)
    lamm = np.array([lam.get(v, [0] * H) for v in range(g)], dtype=float)  # g x H
    phi = np.linalg.solve(L, np.broadcast_to(lamm, (P, g, H)))
    y = np.einsum('ev,pvh->peh', B, phi) / t[:, :, None]  # P x N x H
    M = np.einsum('re,pe,se->prs', Cm, t, Cm)
    U = np.linalg.det(M)
    Mi = np.linalg.inv(M)
    # rho_ef = -1/2 sum_ab Mi_ab (u^a_e v^b_f - u^a_f v^b_e), u^a_e = C_0e C_ae, v^b_f = C_bf C_1f
    u = Cm[0][None, :] * Cm  # a x N
    v = Cm * Cm[1][None, :]  # b x N
    X = np.einsum('pab,ae,bf->pef', Mi, u, v)
    rho = -0.5 * (X - np.transpose(X, (0, 2, 1)))
    Psi = np.einsum('re,peh->phre', Cm, y).reshape(P, 2 * H, N)  # rows (i outer, r inner)
    e0 = 0
    S = [e for e in range(N) if e != e0]
    pos = {e: k for k, e in enumerate(S)}
    cf = np.zeros(P)
    for e, f in itertools.combinations(S, 2):
        rest = [x for x in S if x not in (e, f)]
        sgn = (-1) ** (pos[e] + pos[f] - 1)
        cf += sgn * rho[:, e, f] * np.linalg.det(Psi[:, :, rest])
    return cf / t[:, e0] / U ** H / np.sqrt(U)


def numeric_theta_T1(ns, lam, m):
    """int over Delta^2 x prod Delta^{n_P-1} of D(t(xi)) |det J| by Gauss-Legendre after maps:
    w: w_A = 1-u, (w_B, w_C) = u (sin^2 phi, cos^2 phi), u = 1-(1-s)^2 ... (smooth maps);
    x^P (n_P = 2): x_2 = sin^2;  (n_P = 3): Duffy x = (1-a, a b, a (1-b))."""
    V, edges, interior = spc.paths_graph(ns)
    C = vk.positive_cycles(V, edges, vk.dfs_tree(V, edges))
    gx, gw = np.polynomial.legendre.leggauss(m)
    gx = 0.5 * (gx + 1)
    gw = 0.5 * gw
    # 1-dim nodes on (0,1) for each coordinate
    axes = []
    # w block: u in (0,1) with u = 1 - (1-s)^4 (s in (0,1)), phi = pi/2 * r
    su = gx
    u = 1 - (1 - su) ** 4
    du = 4 * (1 - su) ** 3
    r = gx
    sph = np.sin(np.pi / 2 * r) ** 2
    dsph = np.pi / 2 * 2 * np.sin(np.pi / 2 * r) * np.cos(np.pi / 2 * r)
    grids = [(u, du * gw), (sph, dsph * gw)]
    for n_ in ns:
        if n_ == 2:
            grids.append((gx, gw))
        elif n_ == 3:
            grids.append((gx, gw))
            grids.append((gx, gw))
    nd = len(grids)
    idx = np.array(list(itertools.product(range(m), repeat=nd)))
    total = 0.0
    chunk = 200000
    for c0 in range(0, len(idx), chunk):
        I = idx[c0:c0 + chunk]
        vals = [grids[d][0][I[:, d]] for d in range(nd)]
        wts = np.prod([grids[d][1][I[:, d]] for d in range(nd)], axis=0)
        uu, vv = vals[0], vals[1]
        wA = 1 - uu
        wB = uu * vv
        wC = uu * (1 - vv)
        jw = uu  # Duffy Jacobian of (w_B, w_C) in (u, v)
        k = 2
        xs = []
        jx = np.ones(len(I))
        for n_ in ns:
            if n_ == 1:
                xs.append([np.ones(len(I))])
            elif n_ == 2:
                a = vals[k]; k += 1
                xs.append([1 - a, a])
            else:
                a = vals[k]; b = vals[k + 1]; k += 2
                xs.append([1 - a, a * b, a * (1 - b)])
                jx = jx * a
        wv = [wA, wB, wC]
        t = np.stack([xs[P][j] / wv[P] for P in range(3) for j in range(ns[P])], axis=1)
        # |det[t; d t / d xi]| for xi = (w_B, w_C, x^P_2..) : build numerically per row
        Nn = t.shape[1]
        J = np.zeros((len(I), Nn, Nn))
        J[:, 0, :] = t
        row = 1
        col = 0
        cols = {}
        for P in range(3):
            for j in range(ns[P]):
                cols[(P, j)] = col
                col += 1
        for Q in (1, 2):
            for P in range(3):
                dP = (1 if P == Q else 0) - (1 if P == 0 else 0)
                for j in range(ns[P]):
                    J[:, row, cols[(P, j)]] = -xs[P][j] / wv[P] ** 2 * dP
            row += 1
        for P in range(3):
            for kk in range(1, ns[P]):
                for j in range(ns[P]):
                    J[:, row, cols[(P, j)]] = ((1 if j == kk else 0) - (1 if j == 0 else 0)) / wv[P]
                row += 1
        dJ = np.abs(np.linalg.det(J))
        Dv = batch_density_T1(V, edges, C, t, lam, H)
        total += np.sum(Dv * dJ * jw * jx * wts)
    return total


for ns, Kc in (((1, 3, 3), -mp.pi / 20), ((2, 2, 3), mp.pi / 30)):
    V, edges, interior = spc.paths_graph(ns)
    Dt = 0
    while Dt == 0:
        lam = {v: [rng.randint(-5, 5) for _ in range(H)] for v in range(V - 1)}
        Dt = Dtheta(ns, lam, interior)
    exact = float(Kc * Dt)
    for m in (8, 12):
        t0 = time.time()
        val = numeric_theta_T1(ns, lam, m)
        print("theta%s T=1 sources %s: D_theta = %s, exact K = %.15f, Gauss-Legendre m=%d: %.15f, rel.dev %.1e (%.0fs)" % (
            ns, [lam[v] for v in range(V - 1)], Dt, exact, m, val, abs(val - exact) / abs(exact), time.time() - t0), flush=True)
