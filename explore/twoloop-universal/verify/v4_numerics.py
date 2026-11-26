"""Direct numerical integration of the Schwinger integral (C12, and C3 at z != 0), own chart and own integrand.

Integrand: pullback of Omega = U^{-H} rho_2^T ^ prod psi (own float Kirchhoff solution of KCL + KVL with a BFS
cycle basis, own forms, cycle matrix C3) to the section R_C = 1 of R^N_{>0} -> RP^{N-1}_{>0}:
  K = int dxi  Omega(dt/dxi_1, ..., dt/dxi_{N-1}) sgn det[t; dt/dxi]       (no closed form, no explorer chart)
Chart xi = (v, th, u^A, u^B, u^C):  alpha_C = 1, rho = (1 - v^2)/v^2, c = sin^2(pi th / 2),
  alpha_A = rho c, alpha_B = rho (1 - c), R_P = 1/alpha_P, x^P = stick-breaking(u^P), t_{P_j} = R_P x^P_j.
Gauss-Legendre in every coordinate.  At z != 0 the factor exp(-sum_e y_e . z_e) uses the same currents.
Representation C3 at z != 0 is integrated separately on the same (v, th, u) grid mapped to (w, x).
"""
import sys, time, json, math, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.dont_write_bytecode = True
from vcore import theta_std, bfs_cycles, incidence_full, full_sources, D_theta, Qform, ghalf, sg, sigma_moments

DIR = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/"


def gl01(m):
    x, w = np.polynomial.legendre.leggauss(m)
    return (x + 1) / 2, w / 2


def stick(u):
    """u: (P, n-1) -> x: (P, n), dx: (P, n, n-1) = d x_j / d u_k."""
    P, k = u.shape
    n = k + 1
    x = np.empty((P, n))
    dx = np.zeros((P, n, k))
    rem = np.ones(P)
    drem = np.zeros((P, k))          # d rem / d u
    for j in range(k):
        x[:, j] = rem * u[:, j]
        dx[:, j, :] = drem * u[:, j:j + 1]
        dx[:, j, j] += rem
        newrem = rem * (1 - u[:, j])
        dnew = drem * (1 - u[:, j:j + 1])
        dnew[:, j] -= rem
        rem, drem = newrem, dnew
    x[:, n - 1] = rem
    dx[:, n - 1, :] = drem
    return x, dx


def integrand_batch(th, H, T, lam, z, v, tht, U_, C3, KCL, Cb):
    """v, tht: (P,), U_: dict P -> (P, n_P - 1) arrays of u's.  Returns Omega(J) sgn det[t;J] (P,)."""
    P = v.shape[0]
    ne = th.ne
    rho = (1 - v ** 2) / v ** 2
    drho = -2 / v ** 3
    c = np.sin(np.pi * tht / 2) ** 2
    dc = np.pi / 2 * np.sin(np.pi * tht)
    alpha = {"A": rho * c, "B": rho * (1 - c), "C": np.ones(P)}
    dal_v = {"A": drho * c, "B": drho * (1 - c), "C": np.zeros(P)}
    dal_t = {"A": rho * dc, "B": -rho * dc, "C": np.zeros(P)}
    R = {Pn: 1 / alpha[Pn] for Pn in "ABC"}
    t = np.empty((P, ne))
    J = np.zeros((P, ne, ne - 1))
    col = 2
    for Pn in "ABC":
        es = [e for e, _ in th.paths[Pn]]
        n = len(es)
        if n > 1:
            x, dx = stick(U_[Pn])
        else:
            x, dx = np.ones((P, 1)), np.zeros((P, 1, 0))
        for j, e in enumerate(es):
            t[:, e] = R[Pn] * x[:, j]
            J[:, e, 0] = -R[Pn] ** 2 * dal_v[Pn] * x[:, j]
            J[:, e, 1] = -R[Pn] ** 2 * dal_t[Pn] * x[:, j]
            J[:, e, col:col + n - 1] = R[Pn][:, None] * dx[:, j, :]
        col += n - 1
    # projective rescaling (Omega_{lam t}(lam J) = Omega_t(J))
    sc = t.max(axis=1)
    t = t / sc[:, None]
    J = J / sc[:, None, None]
    sgn = np.sign(np.linalg.det(np.concatenate([t[:, None, :], np.transpose(J, (0, 2, 1))], axis=1)))
    Ct = np.asarray(C3, float)
    M = np.einsum("ae,pe,be->pab", Ct, t, Ct)
    Ud = M[:, 0, 0] * M[:, 1, 1] - M[:, 0, 1] * M[:, 1, 0]
    expo = np.zeros(P)
    if H:
        A = np.empty((P, ne, ne))
        nk = KCL.shape[0]
        A[:, :nk, :] = KCL[None]
        A[:, nk:, :] = Cb[None] * t[:, None, :]
        rhs = np.zeros((P, ne, H))
        rhs[:, :nk, :] = lam[None]
        y = np.linalg.solve(A, rhs)                         # (P, ne, H)
        F = np.empty((P, 2 * H, ne))
        for i in range(H):
            for r in range(2):
                F[:, 2 * i + r, :] = Ct[r][None, :] * y[:, :, i]
        FJ = F @ J                                           # (P, 2H, ne-1)
        if z is not None:
            expo = -np.einsum("pei,ei->p", y, z)
    if T == 0:
        om = np.linalg.det(FJ) / Ud ** H
    else:
        Minv = np.empty_like(M)
        Minv[:, 0, 0], Minv[:, 1, 1] = M[:, 1, 1] / Ud, M[:, 0, 0] / Ud
        Minv[:, 0, 1], Minv[:, 1, 0] = -M[:, 0, 1] / Ud, -M[:, 1, 0] / Ud
        dMJ = np.einsum("ae,be,pel->pabl", Ct, Ct, J)       # dM_ab (J_l)
        m = ne - 1
        om = np.zeros(P)
        for k in range(m):
            for l in range(k + 1, m):
                a12 = np.zeros(P)
                for a in range(2):
                    for b in range(2):
                        a12 += Minv[:, a, b] * (dMJ[:, 0, a, k] * dMJ[:, b, 1, l] - dMJ[:, 0, a, l] * dMJ[:, b, 1, k])
                a12 *= -0.5
                rest = [x for x in range(m) if x not in (k, l)]
                d = np.linalg.det(FJ[:, :, rest]) if H else np.ones(P)
                om += (-1) ** (k + l + 3) * a12 * d
        om = om / Ud ** (H + 0.5)
    return om * sgn * np.exp(expo)


def kernel_numeric(th, H, T, lam_full, z=None, mv=40, mt=40, mu=4, batch=60000):
    ne = th.ne
    C3 = th.C3()
    Bf = incidence_full(th.nv, th.edges)
    KCL = np.array([[Bf[e][v] for e in range(ne)] for v in range(th.nv) if v != th.Np], float)
    Cb = np.array(bfs_cycles(th.nv, th.edges, th.Np), float)
    lam = np.array([[float(lam_full[v][i]) for i in range(H)] for v in range(th.nv) if v != th.Np]).reshape(-1, H) if H else None
    zz = None if z is None else np.array([[float(a) for a in ze] for ze in z])
    xv, wv = gl01(mv)
    xt, wt = gl01(mt)
    xu, wu = gl01(mu)
    nfib = {Pn: th.n[Pn] - 1 for Pn in "ABC"}
    nf = sum(nfib.values())
    fib_nodes = np.array(list(itertools.product(range(mu), repeat=nf))).reshape(-1, nf) if nf else np.zeros((1, 0), int)
    fib_u = xu[fib_nodes]
    fib_w = np.prod(wu[fib_nodes], axis=1) if nf else np.ones(1)
    tot = 0.0
    npts = 0
    VV, TT = np.meshgrid(np.arange(mv), np.arange(mt), indexing="ij")
    base = np.stack([VV.ravel(), TT.ravel()], 1)
    nb = len(fib_w)
    step = max(1, batch // nb)
    for s0 in range(0, len(base), step):
        bb = base[s0:s0 + step]
        P = len(bb) * nb
        v = np.repeat(xv[bb[:, 0]], nb)
        tht = np.repeat(xt[bb[:, 1]], nb)
        w = np.repeat(wv[bb[:, 0]] * wt[bb[:, 1]], nb) * np.tile(fib_w, len(bb))
        uall = np.tile(fib_u, (len(bb), 1))
        U_, k = {}, 0
        for Pn in "ABC":
            U_[Pn] = uall[:, k:k + nfib[Pn]]
            k += nfib[Pn]
        f = integrand_batch(th, H, T, lam, zz, v, tht, U_, C3, KCL, Cb)
        tot += float(np.sum(f * w))
        npts += P
    return tot, npts


def representation_numeric(th, H, T, lam_full, z, mv=40, mt=40, mu=4):
    """Claimed representation C3:  s (D_theta/2) int d^2w prod w^{H+1/2-n} int dx exp(sum lam.X)  (T = 1),
    s int d^2w prod w^{H-n} int dx Q(L(x)) exp(sum lam.X)  (T = 0), on the grid w_C = v^2,
    (w_A, w_B) = (1 - v^2)(c, 1 - c), c = sin^2(pi th/2), x = stick(u) (vectorized over fibre nodes)."""
    xv, wv = gl01(mv)
    xt, wt = gl01(mt)
    xu, wu = gl01(mu)
    lam = {vv: np.array([float(a) for a in lam_full[vv]]) for vv in range(th.nv)}
    zf = np.array([[float(a) for a in ze] for ze in z])
    nfib = {Pn: th.n[Pn] - 1 for Pn in "ABC"}
    nf = sum(nfib.values())
    fib = np.array(list(itertools.product(range(mu), repeat=nf))).reshape(-1, nf) if nf else np.zeros((1, 0), int)
    F = fib.shape[0]
    wf = np.prod(wu[fib], axis=1) if nf else np.ones(1)
    U = xu[fib]
    s = sg(H * (H + 1) // 2) if T == 1 else sg(H * (H + 1) // 2 + 1)
    Dth = float(D_theta(th, lam_full, H)) if T == 1 else None
    if T == 0:
        basis = np.eye(H)
        Bq = np.array([[float(Qbil(th, lam_full, H, basis[a], basis[b])) for b in range(H)] for a in range(H)])
    # fibre data: sigma^P_k (F,), tails of z
    sig, k = {}, 0
    for Pn in "ABC":
        n = th.n[Pn]
        if n > 1:
            u = U[:, k:k + n - 1]
            x, _ = stick(u)
            jx = np.ones(F)
            for kk in range(n - 2):
                jx *= (1 - u[:, kk]) ** (n - 2 - kk)
            wf = wf * jx
            sig[Pn] = [np.sum(x[:, kk:], axis=1) for kk in range(1, n)]
            k += n - 1
        else:
            sig[Pn] = []
    Zp = {Pn: sum(sgn_e * zf[e] for e, sgn_e in th.paths[Pn]) for Pn in "ABC"}
    tails = {Pn: [sum(sgn_e * zf[e] for e, sgn_e in th.paths[Pn][kk:]) for kk in range(1, th.n[Pn])] for Pn in "ABC"}
    Lfib = np.tile(lam[th.N], (F, 1)) if H else None
    for Pn in "ABC":
        for kk, vv in enumerate(th.interior[Pn], start=1):
            Lfib = Lfib + sig[Pn][kk - 1][:, None] * lam[vv][None, :]
    val = np.full(F, Dth / 2) if T == 1 else np.einsum("fa,ab,fb->f", Lfib, Bq, Lfib)
    tot = 0.0
    for iv in range(mv):
        v = xv[iv]
        for it in range(mt):
            c = math.sin(math.pi * xt[it] / 2) ** 2
            jac = 2 * v * (1 - v ** 2) * (math.pi / 2) * math.sin(math.pi * xt[it])
            w = {"A": (1 - v ** 2) * c, "B": (1 - v ** 2) * (1 - c), "C": v ** 2}
            W = sum(w[Pn] * Zp[Pn] for Pn in "ABC")
            mono = 1.0
            for Pn in "ABC":
                mono *= w[Pn] ** ((H + 0.5 - th.n[Pn]) if T == 1 else (H - th.n[Pn]))
            ex = np.full(F, -float(np.dot(lam[th.N], W)))
            for Pn in "ABC":
                for kk, vv in enumerate(th.interior[Pn], start=1):
                    ex += sig[Pn][kk - 1] * float(np.dot(lam[vv], Zp[Pn] - W)) - float(np.dot(lam[vv], tails[Pn][kk - 1]))
            tot += wv[iv] * wt[it] * jac * mono * float(np.sum(wf * val * np.exp(ex)))
    return s * tot


def Qbil(th, lam, H, a, b):
    from vcore import col, det
    zz = [0] * H
    aa = [Fr(float(x)).limit_denominator(1) for x in a]
    bb = [Fr(float(x)).limit_denominator(1) for x in b]
    cols = [col(aa, zz, H), col(zz, bb, H)] + [col(lam[w], zz, H) for w in th.interior["A"]] + \
           [col(lam[w], lam[w], H) for w in th.interior["B"]] + [col(zz, lam[w], H) for w in th.interior["C"]]
    return det([[cc[i] for cc in cols] for i in range(2 * H)])


def exact_K_C3(th, H, T, lam):
    """own closed z = 0 value of K[C3]: (rational, power of pi)."""
    fac = Fr(1)
    for n_ in th.pqr:
        fac /= math.factorial(n_ - 1)
    if T == 1:
        c1 = fac * ghalf(H + 1 - th.pqr[0]) * ghalf(H + 1 - th.pqr[1]) * ghalf(H + 1 - th.pqr[2]) / ghalf(H + 1) / 2
        return sg(H * (H + 1) // 2) * c1 * D_theta(th, lam, H), 1
    c0 = fac * Fr(math.prod(math.factorial(H - n_) for n_ in th.pqr), math.factorial(H + 1))
    mu = {th.N: Fr(1), th.Np: Fr(0)}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            mu[u] = Fr(th.n[P] - k, th.n[P])
    Lbar = [sum(mu[v] * Fr(lam[v][i]) for v in mu) for i in range(H)]
    return sg(H * (H + 1) // 2 + 1) * c0 * Qform(th, lam, H, Lbar), 0


if __name__ == "__main__":
    t0 = time.time()
    cases = json.load(open(DIR + "numerics_own.json"))
    ZONLY = len(sys.argv) > 1 and sys.argv[1] == "z"
    out = []
    print("=== C12: z = 0, explorer's sources (read from numerics_own.json as data) ===")
    for cs in ([] if ZONLY else cases):
        H, T, pqr = cs["H"], cs["T"], tuple(cs["pqr"])
        th = theta_std(*pqr)
        src = {v: [Fr(a) for a in cs["sources_nonroot"][v]] for v in range(th.nv - 1)}
        lam = full_sources(th.nv, th.Np, src)
        q, k = exact_K_C3(th, H, T, lam)
        ex = float(q) * math.pi ** k
        print("H=%d T=%d %s  own exact K[C3] = %s%s = %.15f   explorer K_C3_exact = %s" %
              (H, T, pqr, q, "*pi" if k else "", ex, cs["K_C3_exact"]), flush=True)
        for (mv, mt, mu) in [(12, 12, 3), (20, 20, 3), (30, 30, 4)]:
            ts = time.time()
            val, npts = kernel_numeric(th, H, T, lam, None, mv, mt, mu)
            dig = -math.log10(abs(val / ex - 1) + 1e-17)
            print("    m=(%d,%d,%d) points %d  numeric %.15f  digits %.1f  (%.1f s)" % (mv, mt, mu, npts, val, dig, time.time() - ts), flush=True)
            out.append({"H": H, "T": T, "pqr": pqr, "m": (mv, mt, mu), "value": val, "exact": ex, "digits": round(dig, 1)})
    print("=== C3 at z != 0: direct Schwinger integral vs claimed representation (same grids), and convergence ===")
    import random
    for H, T, pqr in [(1, 1, (1, 2, 2)), (2, 0, (1, 2, 2)), (2, 1, (2, 2, 3)), (3, 0, (2, 2, 3)), (3, 0, (1, 3, 3))]:
        th = theta_std(*pqr)
        rg = random.Random(31 + H)
        lam = full_sources(th.nv, th.Np, {v: [Fr(rg.randint(-6, 6), rg.randint(1, 3)) for _ in range(H)] for v in range(th.nv - 1)})
        z = [[Fr(rg.randint(-4, 4), 40) for _ in range(H)] for _ in range(th.ne)]
        q, k = exact_K_C3(th, H, T, lam)
        print("H=%d T=%d %s   K(lam;0) = %.12f" % (H, T, pqr, float(q) * math.pi ** k), flush=True)
        for (mv, mt, mu) in [(24, 24, 4), (30, 30, 5), (30, 30, 6), (36, 36, 7)]:
            ts = time.time()
            val, npts = kernel_numeric(th, H, T, lam, z, mv, mt, mu)
            rep = representation_numeric(th, H, T, lam, z, mv, mt, mu)
            print("    m=(%d,%d,%d) points %8d  direct %.15f  representation %.15f  digits(direct vs rep, same grid) %.1f  (%.1f s)"
                  % (mv, mt, mu, npts, val, rep, -math.log10(abs(val / rep - 1) + 1e-17), time.time() - ts), flush=True)
    if not ZONLY:
        json.dump(out, open(DIR + "verify/v4_numerics.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))
