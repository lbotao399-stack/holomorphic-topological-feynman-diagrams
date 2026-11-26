"""psm.py -- Kontsevich weights (Poisson sigma model on the upper half plane) as T = 2 circuits with images.

Half plane H = {(x, y): y > 0}; ground vertices L = (0, 0) (f), R = (1, 0) (g); aerial vertices 0..n-1.
Edge e = (p -> h): Kontsevich propagator (1/2pi) d phi(p, h), phi(p, h) = arg(h - p) - arg(h - pbar)
   (Buring-Kiselev 2209.14438 Def. 1 / Kontsevich q-alg/9709040 sec 6.2; phi = 0 for p on the boundary).
Heat-kernel (GKW Higher.tex (3.15), T = 2) representation, s = Z / sqrt t:
   (1/2pi) d arg Z = int_0^oo psi(Z; t) dt,   psi = e^{-|Z|^2/t} (Z^1 dZ^2 - Z^2 dZ^1) / (2 pi t^2),
   so (1/2pi) d phi_e = int dt_e [psi(Z_e) - psi(Z_e^img)],  Z_e = h - p,  Z_e^img = h - pbar  (method of images).
Half plane -> plane: d phi(pbar, h) = -d phi(p, h), d phi(p, hbar) = d phi(p, h); every aerial vertex is the tail
   of two edges, so the integrand is even under y_v -> -y_v and  int_H w = int_{R^2} (sgn y_v / 2) w;
   sgn(y)/2 = int_0^oo chi(y; tau) dtau,  chi = y e^{-y^2/tau} / (2 sqrt(pi) tau^{3/2})   (T = 1 edge to the ground).
Hence  I_Gamma := (2 pi)^{-2n} int_{H^n} wedge_e d phi_e = int_{R_+^{3n}} d^{2n}t d^n tau  J(t, tau),
   J = int_{R^{2n}} prod_v chi_v  wedge_e [psi_e - psi_e^img]   (orientation dx_0 dy_0 dx_1 dy_1 ...),
   a sum over image patterns of Gaussian circuit integrals: parallel direction = ordinary incidence with
   Dirichlet data x_L = 0, x_R = 1; perpendicular direction = signed incidence (+1,+1 for image edges),
   ground y = 0, auxiliary edges v -> ground with times tau_v.
"""
import itertools, math
import numpy as np


def _pairings_moment(lin, mean, cov):
    """E[prod_k (lin_k . X)] for X ~ N(mean, cov); lin: list of (coef array (M, d), const (M,)); vectorized."""
    # shift: factor = lin_k.(X - mean) + (lin_k.mean + const)
    fac = []
    for a, c in lin:
        mu = np.einsum('mi,mi->m', a, mean) + c
        fac.append((a, mu))
    K = len(fac)

    def rec(idx):
        if not idx:
            return 1.0
        i0, rest = idx[0], idx[1:]
        a0, mu0 = fac[i0]
        tot = mu0 * rec(rest)
        for j in range(len(rest)):
            aj, muj = fac[rest[j]]
            cv = np.einsum('mi,mij,mj->m', a0, cov, aj)
            tot = tot + cv * rec(rest[:j] + rest[j + 1:])
        return tot
    return rec(tuple(range(K)))


def J_batch(n, edges, t, tau):
    """edges: list of (p, h), p in 0..n-1, h in 0..n-1 or 'L' or 'R'.  t: (M, 2n), tau: (M, n).  Returns J (M,)."""
    M = t.shape[0]
    E = len(edges)
    assert E == 2 * n
    total = np.zeros(M)
    for sigma in itertools.product((0, 1), repeat=E):
        # affine data: Z1_e = a_e . x + c_e ; Z2_e = b_e . y ; sign c_e
        a = np.zeros((E, n))
        c = np.zeros(E)
        b = np.zeros((E, n))
        sgn = 1.0
        for e, (p, h) in enumerate(edges):
            a[e, p] -= 1
            if h == 'R':
                c[e] += 1.0
            elif h != 'L':
                a[e, h] += 1
            if h not in ('L', 'R'):
                b[e, h] += 1
            if sigma[e] == 0:
                b[e, p] -= 1
            else:
                b[e, p] += 1
                sgn = -sgn
        # parallel Gaussian: sum_e (a_e.x + c_e)^2 / t_e
        w = 1.0 / t                                       # (M, E)
        Apar = np.einsum('me,ei,ej->mij', w, a, a)
        bpar = np.einsum('me,ei,e->mi', w, a, c)
        cpar = np.einsum('me,e->m', w, c * c)
        Ainv = np.linalg.inv(Apar)
        mean = -np.einsum('mij,mj->mi', Ainv, bpar)
        Epar = cpar - np.einsum('mi,mij,mj->m', bpar, Ainv, bpar)
        covpar = 0.5 * Ainv
        detpar = np.linalg.det(Apar)
        # perpendicular: sum_e (b_e.y)^2 / t_e + sum_v y_v^2 / tau_v
        Aperp = np.einsum('me,ei,ej->mij', w, b, b) + np.einsum('mv,vi,vj->mij', 1.0 / tau, np.eye(n), np.eye(n))
        Bp_inv = np.linalg.inv(Aperp)
        covperp = 0.5 * Bp_inv
        detperp = np.linalg.det(Aperp)
        zero_mean = np.zeros((M, n))
        pref = sgn * np.prod(1.0 / (2 * np.pi * t ** 2), axis=1) * np.prod(1.0 / (2 * np.sqrt(np.pi) * tau ** 1.5), axis=1)
        pref = pref * (np.pi ** (n / 2) / np.sqrt(detpar)) * np.exp(-Epar) * (np.pi ** (n / 2) / np.sqrt(detperp))
        acc = np.zeros(M)
        for A in itertools.combinations(range(E), n):
            # coefficient of dx_0 dy_0 dx_1 dy_1 ... in wedge_e theta_e, theta_e = Z1 dZ2 (e in A), -Z2 dZ1 (else)
            Mat = np.zeros((E, 2 * n))
            for e in range(E):
                if e in A:
                    Mat[e, 1::2] = b[e]
                else:
                    Mat[e, 0::2] = a[e]
            dcoef = np.linalg.det(Mat)
            if abs(dcoef) < 1e-12:
                continue
            linpar = [(np.broadcast_to(a[e], (M, n)), np.full(M, c[e])) for e in A]
            linperp = [(np.broadcast_to(b[e], (M, n)), np.zeros(M)) for e in range(E) if e not in A]
            linperp += [(np.broadcast_to(np.eye(n)[v], (M, n)), np.zeros(M)) for v in range(n)]
            mpar = _pairings_moment(linpar, mean, covpar)
            mperp = _pairings_moment(linperp, zero_mean, covperp)
            acc = acc + dcoef * ((-1) ** (E - n)) * mpar * mperp
        total = total + pref * acc
    return total


def weight_mc(n, edges, M, rng, scale=2.5, batch=200000):
    """importance-sampled Monte Carlo of I = int J over R_+^{3n} in log coordinates, logistic proposal."""
    d = 3 * n
    acc, acc2, cnt = 0.0, 0.0, 0
    while cnt < M:
        m = min(batch, M - cnt)
        u = rng.random((m, d))
        s = scale * np.log(u / (1 - u))                   # logistic(0, scale)
        logp = np.sum(-np.abs(s) / scale - 2 * np.log1p(np.exp(-np.abs(s) / scale)) - np.log(scale), axis=1)
        tt = np.exp(s)
        J = J_batch(n, edges, tt[:, :2 * n], tt[:, 2 * n:])
        val = J * np.prod(tt, axis=1) * np.exp(-logp)
        acc += val.sum()
        acc2 += (val ** 2).sum()
        cnt += m
    mean = acc / M
    err = math.sqrt(max(acc2 / M - mean ** 2, 0) / M)
    return mean, err


def weight_mc_radial(n, edges, M, rng, b=2.0, alpha=0.5, batch=100000):
    """I = int J over R_+^{3n};  proposal t = lambda u, log lambda ~ logistic(0, b), u ~ Dirichlet(alpha)."""
    d = 3 * n
    acc, acc2, cnt = 0.0, 0.0, 0
    lg = math.lgamma
    logB = d * lg(alpha) - lg(d * alpha)
    while cnt < M:
        m = min(batch, M - cnt)
        u01 = rng.random(m)
        sig = b * np.log(u01 / (1 - u01))
        logp_sig = -np.abs(sig) / b - 2 * np.log1p(np.exp(-np.abs(sig) / b)) - math.log(b)
        u = rng.dirichlet([alpha] * d, size=m)
        u = np.clip(u, 1e-300, None)
        logp_u = np.sum((alpha - 1) * np.log(u), axis=1) - logB
        lam = np.exp(sig)
        tt = lam[:, None] * u
        # p(t) = p_sig(log lam) * lam^{-d} * p_u(u)
        logp = logp_sig - d * sig + logp_u
        J = J_batch(n, edges, tt[:, :2 * n], tt[:, 2 * n:])
        val = J * np.exp(-logp)
        val = np.where(np.isfinite(val), val, 0.0)
        acc += val.sum()
        acc2 += (val ** 2).sum()
        cnt += m
    mean = acc / M
    err = math.sqrt(max(acc2 / M - mean ** 2, 0) / M)
    return mean, err
