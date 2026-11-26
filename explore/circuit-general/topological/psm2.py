"""psm2.py -- Kontsevich weights from the circuit representation with images, deterministic quadrature.

Same representation as psm.py; the auxiliary T = 1 edges are integrated first,
   int_0^oo chi(y; tau) dtau = sgn(y)/2,
so for each image pattern the perpendicular Gaussian carries prod_v sgn(y_v)/2 (closed forms, n <= 2:
   E|y| = sqrt(2 S/pi);  E[s] = (2/pi) arcsin rho,  E[y_i y_j s] = S_ij E[s] + (Stein terms) for s = sgn y0 sgn y1),
and the parallel Gaussian has Dirichlet data x_L = 0, x_R = 1:  e^{-E/lambda} with E = 1/R_LR(u).
Homogeneity in the 2n edge times t = lambda u:  J(lambda u) = lambda^{-a} e^{-E(u)/lambda} (P0(u) + lambda P1(u)),
P0 from the products of parallel means, P1 from the parallel covariance;  a = 5/2 (n = 1), 5 (n = 2).
Claim: P1 = 0 (checked pointwise); then I = Gamma(a - 2n) int_{Delta^{2n-1}} P0(u) E(u)^{2n-a} du.
"""
import itertools, math
import numpy as np


def _inv_det(A):
    """explicit inverse and determinant for (M, k, k), k = 1, 2 (no exception on near-singular input)"""
    k = A.shape[1]
    if k == 1:
        d = A[:, 0, 0]
        return (1.0 / d)[:, None, None], d
    a, b, c, dd = A[:, 0, 0], A[:, 0, 1], A[:, 1, 0], A[:, 1, 1]
    det = a * dd - b * c
    inv = np.empty_like(A)
    inv[:, 0, 0] = dd / det
    inv[:, 0, 1] = -b / det
    inv[:, 1, 0] = -c / det
    inv[:, 1, 1] = a / det
    return inv, det


def Jparts(n, edges, t, noexp=False):
    """t: (M, 2n) edge times.  Returns (J0, J1, Epar) with J = J0 + J1 (tau already integrated);
    noexp=True: the common factor exp(-Epar) is omitted (P0 = J0, P1 = J1)."""
    M = t.shape[0]
    E = len(edges)
    J0 = np.zeros(M)
    J1 = np.zeros(M)
    Epar_out = None
    for sigma in itertools.product((0, 1), repeat=E):
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
        w = 1.0 / t
        Apar = np.einsum('me,ei,ej->mij', w, a, a)
        bpar = np.einsum('me,ei,e->mi', w, a, c)
        cpar = np.einsum('me,e->m', w, c * c)
        Ainv, detpar = _inv_det(Apar)
        mean = -np.einsum('mij,mj->mi', Ainv, bpar)
        Epar = cpar - np.einsum('mi,mij,mj->m', bpar, Ainv, bpar)
        Epar_out = Epar
        covpar = 0.5 * Ainv
        Aperp = np.einsum('me,ei,ej->mij', w, b, b)
        Pinv, detperp = _inv_det(Aperp)
        S = 0.5 * Pinv                                  # perpendicular covariance
        pref = sgn * np.prod(1.0 / (2 * np.pi * t ** 2), axis=1)
        pref = pref * (np.pi ** (n / 2) / np.sqrt(detpar)) * (1.0 if noexp else np.exp(-Epar)) * (np.pi ** (n / 2) / np.sqrt(detperp))
        if n == 2:
            S00, S01, S11 = S[:, 0, 0], S[:, 0, 1], S[:, 1, 1]
            detS = S00 * S11 - S01 ** 2
            rho = S01 / np.sqrt(S00 * S11)
            Es = (2 / np.pi) * np.arcsin(np.clip(rho, -1, 1))
            sq = np.sqrt(np.maximum(detS, 0))
            Eyy = np.empty((M, 2, 2))
            Eyy[:, 0, 0] = S00 * Es + 2 * S01 * sq / (np.pi * S11)
            Eyy[:, 1, 1] = S11 * Es + 2 * S01 * sq / (np.pi * S00)
            Eyy[:, 0, 1] = S01 * Es + 2 * sq / np.pi
            Eyy[:, 1, 0] = Eyy[:, 0, 1]
        for A in itertools.combinations(range(E), n):
            Mat = np.zeros((E, 2 * n))
            for e in range(E):
                if e in A:
                    Mat[e, 1::2] = b[e]
                else:
                    Mat[e, 0::2] = a[e]
            dcoef = np.linalg.det(Mat)
            if abs(dcoef) < 1e-12:
                continue
            Ac = [e for e in range(E) if e not in A]
            # parallel moments: E[prod_{e in A} (a_e.x + c_e)]
            mu = [np.einsum('mi,i->m', mean, a[e]) + c[e] for e in A]
            m0 = np.prod(mu, axis=0) if n > 0 else np.ones(M)
            if n == 2:
                m1 = np.einsum('i,mij,j->m', a[A[0]], covpar, a[A[1]])
            else:
                m1 = np.zeros(M)
            # perpendicular: E[prod_{e in Ac} (b_e.y) prod_v sgn(y_v)/2]
            if n == 1:
                be = b[Ac[0], 0]
                mperp = be * np.sqrt(2 * S[:, 0, 0] / np.pi) / 2
            else:
                b1, b2 = b[Ac[0]], b[Ac[1]]
                mperp = np.einsum('i,mij,j->m', b1, Eyy, b2) / 4
            fac = dcoef * ((-1) ** (E - n)) * pref * mperp
            J0 = J0 + fac * m0
            J1 = J1 + fac * m1
    return J0, J1, Epar_out


def tanh_sinh(h, K):
    k = np.arange(-K, K + 1) * h
    x = 0.5 * (1 + np.tanh(0.5 * np.pi * np.sinh(k)))
    wgt = 0.5 * h * 0.5 * np.pi * np.cosh(k) / np.cosh(0.5 * np.pi * np.sinh(k)) ** 2
    good = (x > 1e-15) & (x < 1 - 1e-15)
    return x[good], wgt[good]


def weight(n, edges, h=0.1, K=40, check=True):
    a = 2.5 if n == 1 else 5.0
    d = 2 * n
    x, wq = tanh_sinh(h, K)
    if n == 1:
        u = np.stack([x, 1 - x], axis=1)
        W = wq
    else:
        S1, S2, S3 = np.meshgrid(x, x, x, indexing='ij')
        W1, W2, W3 = np.meshgrid(wq, wq, wq, indexing='ij')
        S1, S2, S3, W1, W2, W3 = [z.ravel() for z in (S1, S2, S3, W1, W2, W3)]
        u1 = S1
        u2 = (1 - S1) * S2
        u3 = (1 - S1) * (1 - S2) * S3
        u4 = (1 - S1) * (1 - S2) * (1 - S3)
        u = np.stack([u1, u2, u3, u4], axis=1)
        W = W1 * W2 * W3 * (1 - S1) ** 2 * (1 - S2)
    out = 0.0
    p1max = 0.0
    p0max = 0.0
    B = 200000
    for i in range(0, u.shape[0], B):
        uu = u[i:i + B]
        P0, P1, E = Jparts(n, edges, uu, noexp=True)
        p1max = max(p1max, np.max(np.abs(P1) / (np.abs(P0) + np.abs(P1) + 1e-300)))
        val = P0 * E ** (d - a)
        out += np.sum(np.where(np.isfinite(val), val, 0.0) * W[i:i + B])
    return math.gamma(a - d) * out, p1max


def weight_sectors(n, edges, h=0.15, K=None):
    """same integral, Hepp-sector coordinates: in the sector t_pi(1) >= ... >= t_pi(d), t_pi(1) = 1,
    t_pi(k+1) = t_pi(k) s_k, s in (0,1)^{d-1};  int_Delta G du = sum_pi int G(t(s)) prod_k s_k^{d-1-k} ds,
    G = P0 E^{d-a} homogeneous of degree -d."""
    a = 2.5 if n == 1 else 5.0
    d = 2 * n
    K = int(3.5 / h) if K is None else K
    x, wq = tanh_sinh(h, K)
    x = np.clip(x, 1e-10, 1 - 1e-10)
    grids = np.meshgrid(*([x] * (d - 1)), indexing='ij')
    wgr = np.meshgrid(*([wq] * (d - 1)), indexing='ij')
    S = np.stack([g.ravel() for g in grids], axis=1)
    W = np.prod(np.stack([g.ravel() for g in wgr], axis=1), axis=1)
    cum = np.ones((S.shape[0], d))
    for k in range(1, d):
        cum[:, k] = cum[:, k - 1] * S[:, k - 1]
    jac = np.ones(S.shape[0])
    # d^{d-1}t = prod_k s_k^{d-1-k} ds  (k = 1..d-1)
    for k in range(1, d):
        jac = jac * S[:, k - 1] ** (d - 1 - k)
    total = 0.0
    for perm in itertools.permutations(range(d)):
        t = np.empty_like(cum)
        for pos, e in enumerate(perm):
            t[:, e] = cum[:, pos]
        P0, P1, E = Jparts(n, edges, t, noexp=True)
        G = P0 * E ** (d - a)
        total += np.sum(np.where(np.isfinite(G), G, 0.0) * W * jac)
    return math.gamma(a - d) * total
