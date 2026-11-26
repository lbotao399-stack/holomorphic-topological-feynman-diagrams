"""(d) Holomorphic derivative on a propagator end = current insertion; end-to-end numerical check.

Pointwise (normalization b2, edge form):  Theta_e ~ exp(-eta_e^2 - y_e.(z_{e(0)} - z_{e(1)} + z_e)),
y_e = (zb_{e(0)} - zb_{e(1)})/(4 t_e), so
    d/dz^i_{e(0)} Theta_e = -y^i_e Theta_e,     d/dz^i_{e(1)} Theta_e = +y^i_e Theta_e,
and after the holomorphic Gaussian (normalization b6: G(zb) -> G(4 Lambda^{-1} lambda)) y_e becomes the htk
Kirchhoff current.  Prediction for the raw integral with derivatives:
    G[prod d] = vs_ext s1 s2 c^{-n} K_htk(powers) prod_(e,end) (-1 if end = 0 else +1).
Check: (1,0) bananas (X = C, Omega = dz, v = -2i), plane wave e^{lam z} at vertex 0, root = vertex 1 at 0,
the configuration of normalization n4 with derivatives added:
    k_t(z) = e^{-|z|^2/4t}/(4 pi t)  (density of Omega ^ [K_t]),  P(z) = (1 - e^{-|z|^2/4eps})/(pi v z),
    G[ins] = int d^2x e^{lam z} sum_j prod_i X_i,   X_j = D_j k_eps, X_i = D_i P (i != j),
    D_i = d_z (derivative at vertex 0 on edge i), -d_z (derivative at the root on edge i), 1 (none).
"""
import sys
import numpy as np
import sympy as sp
from scipy import integrate
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import conventions

v = -2j


def pieces(eps):
    def k(z, r2):
        return np.exp(-r2 / (4 * eps)) / (4 * np.pi * eps)

    def dk(z, r2):
        return -(np.conj(z) / (4 * eps)) * k(z, r2)

    def P(z, r2):
        return (1 - np.exp(-r2 / (4 * eps))) / (np.pi * v * z)

    def dP(z, r2):
        E = np.exp(-r2 / (4 * eps))
        return ((np.conj(z) / (4 * eps)) * E * z - (1 - E)) / (np.pi * v * z * z)
    return k, dk, P, dP


def position_integral(k_edges, ders, lam, eps):
    """ders[i] in {None, 'v0', 'root'}."""
    k, dk, P, dP = pieces(eps)

    def X(i, heat, z, r2):
        d = ders[i]
        if heat:
            base = dk(z, r2) if d else k(z, r2)
        else:
            base = dP(z, r2) if d else P(z, r2)
        return -base if d == 'root' else base

    def f(rr, ph, part):
        z = rr * np.exp(1j * ph)
        r2 = rr * rr
        tot = 0
        for j in range(k_edges):
            term = 1
            for i in range(k_edges):
                term = term * X(i, i == j, z, r2)
            tot = tot + term
        val = np.exp(lam * z) * tot * rr
        return val.real if part == 0 else val.imag
    R = 14 * np.sqrt(eps)
    re = integrate.dblquad(lambda ph, rr: f(rr, ph, 0), 0, R, 0, 2 * np.pi, epsabs=1e-13, epsrel=1e-12)[0]
    im = integrate.dblquad(lambda ph, rr: f(rr, ph, 1), 0, R, 0, 2 * np.pi, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j * im


lam_s = sp.Symbol('l')
cases = [
    (2, ['v0', None], {(0, 0): 1}, -1),
    (2, [None, 'root'], {(1, 0): 1}, +1),
    (2, ['v0', 'v0'], {(0, 0): 1, (1, 0): 1}, +1),
    (3, ['v0', None, None], {(0, 0): 1}, -1),
    (3, ['v0', 'root', None], {(0, 0): 1, (1, 0): 1}, -1),
    (3, [None, 'root', 'root'], {(1, 0): 1, (2, 0): 1}, +1),
]
print(__doc__)
for k_edges, ders, powers, endsign in cases:
    G = htk.banana(k_edges)
    n = k_edges - 1
    K = sp.expand(htk.kernel_exact(G, [[lam_s]], T=0, powers=powers))
    vs = conventions.extraction_sign_formula(G, 1) * (-1) ** (n * 1)
    for lam in (0.8 + 0.3j, -0.5 + 1.1j):
        Knum = complex(K.subs(lam_s, lam))
        pred = endsign * vs * (v * np.pi) ** (-n) * Knum
        for eps in (0.2, 0.05):
            val = position_integral(k_edges, ders, lam, eps)
            dig = -np.log10(abs(val - pred) / abs(pred))
            ctrl = -pred
            print("banana(%d) ders=%s powers=%s end-sign %+d lam=%s eps=%.2f: position %s ; K_htk = %s ; prediction %s ; digits %.1f ; "
                  "control (opposite end sign) rel.diff %.2f"
                  % (k_edges, ders, powers, endsign, lam, eps, np.round(val, 13), K, np.round(pred, 13), dig,
                     abs(val - ctrl) / abs(pred)))
