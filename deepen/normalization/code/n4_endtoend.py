"""(b) End-to-end numerical checks of  G = vs_ext s1 s2 c^{-n} K  on the graphs where the
position integral is low-dimensional (same parity structure as (H,T) = (2,1): P even, K and Theta odd).

s1 = (-1)^{(Dn)(D(|V|-1))} (Psi past [dX top]), s2 = (-1)^{C(|V|-1,2) C(D,2)} (reordering to prod_v),
T = 1:  c = 2 v_H pi^{H+1/2} ;  T = 0:  c = v_H pi^H,
v_H defined by Omega ^ dtau^T ^ d^H zb = v_H d^{2H+T}x.   (H,T) = (2,1): v = 4, (8 pi^{5/2})^{-n}.

(0,1) banana(k): X = R, [K_t] = e^{-tau^2/4t}(4 pi t)^{-1/2} dtau, P_{eps,oo} = (1/2) erf(tau/(2 sqrt eps)),
      G = int_R D(P^k)|_UV = k int P^{k-1} K_eps.
(1,0) banana(k): X = C, Omega = dz, v_1 = -2i, [K_t] = (1/v) e^{-|z|^2/4t} (4 pi t)^{-1} dzb, D* = -4 iota_{d zb} d_z,
      P_{eps,oo} = (1 - e^{-|z|^2/4 eps}) / (pi v z),  G = int_C Omega e^{lam z} D(P^k)|_UV = k int Omega e^{lam z} P^{k-1} K_eps.
"""
import sys
import numpy as np
import sympy as sp
from scipy import integrate, special
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import conventions

print("(0,1) bananas, root = vertex 1 at 0, edges 0 -> 1:")
for k in (3, 5):
    n = k - 1
    G = htk.banana(k)
    K = htk.kernel_exact(G, None, T=1)
    Kf = float(K)
    vs = conventions.extraction_sign_formula(G, 1)
    s1 = (-1) ** (n * 1)
    vs = vs * s1
    pred = vs * (2 * np.sqrt(np.pi)) ** (-n) * Kf
    for eps in (0.3, 0.01):
        P = lambda x: 0.5 * special.erf(x / (2 * np.sqrt(eps)))
        Kd = lambda x: np.exp(-x * x / (4 * eps)) / np.sqrt(4 * np.pi * eps)
        val, err = integrate.quad(lambda x: k * P(x) ** (k - 1) * Kd(x), -np.inf, np.inf, epsabs=1e-15, epsrel=1e-14)
        print("  banana(%d), n=%d, eps=%.2f: position integral %.15f ; K_htk = %s, vs_ext*s1*s2 = %+d, prediction %.15f ; digits %.1f"
              % (k, n, eps, val, K, vs, pred, -np.log10(abs(val - pred) / abs(pred))))

print("(1,0) bananas, root = vertex 1 at 0, edges 0 -> 1, plane wave e^{lam z} at vertex 0:")
v = -2j
lam_s = sp.Symbol('l')
for k in (2, 3):
    n = k - 1
    G = htk.banana(k)
    K = sp.expand(htk.kernel_exact(G, [[lam_s]], T=0))
    vs = conventions.extraction_sign_formula(G, 1) * (-1) ** (n * 1)
    for lam in (0.8 + 0.3j, -0.5 + 1.1j):
        Knum = complex(K.subs(lam_s, lam))
        pred = vs * (v * np.pi) ** (-n) * Knum
        for eps in (0.2, 0.05):
            def f(rr, ph, part):
                zz = rr * np.exp(1j * ph)
                P = (1 - np.exp(-rr * rr / (4 * eps))) / (np.pi * v * zz)
                val = k * np.exp(lam * zz) * P ** (k - 1) * np.exp(-rr * rr / (4 * eps)) / (4 * np.pi * eps) * rr
                return val.real if part == 0 else val.imag
            re = integrate.dblquad(lambda ph, rr: f(rr, ph, 0), 0, 12 * np.sqrt(eps), 0, 2 * np.pi, epsabs=1e-13, epsrel=1e-12)[0]
            im = integrate.dblquad(lambda ph, rr: f(rr, ph, 1), 0, 12 * np.sqrt(eps), 0, 2 * np.pi, epsabs=1e-13, epsrel=1e-12)[0]
            val = re + 1j * im
            print("  banana(%d), n=%d, lam=%s, eps=%.2f: position integral %s ; K_htk(lam) = %s, vs_ext*s1*s2 = %+d, prediction %s ; digits %.1f"
                  % (k, n, lam, eps, np.round(val, 13), K, vs, np.round(pred, 13), -np.log10(abs(val - pred) / abs(pred))))
            pred0 = conventions.extraction_sign_formula(G, 1) * (v * np.pi) ** (-n) * Knum
            print("      control without s1: prediction %s ; |val - pred|/|pred| = %.2e" % (np.round(pred0, 13), abs(val - pred0) / abs(pred0)))
