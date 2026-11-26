"""V4: end-to-end position-space integrals (heat-edge sum D wedge P |_UV computed directly, no Schwinger
fiber integration) for graphs with the parity structure of the flagship (P even, K and Theta odd), compared with
c * K_htk where c is the pointwise constant found in V2.
(0,1): P(tau) = int_eps^oo D*K_t dt = (1/2) erf(tau/(2 sqrt eps)), K_eps = (4 pi eps)^{-1/2} e^{-tau^2/4eps} dtau.
       banana(k), root at 0: G = int dtau D(P^k)|_UV = k int P^{k-1} K_eps dtau.  Exact value 2^{-k}(1-(-1)^k).
(1,0): Omega = dz, P = (1 - e^{-|z|^2/4eps})/((-2i) pi z), [K_eps] = e^{-|z|^2/4eps}/((-2i) 4 pi eps) dzb,
       G = int Omega e^{lam z} D(P^k)|_UV = k (-2i) int d^2x e^{lam z} P^{k-1} K.  Exact value lam^{k-1}/((k-1)! (-2 pi i)^{k-1}).
"""
import sys
import math
import cmath
import numpy as np
from scipy import integrate, special
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk

print("V4  end-to-end calibration integrals")
c_V2 = {('(0,1)', 3): 1 / (4 * math.pi), ('(0,1)', 5): 1 / (16 * math.pi**2),
        ('(1,0)', 2): 1j / (2 * math.pi), ('(1,0)', 3): -1 / (4 * math.pi**2)}
for k in (3, 5):
    K = htk.kernel_exact(htk.banana(k), None, T=1)
    Kf = float(sp.sympify(str(K)).subs('pi', sp.pi)) if not isinstance(K, (int, float)) else K
    for eps in (0.3, 0.01):
        s = 2 * math.sqrt(eps)
        f = lambda tau: k * (0.5 * special.erf(tau / s)) ** (k - 1) * math.exp(-tau**2 / (4 * eps)) / math.sqrt(4 * math.pi * eps)
        val, err = integrate.quad(f, -np.inf, np.inf, epsabs=1e-15, epsrel=1e-14)
        pred = c_V2[('(0,1)', k)] * Kf
        print("  (0,1) banana(%d) eps=%.2f: end-to-end %.15f ; exact 2^-k(1-(-1)^k) = %.15f ; c_V2 * K_htk = %.15f (K_htk = %s) ; agreeing digits %.1f"
              % (k, eps, val, 2.0**-k * 2, pred, K, -math.log10(abs(val - pred) / abs(pred) + 1e-17)))
for k in (2, 3):
    lam = sp.Symbol('l')
    K = htk.kernel_exact(htk.banana(k), [[lam]], T=0)
    for lv in (0.8 + 0.3j, -0.5 + 1.1j):
        for eps in (0.2, 0.05):
            # angular integral done analytically: only z^{k-1} of e^{lam z} survives against z^{-(k-1)}
            def radial(r):
                return (1 - math.exp(-r * r / (4 * eps))) ** (k - 1) * math.exp(-r * r / (4 * eps)) * 2 * math.pi * r
            rad, err = integrate.quad(radial, 0, np.inf, epsabs=1e-15, epsrel=1e-14)
            val = k * (-2j) * lv ** (k - 1) / math.factorial(k - 1) * rad / (((-2j) * math.pi) ** (k - 1) * (-2j) * 4 * math.pi * eps)
            # brute 2d check of the angular step at one point
            if eps == 0.2 and lv == 0.8 + 0.3j:
                def integrand(th, r, part):
                    z = r * cmath.exp(1j * th)
                    P = (1 - math.exp(-r * r / (4 * eps))) / ((-2j) * math.pi * z)
                    Kc = math.exp(-r * r / (4 * eps)) / ((-2j) * 4 * math.pi * eps)
                    v = k * (-2j) * cmath.exp(lv * z) * P ** (k - 1) * Kc * r
                    return v.real if part == 0 else v.imag
                re_ = integrate.dblquad(lambda th, r: integrand(th, r, 0), 0, 12, 0, 2 * math.pi, epsabs=1e-13, epsrel=1e-12)[0]
                im_ = integrate.dblquad(lambda th, r: integrand(th, r, 1), 0, 12, 0, 2 * math.pi, epsabs=1e-13, epsrel=1e-12)[0]
                print("     (1,0) banana(%d) 2d quadrature without the angular shortcut: %s vs %s" % (k, complex(re_, im_), val))
            exact = lv ** (k - 1) / (math.factorial(k - 1) * (-2j * math.pi) ** (k - 1))
            Kv = complex(sp.sympify(K).subs(lam, sp.sympify(lv)))
            pred = c_V2[('(1,0)', k)] * Kv
            print("  (1,0) banana(%d) lam=%s eps=%.2f: end-to-end %s ; exact %s ; c_V2 * K_htk = %s (K_htk = %s) ; digits %.1f"
                  % (k, lv, eps, val, exact, pred, K, -math.log10(abs(val - pred) / abs(pred) + 1e-17)))
