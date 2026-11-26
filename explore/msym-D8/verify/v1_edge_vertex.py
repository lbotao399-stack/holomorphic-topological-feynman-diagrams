"""C2 (edge and vertex factors on C^4) and the measure convention of C1; own exterior algebra (sympy coefficients).

(a) dx^1..dx^d ^ dxbar^1..dxbar^d in real coordinates (x^a = u^a + i v^a), d = 2, 3, 4.
(b) [K_t]_Omega = e^{-x.xbar/2t} / (2^d (2 pi t)^d) dxbar^{1..d}  (S2 eq-heat for d = 3, 2^3 = 8);
    int Omega ^ [K_t] = 1 with Omega ^ dxbar^{1..d} = 2^d d^{2d}x_E (Gaussian integral, exact).
(c) [K_t] - dt ^ dbar^*[K_t] == (2 pi)^{-d} e^{-y.x} d^d y, y = xbar/(2t), dbar^* = -2 sum_a iota_{d/dxbar^a} d/dx^a,
    d = 3 (S2 eq-y: 1/(8 pi^3)) and d = 4 (claim: 1/(16 pi^4)).
(d) vertex factor: Omega = c_Omega dx^{1..d}, x = i s (S2 eq-vertex recipe): c_Omega i^d (2 pi)^d.
"""
import sympy as sp
import itertools


class Form:
    def __init__(self, terms=None):
        self.t = {k: v for k, v in (terms or {}).items() if v != 0}

    @staticmethod
    def gen(i, c=1):
        return Form({(i,): sp.sympify(c)})

    @staticmethod
    def scal(c):
        return Form({(): sp.sympify(c)})

    def __add__(self, o):
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = sp.simplify(r.get(k, 0) + v)
        return Form(r)

    def __neg__(self):
        return Form({k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-o)

    def __xor__(self, o):
        r = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                if set(k1) & set(k2):
                    continue
                seq = list(k1 + k2)
                sign = 1
                for i in range(len(seq)):
                    for j in range(i + 1, len(seq)):
                        if seq[i] > seq[j]:
                            sign = -sign
                key = tuple(sorted(seq))
                r[key] = r.get(key, 0) + sign * v1 * v2
        return Form(r)

    def interior(self, i):
        """iota_{d/d(gen i)} (left contraction)."""
        r = {}
        for k, v in self.t.items():
            if i in k:
                p = k.index(i)
                r[k[:p] + k[p + 1:]] = r.get(k[:p] + k[p + 1:], 0) + (-1) ** p * v
        return Form(r)

    def coeff_map(self, f):
        return Form({k: f(v) for k, v in self.t.items()})

    def equals(self, o):
        d = self - o
        return all(sp.simplify(v) == 0 for v in d.t.values())


# (a)
for d in (2, 3, 4):
    # generators: du^a = 2a, dv^a = 2a+1
    dx = [Form.gen(2 * a) + Form.gen(2 * a + 1, sp.I) for a in range(d)]
    dxb = [Form.gen(2 * a) + Form.gen(2 * a + 1, -sp.I) for a in range(d)]
    w = Form.scal(1)
    for f in dx + dxb:
        w = w ^ f
    vol = tuple(range(2 * d))
    print('(a) d = %d: dx^{1..d} ^ dxbar^{1..d} = %s d^{%d}x_E' % (d, sp.simplify(w.t.get(vol, 0)), 2 * d))

# (b), (c)
t = sp.Symbol('t', positive=True)
for d in (3, 4):
    x = sp.symbols('x1:%d' % (d + 1))
    xb = sp.symbols('xb1:%d' % (d + 1))
    # generators: dt = 0, dxbar^a = a+1
    E = sp.exp(-sum(x[a] * xb[a] for a in range(d)) / (2 * t))
    K = Form({tuple(range(1, d + 1)): E / (2 ** d * (2 * sp.pi * t) ** d)})
    dstar = Form()
    for a in range(d):
        dstar = dstar + K.interior(a + 1).coeff_map(lambda c, a=a: -2 * sp.diff(c, x[a]))
    lhs = K - (Form.gen(0) ^ dstar)
    dy = [Form.gen(a + 1, 1 / (2 * t)) + Form.gen(0, -xb[a] / (2 * t ** 2)) for a in range(d)]
    ddy = Form.scal(1)
    for f in dy:
        ddy = ddy ^ f
    Ey = sp.exp(-sum(xb[a] / (2 * t) * x[a] for a in range(d)))
    rhs = ddy.coeff_map(lambda c: c * Ey / (2 * sp.pi) ** d)
    print('(c) d = %d: [K_t] - dt ^ dbar^*[K_t] == (2 pi)^{-%d} e^{-y.x} d^%dy, (2 pi)^{-%d} = %s: %s'
          % (d, d, d, d, sp.simplify(1 / (2 * sp.pi) ** d), lhs.equals(rhs)))
    # normalization: x.xbar = |x|^2_E ; int e^{-|x|^2/2t} d^{2d}x = (2 pi t)^d
    s = sp.Symbol('s', real=True)
    g1 = sp.integrate(sp.exp(-s ** 2 / (2 * t)), (s, -sp.oo, sp.oo))
    gauss = sp.simplify(g1 ** (2 * d))
    norm = sp.simplify(2 ** d * gauss / (2 ** d * (2 * sp.pi * t) ** d))
    print('(b) d = %d: int Omega ^ [K_t] = 2^%d * (2 pi t)^%d / (2^%d (2 pi t)^%d) = %s' % (d, d, d, d, d, norm))
    # parity: dbar^*[K_t] is a (0, d-1)-form
    degs = sorted({len(k) for k in dstar.t})
    print('    dbar^*[K_t] has form degree %s -> propagator parity %s' % (degs, 'odd' if degs[0] % 2 else 'even'))

# (d)
for d, cO in ((3, sp.I), (4, 1)):
    print('(d) d = %d: c_Omega i^d (2 pi)^d = %s' % (d, sp.simplify(cO * sp.I ** d * (2 * sp.pi) ** d)))
