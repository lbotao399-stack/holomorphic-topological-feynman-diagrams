"""(b) Gaussian integrals over vertex positions: numerical quadrature on the smallest pieces, exact
identities for the theta graphs.

(F) one Fourier step, one holomorphic component:
    int d^2x (4 pi t)^{-1} e^{-|z-w|^2/4t} e^{lam z} zb^k = e^{lam w} (wb + 4 t lam)^k .
(V) one vertex v joined to fixed points a_e (e = 1..3) by heat edges t_e, Lam = sum 1/t_e:
    int d^4x e^{-sum_e |z-a_e|^2/4t_e + lam.z} zb1^k1 zb2^k2 = (4pi/Lam)^2 e^{lam.c + E0} prod_i (cb_i + 4 lam_i/Lam)^{k_i},
    int dtau e^{-sum_e (tau-al_e)^2/4t_e} tau^k = (4pi/Lam)^{1/2} e^{E0'} E[(taub + xi)^k],  xi ~ N(0, 2/Lam).
(M) two integrated vertices, one component, generic Lam (2x2):
    int d^2x1 d^2x2 e^{-(1/4) z^dag Lam z + lam.z} G(zb) = (4 pi)^2 / det Lam * G(4 Lam^{-1} lam).
(T) theta graphs, exact rationals: det Lam = U / prod t;  <tau tau^T> = 2 Lam^{-1} gives
    <w_e w_f> = (1/2)[delta_ef/t_e - (C^T M^{-1} C)_ef];  zb = 4 Lam^{-1} lam gives the Kirchhoff currents of htk;
    <psi^tau_1 psi^tau_2> = -(1/2) sum_ab (M^{-1})_ab dM_1a ^ dM_b2  (S3 prop-rho).
"""
import sys
import random
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import brute

rng = np.random.default_rng(7)


def gh2(f, cx, cy, s, m=80):
    """int_{R^2} f(x, y) dx dy, Gauss-Hermite with weight e^{-(u^2+v^2)}, x = cx + s u, y = cy + s v."""
    u, w = np.polynomial.hermite.hermgauss(m)
    U, V = np.meshgrid(u, u, indexing='ij')
    W = np.outer(w, w)
    x, y = cx + s * U, cy + s * V
    return np.sum(W * f(x, y) * np.exp(U ** 2 + V ** 2)) * s * s


def digits(a, b):
    return -np.log10(abs(a - b) / abs(b)) if a != b else 16.0


print("(F) one Fourier step")
t = 0.45
w0 = 0.3 + 0.2j
lam = 0.7 - 0.4j
for k in range(4):
    f = lambda x, y: np.exp(-abs(x + 1j * y - w0) ** 2 / (4 * t)) / (4 * np.pi * t) * np.exp(lam * (x + 1j * y)) * (x - 1j * y) ** k
    num = gh2(f, w0.real, w0.imag, 2 * np.sqrt(t))
    ex = np.exp(lam * w0) * (np.conj(w0) + 4 * t * lam) ** k
    print("  k=%d: quadrature %s  exact %s  agreeing digits %.1f" % (k, np.round(num, 13), np.round(ex, 13), digits(num, ex)))

print("(V) one vertex with three heat edges")
ts = np.array([0.3, 0.55, 0.8])
a1 = np.array([0.2 + 0.1j, -0.4 + 0.3j, 0.1 - 0.5j])     # component 1 positions of the neighbours
a2 = np.array([-0.3 + 0.2j, 0.5 - 0.1j, 0.2 + 0.4j])     # component 2
lam1, lam2 = 0.6 + 0.3j, -0.5 + 0.2j
Lam = np.sum(1 / ts)
for (k1, k2) in [(0, 0), (1, 0), (1, 1), (2, 1)]:
    tot = 1.0 + 0j
    ex = (4 * np.pi / Lam) ** 2
    for (a, lm, k) in [(a1, lam1, k1), (a2, lam2, k2)]:
        c = np.sum(a / ts) / Lam
        E0 = -0.25 * np.sum(abs(a) ** 2 / ts) + 0.25 * Lam * abs(c) ** 2
        f = lambda x, y, a=a, lm=lm, k=k: np.exp(-sum(abs(x + 1j * y - a[e]) ** 2 / (4 * ts[e]) for e in range(3))
                                                 + lm * (x + 1j * y)) * (x - 1j * y) ** k
        tot *= gh2(f, c.real, c.imag, 2 / np.sqrt(Lam))
        ex *= np.exp(lm * c + E0) * (np.conj(c) + 4 * lm / Lam) ** k
    print("  (k1,k2)=(%d,%d): quadrature %s exact %s digits %.1f" % (k1, k2, np.round(tot, 12), np.round(ex, 12), digits(tot, ex)))
al = np.array([0.4, -0.7, 0.25])
taub = np.sum(al / ts) / Lam
E0 = -0.25 * np.sum(al ** 2 / ts) + 0.25 * Lam * taub ** 2
u, wq = np.polynomial.hermite.hermgauss(80)
for k in range(3):
    s = 2 / np.sqrt(Lam)
    xx = taub + s * u
    num = np.sum(wq * np.exp(u ** 2) * np.exp(-sum((xx - al[e]) ** 2 / (4 * ts[e]) for e in range(3))) * xx ** k) * s
    mom = [1, taub, taub ** 2 + 2 / Lam][k]
    ex = np.sqrt(4 * np.pi / Lam) * np.exp(E0) * mom
    print("  tau, k=%d: quadrature %.15f exact %.15f digits %.1f" % (k, num, ex, digits(num, ex)))

print("(M) two integrated vertices, one holomorphic component")
L = np.array([[1 / 0.3 + 1 / 0.6 + 1 / 0.9, -1 / 0.6], [-1 / 0.6, 1 / 0.6 + 1 / 0.45]])   # nodal Laplacian
lv = np.array([0.5 - 0.2j, -0.3 + 0.4j])
m = 34
u, wq = np.polynomial.hermite.hermgauss(m)
# reference variables: z_v = s_v (u + i v)
evals = np.linalg.eigvalsh(L)
s = 2 / np.sqrt(evals.min())
U1, V1, U2, V2 = np.meshgrid(u, u, u, u, indexing='ij')
W = np.einsum('i,j,k,l->ijkl', wq, wq, wq, wq) * np.exp(U1 ** 2 + V1 ** 2 + U2 ** 2 + V2 ** 2)
z1, z2 = s * (U1 + 1j * V1), s * (U2 + 1j * V2)
quad = -0.25 * (L[0, 0] * abs(z1) ** 2 + 2 * L[0, 1] * (z1 * np.conj(z2)).real + L[1, 1] * abs(z2) ** 2)
base = np.exp(quad + lv[0] * z1 + lv[1] * z2)
phi = 4 * np.linalg.solve(L, lv)
for name, G, Gex in [("1", lambda a, b: 1 + 0 * a, 1), ("zb1", lambda a, b: a, phi[0]), ("zb1 zb2", lambda a, b: a * b, phi[0] * phi[1]),
                     ("zb2^2 zb1", lambda a, b: b * b * a, phi[1] ** 2 * phi[0])]:
    num = np.sum(W * base * G(np.conj(z1), np.conj(z2))) * s ** 4
    ex = (4 * np.pi) ** 2 / np.linalg.det(L) * Gex
    print("  G = %-9s: quadrature %s exact %s digits %.1f" % (name, np.round(num, 11), np.round(ex, 11), digits(num, ex)))

print("(T) theta graphs, exact rational identities")
R = random.Random(5)


def fr_inv(Mx):
    return sp.Matrix(Mx).inv()


for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, inter = htk.theta(*pqr)
    N, nv = G.ne, G.nv
    C = sp.Matrix(G.cycle_matrix())
    B = sp.Matrix(G.incidence())          # N x (nv-1)
    oks = []
    for trial in range(3):
        tt = [sp.Rational(R.randint(1, 30), R.randint(1, 30)) for _ in range(N)]
        Tm = sp.diag(*tt)
        Lm = B.T * Tm.inv() * B
        M = C * Tm * C.T
        U = M.det()
        ok1 = sp.simplify(Lm.det() - U / sp.prod(tt)) == 0
        lhs = sp.Rational(1, 2) * Tm.inv() * B * Lm.inv() * B.T * Tm.inv()   # <w w^T> from <tau tau^T> = 2 Lam^{-1}
        rhs = sp.Rational(1, 2) * (Tm.inv() - C.T * M.inv() * C)
        ok2 = (lhs - rhs) == sp.zeros(N, N)
        lam = [[sp.Rational(R.randint(-9, 9)), sp.Rational(R.randint(-9, 9))] for _ in range(nv - 1)]
        y_engine = brute.currents(G, tt, lam)
        ok3 = True
        for i in range(2):
            phi = 4 * Lm.inv() * sp.Matrix([lam[v][i] for v in range(nv - 1)])   # zb_v
            y = [sum(B[e, v] * phi[v] for v in range(nv - 1)) / (4 * tt[e]) for e in range(N)]
            ok3 &= all(sp.simplify(y[e] - y_engine[e][i]) == 0 for e in range(N))
        # <psi^tau_1 psi^tau_2>: coefficient of dt_e ^ dt_f (e<f)
        W2 = rhs
        ok4 = True
        Mi = M.inv()
        for e in range(N):
            for f in range(e + 1, N):
                a = C[0, e] * C[1, f] * W2[e, f] - C[0, f] * C[1, e] * W2[f, e]
                # -(1/2) sum_ab Mi_ab dM_1a ^ dM_b2, dM_ra = sum_e C_re C_ae dt_e
                b = 0
                for aa in range(2):
                    for bb in range(2):
                        b += Mi[aa, bb] * (C[0, e] * C[aa, e] * C[bb, f] * C[1, f] - C[0, f] * C[aa, f] * C[bb, e] * C[1, e])
                b = -sp.Rational(1, 2) * b
                ok4 &= sp.simplify(a - b) == 0
        oks.append((ok1, ok2, ok3, ok4))
    print("  theta%s: [det Lam = U/prod t, <ww> formula, zb = 4 Lam^-1 lam gives htk currents, <psi psi> = rho-form] at 3 random t:" % (pqr,), oks)
