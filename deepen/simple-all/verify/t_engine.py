"""Test of vpoint.Engine against the explicit index network (exact Fractions, then reduced mod p)."""
import sys, numpy as np, sympy as sp, itertools
from fractions import Fraction as Fr
import vlie, vpoint

def red(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator, p - 2, p) % p

def network(L, A):
    d = L.d
    F = L.f_dense()
    Ki = np.array(sp.Matrix(L.K).inv().tolist(), dtype=object)
    Ki = np.vectorize(lambda z: Fr(int(sp.Rational(z).p), int(sp.Rational(z).q)))(Ki)
    fA = [np.einsum('a,abc->bc', a, F) for a in A]          # fA[i][b,c] = f(A_i, e_b, e_c)
    # W133 = - f_xyz f_{A1 y'u} f_{A2 u'w} f_{A3 z's} f_{A4 s'r} f_{x'w'r'}
    Fu = np.einsum('xyz,xa,yb,zc->abc', F, Ki, Ki, Ki)        # f with all indices raised
    P_yw = Ki.dot(fA[0]).dot(Ki).dot(fA[1])                    # y(raised by contraction with f^{xyz}) ... careful below
    # explicit: sum f^{x y z}_raised? do it literally with einsum over all indices
    W133 = -np.einsum('xyz,Yu,Uw,Zs,Sr,XWR,xX,yY,zZ,uU,wW,sS,rR->', F, fA[0], fA[1], fA[2], fA[3], F, Ki, Ki, Ki, Ki, Ki, Ki, Ki, optimize=True)
    W223 = -np.einsum('xyz,Xu,Yw,Zs,Sr,UWR,xX,yY,zZ,uU,wW,sS,rR->', F, fA[0], fA[1], fA[2], fA[3], F, Ki, Ki, Ki, Ki, Ki, Ki, Ki, optimize=True)
    return W133, W223

p = 1048573
for key in sys.argv[1:]:
    L = vlie.build(key)
    E = vpoint.Engine(L, p)
    rng = np.random.default_rng(7)
    ok = True
    for trial in range(2):
        A = [np.array([Fr(int(v)) for v in rng.integers(-3, 4, L.d)], dtype=object) for _ in range(4)]
        w133, w223 = network(L, A)
        Ap = [np.array([red(v, p) for v in a], dtype=np.float64) for a in A]
        e133, e223 = E.W(Ap)
        ok &= (red(w133, p) == e133) and (red(w223, p) == e223)
        print(key, 'trial', trial, 'network W133 =', w133, ' W223 =', w223, ' engine mod p equal:', red(w133, p) == e133, red(w223, p) == e223)
    print(key, 'engine == explicit index network:', ok)
