"""Two independent massless two-loop reductions, with the kinetic coefficient in (1.1.2)."""
from itertools import product
import sympy as s


def color():
    # For SU(2) in an orthonormal real basis check the Gram identities derived
    # from f_{amn} f_b^{mn}=C_A g_ab and C_1+C_2+C_3=0 for any simple algebra.
    def f(i,j,k):
        return s.LeviCivita(i,j,k)
    gram=s.zeros(3)
    for r in product(range(3),repeat=3):
        v=s.Matrix([sum(f(0,r[i],m)*f(r[(i+1)%3],r[(i+2)%3],m) for m in range(3))
                    for i in range(3)])
        assert sum(v)==0
        gram+=v*v.T
    assert gram==6*s.eye(3)-2*s.ones(3)


def moments():
    nu=s.symbols('nu',positive=True)
    # G(x) is proportional to (x^2)^(-nu), with nu=d/2-1. For its Fourier
    # transform J_alpha, strip the explicit p^2 in J_(alpha+1)/J_alpha.
    def fourier_ratio(a):
        return 1/(4*a*(nu-a))
    q=2*nu+1
    coordinate=s.factor(4*nu*fourier_ratio(3*nu)*(1+4*q*(2*q-1)*fourier_ratio(3*nu+1)))
    # Perform the (k_2,k_3) tensor bubble first, then reduce powers of its momentum.
    d=2*nu+2
    def bubble_ratio(n):
        return (n+1)/(n+d-1)
    n=nu-1
    momentum=s.factor((1-2*bubble_ratio(n)+bubble_ratio(n)*bubble_ratio(n+1))/(4*(d-1)))
    assert coordinate==momentum and s.simplify(momentum-1/(6*(3*nu+1)))==0
    assert momentum.subs(nu,1)==s.Rational(1,24)
    # The differentiated kinetic quartic is (1/6) sum C_i X_i. After the
    # common integral its squared color norm is C_A^2/8 times the X_1 moment.
    vertex_norm=s.Rational(1,6)**2*s.Rational(3,2)*3
    scalar_residue=-s.Rational(1,4)  # in units of p^2/(4*pi)^4
    scalar_x=4*momentum.subs(nu,1)  # X_1=2p.(k_2-k_3) modulo canceled propagators
    inverse_propagator=-s.Rational(1,6)*s.Rational(1,2)**3*vertex_norm*scalar_x*scalar_residue
    assert inverse_propagator==s.Rational(1,9216)
    return inverse_propagator


def cubic_test():
    # Complexified momenta are given by independent holomorphic and
    # antiholomorphic components; no reality condition is imposed in this test.
    v=[s.Matrix(i) for i in ((1,0),(0,1),(-1,-1))]
    w=[s.Matrix(i) for i in ((0,1),(-1,0),(1,-1))]
    assert sum(v,s.zeros(2,1))==sum(w,s.zeros(2,1))==s.zeros(2,1)
    assert all(v[i].dot(w[j])+v[j].dot(w[i])==0 for i,j in product(range(3),repeat=2))
    omega=lambda i,j:v[i].dot(w[j])-v[j].dot(w[i])
    assert all(omega(i,(i+1)%3)==-2 for i in range(3))
    M=[a*b.T for a,b in zip(v,w)]
    assert all(m.trace()==0 for m in M)
    assert (M[0]*(M[1]*M[2]-M[2]*M[1])).trace()==2


if __name__=='__main__':
    color()
    cubic_test()
    print('Negative scalar sunset pole; two tensor reductions agree. Inverse-propagator residue:',moments())
