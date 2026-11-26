"""Commuting normal deformations of the regular D3 Koszul factorization."""

from itertools import permutations

import sympy as sp

from quadratic_d3_factorization import exterior_matrices


def brane_obstruction():
    w = sp.symbols("x y z")
    u = sp.symbols("u1:4")
    b = sp.Symbol("b")
    potential = sum(x**3 for x in w)/6 + b*sp.prod(w)
    displacement = potential.xreplace(dict(zip(w, u)))
    wedge, contraction = exterior_matrices()
    unit = sp.eye(8)

    def verify(matrix):
        assert matrix.applyfunc(sp.expand) == sp.zeros(8)

    base = sum((x*h + potential.diff(x)*e/3
                for x,h,e in zip(w, contraction, wedge)), sp.zeros(8))
    # D_{-t u}, where D_v resolves the point w=v. The coefficient of t
    # has positive normal contraction. Its primitive at order two is even in u.
    first = sum((a*h - a*sum((potential.diff(x,y)*e/6 for y,e in zip(w,wedge)),
                              sp.zeros(8)) for x,a,h in zip(w,u,contraction)), sp.zeros(8))
    second = sum((displacement.diff(a)*e/3 for a,e in zip(u,wedge)), sp.zeros(8))
    verify(base**2-potential*unit)
    verify(base*first+first*base)
    verify(base*second+second*base+first**2)
    verify(first*second+second*first-displacement*unit)
    verify(second**2)
    assert base.subs(dict.fromkeys(w, 0)) == sp.zeros(8)
    for i,ui in enumerate(u):
        eta_i = first.diff(ui)
        verify(base*wedge[i]+wedge[i]*base-w[i]*unit)
        for j,uj in enumerate(u):
            eta_j = first.diff(uj)
            verify(eta_i*eta_j+eta_j*eta_i+potential.diff(w[i],w[j])*unit/3)
    t = sp.Symbol("t")
    moved = base-t*first+t**2*second
    on_ray = dict(zip(w, (t*a for a in u)))
    at_displacement = sum((t**2*displacement.diff(a)*e for a,e in zip(u,wedge)),
                          sp.zeros(8))
    verify(moved.subs(on_ray)-at_displacement)
    primitives = []
    for i in range(3):
        # Integrate (1-v) B''_ij(tu + v(w-tu)), whose Hessian is linear.
        correction = sp.zeros(8)
        for j in range(3):
            hessian = potential.diff(w[i],w[j])
            correction += (hessian/6+hessian.subs(on_ray)/3)*wedge[j]
        primitive = contraction[i]-correction
        verify(moved*primitive+primitive*moved-t**2*displacement.diff(u[i])*unit)
        primitives.append(primitive)
    for i in range(3):
        for j in range(i,3):
            hessian = potential.diff(w[i],w[j])
            homotopy = -sum((hessian.diff(x)*e/3 for x,e in zip(w,wedge)),sp.zeros(8))
            verify(primitives[i]*primitives[j]+primitives[j]*primitives[i]
                   +hessian.subs(on_ray)*unit-moved*homotopy-homotopy*moved)
    print("Curvature coefficients at orders 0 to 4:", [potential, 0, 0, displacement, 0])


def projective_curve():
    lam, X, Y, Z = sp.symbols("lambda X Y Z")
    x, y, z = (Z-lam*X+Y)/2, (Z-lam*X-Y)/2, X
    cubic = x**3+y**3+z**3-3*lam*x*y*z
    a = 4*(lam**3-1)/3
    quadratic, linear, constant = -3*lam**2, 2*lam, -sp.Rational(1,3)
    rhs = a*X**3+quadratic*X**2*Z+linear*X*Z**2+constant*Z**3
    assert sp.expand(4*cubic-3*(Z*Y**2-rhs)) == 0
    p = sp.factor(a*linear-quadratic**2/3)
    q = sp.factor(a**2*constant-a*quadratic*linear/3+2*quadratic**3/27)
    discriminant = sp.factor(-16*(4*p**3+27*q**2))
    assert sp.expand(discriminant-sp.Rational(4096,27)*(lam**3-1)**3) == 0
    b = sp.Symbol("b")
    j = sp.factor((-1728*64*p**3/discriminant).subs(lam, -2*b))
    assert sp.cancel(j-110592*b**3*(1-b**3)**3/(1+8*b**3)**3) == 0
    print("Projectivized unobstructed translations, j =", j)


def normal_disk():
    # Treat the Hessian as an arbitrary symmetric matrix, including its scale.
    xx, xy, xz, yy, yz, zz = sp.symbols("xx xy xz yy yz zz")
    hessian = sp.Matrix([[xx,xy,xz], [xy,yy,yz], [xz,yz,zz]])
    wedge, contraction = exterior_matrices()
    parity = sp.diag(*[(-1)**mask.bit_count() for mask in range(8)])
    rows = [sum((hessian[i,j]*wedge[j] for j in range(3)), sp.zeros(8))
            for i in range(3)]
    derivatives = [contraction[i]+rows[i]/3 for i in range(3)]
    generators = [contraction[i]-rows[i]/6 for i in range(3)]
    ordered = generators[0]*generators[1]*generators[2]
    for indices in permutations(range(3)):
        term = parity
        for i in indices:
            term = term*derivatives[i]
        sign = (-1)**sum(i>j for k,i in enumerate(indices) for j in indices[k+1:])
        assert sp.expand(sign*(term*ordered).trace()-hessian.det()/8) == 0
    x, y, z, b, rho = sp.symbols("x y z b rho")
    w = (x,y,z)
    W = rho*((x**3+y**3+z**3)/6+b*x*y*z)
    det = sp.hessian(W,w).det()
    difference = sum(v*W.diff(v) for v in w)
    assert sp.expand(det-rho**3*(1+8*b**3)*x*y*z+2*rho**2*b**2*difference) == 0
    fermat = (det/8)/sp.prod(W.diff(v) for v in w)
    assert sp.cancel(fermat.subs(b,0)) == 1/(x*y*z)
    print("All six signed normal supertraces equal det(W'')/8; Fermat unit residue is +1.")


def tangent_obstructions():
    x, y, z, b, t, s = sp.symbols("x y z b t s")
    w = (x,y,z)
    potential = (x**3+y**3+z**3)/6+b*x*y*z
    hessian = sp.hessian(potential,w)
    assert sp.expand(hessian.det()+6*b**2*potential-(1+8*b**3)*x*y*z) == 0
    # Moving from tu to tu-sv on the tangent line contributes -B(tu-sv).
    # The fixtures include a family of flexes and an ordinary smooth point.
    for coefficient, u, v, square, cube in (
        (b,(1,-1,0),(2*b,0,1),0,(1+8*b**3)/6),
        (-sp.Rational(5,6),(1,1,2),(1,-1,0),-8*t/3,0),
    ):
        values = dict(zip(w,(t*a-s*c for a,c in zip(u,v))))
        curvature = -potential.subs(b,coefficient).subs(values)
        assert sp.expand(curvature-s**2*square-s**3*cube) == 0
        if coefficient == b:
            flex_cubic = sp.expand(curvature).coeff(s,3)
            bulk_discriminant = (1-20*b**3-8*b**6)**2-64*(b**4-b)**3
            assert sp.expand(216*flex_cubic**3-bulk_discriminant) == 0
    # At the node of the triangular member the three closed generators have
    # a rank-two anticommutator on the generic fiber; the ray is its kernel.
    node = hessian.subs({x:1,y:1,z:1,b:-sp.Rational(1,2)})
    assert node.rank() == 2 and node.nullspace() == [sp.ones(3,1)]
    print("Tangent square at the ordinary point: -8*t/3; flex obstruction at order three: (1+8*b**3)/6.")


if __name__ == "__main__":
    brane_obstruction()
    projective_curve()
    normal_disk()
    tangent_obstructions()
