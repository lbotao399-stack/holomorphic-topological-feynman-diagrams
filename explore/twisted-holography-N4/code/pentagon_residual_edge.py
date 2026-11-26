"""One regular residual edge adjacent to one puncture propagator at the flux primitive.

The heat widths are 4t. Angular moments use the normalized complex unit sphere.
No common BCOV edge, cutoff, or disk normalization is included.
"""

from itertools import permutations
from math import factorial, prod

import sympy as sp


def antisymmetric_sign(indices):
    return (-1)**sum(i > j for n, i in enumerate(indices) for j in indices[n+1:])


def odd_and_dolbeault():
    # Traverse the cycle from the root through p, v, s, q and back to the root.
    # The first and last currents differentiate the two endpoint tangent deltas.
    p = sp.Matrix([1, 2, 1])
    q = sp.Matrix([2, -1, 4])
    nodes = [p, sp.Matrix([-1, 2, 2]), sp.Matrix([3, 1, 0]), q]
    normal = [nodes[e+1]-nodes[e] for e in range(3)]
    currents = [list(p[:2])+[0]*3, *([0, 0]+list(v) for v in normal),
                list(q[:2])+[0]*3]
    edge_volume = sp.zeros(25)
    for e, covector in enumerate(currents):
        for i in range(5):
            for site, sign in ((e, -1), ((e+1) % 5, 1)):
                if site:
                    edge_volume[5*e+i, 5*(site-1)+i] = sign
            edge_volume[5*e+i, 20+e] = covector[i]
    area = p[0]*q[1]-p[1]*q[0]
    # The five auxiliary columns insert the contractions into five degree-five
    # volumes. Each resulting degree-four block has even exterior parity.
    oriented = -area*sp.Matrix.hstack(*normal).det()
    assert edge_volume.det() == sp.Matrix(currents).det() == oriented != 0

    # First row is q dot d(conjugate q). The next three five-form rows use
    # the three conjugate displacements as contraction vectors. The columns
    # integrate tangent v,s, then normal v,s,q, then the three auxiliaries.
    z = [sp.zeros(2, 1), sp.Matrix([2, 1]), sp.Matrix([-1, 3]), sp.zeros(2, 1)]
    variables = [[None]*5, [0, 1, 4, 5, 6], [2, 3, 7, 8, 9],
                 [None, None, 10, 11, 12]]
    pulled = sp.zeros(16)
    for i in range(3):
        pulled[0, 10+i] = q[i]
    for e in range(3):
        diff = list(z[e+1]-z[e])+list(normal[e])
        for i in range(5):
            row = 1+5*e+i
            for site, sign in ((e, -1), (e+1, 1)):
                column = variables[site][i]
                if column is not None:
                    pulled[row, column] = sign
            pulled[row, 13+e] = diff[i]
    reduced = sp.Matrix.hstack(z[1], z[2]).det()*q.dot(q-p)
    assert pulled.det() == reduced != 0


def gaussian_contractions():
    # Wick directly on the independent centered Gaussian in each coordinate.
    # Each letter is (site, coordinate); paired opposite types contribute C_ij.
    c11, c12, c22 = sp.symbols("c11 c12 c22")
    covariance = [[c11, c12], [c12, c22]]
    determinant = c11*c22-c12**2

    def wick(anti, hol):
        return sum(prod(covariance[i][j] if a == b else 0
                        for (i, a), (j, b) in zip(anti, perm))
                   for perm in permutations(hol))

    for first, second in ((0, 1), (1, 0), (0, 0), (1, 1)):
        value = sum(antisymmetric_sign(pair) *
                    wick([(0, pair[0]), (1, pair[1])], [(0, first), (1, second)])
                    for pair in permutations(range(2)))
        assert sp.expand(value-(second-first)*determinant) == 0
    direction = sp.symbols("d1:4")
    normal = sum(antisymmetric_sign(perm)*direction[perm[2]] *
                 (wick([(0, perm[0]), (1, perm[1])], [(0, 0), (1, 1)]) -
                  wick([(0, perm[0]), (1, perm[1])], [(0, 1), (1, 0)]))
                 for perm in permutations(range(3)))
    assert sp.expand(normal-2*direction[2]*determinant) == 0

    t1, t2, t3 = sp.symbols("t1 t2 t3", positive=True)
    total = t1+t2+t3
    covariance = 4*sp.Matrix([[t1*(t2+t3), t1*t3],
                              [t1*t3, (t1+t2)*t3]])/total
    assert sp.factor(covariance.det()-16*t1*t2*t3/total) == 0
    # Ordered odd and Dolbeault tensors have the respective signs - and +.
    # This omits the orientation of both endpoint delta derivatives, BV edges,
    # annular forms and the final Feynman prefactor.


def sphere_integral(poly, coordinates):
    result = sp.S.Zero
    for exponents, coeff in sp.Poly(poly, *coordinates).terms():
        for start in range(0, len(coordinates), 6):
            hol, anti = exponents[start:start+3], exponents[start+3:start+6]
            if hol != anti:
                break
            coeff *= sp.Rational(2*prod(factorial(e) for e in hol), factorial(2+sum(hol)))
        else:
            result += coeff
    return result


def angular_moments():
    x = sp.symbols('x1:4')
    rx = sp.symbols('r1:4')
    p, bp, q, bq = [sp.symbols(letter+'1:4') for letter in ('p', 'P', 'q', 'Q')]
    sphere = p+bp+q+bq
    bracket = p[0]*q[1]-p[1]*q[0]
    exterior = bracket*(1-sum(q[i]*bp[i] for i in range(3)))*(bq[2]-bp[2])
    values = []
    for i, j in ((0, 1), (1, 0)):
        # H=x_j r_3 has squared norm 1/12; the regular dual of x_i is 3 r_i U.
        root = sphere_integral(antisymmetric_sign((2, j, i)) *
                               x[2]*x[i]*rx[2]*rx[i]/6, x+rx)
        far = sphere_integral(3*bp[i]*12*bq[j]*q[2]*exterior, sphere)
        values.append((root, far))
    assert values == [(-sp.Rational(1, 72), 1), (sp.Rational(1, 72), -1)]
    assert sum(root*far for root, far in values) == -sp.Rational(1, 36)

    # In the T->infinity leading coefficient there is only one holomorphic p.
    # Balance of p phases excludes |alpha|>1. The alpha=0 root curl vanishes
    # by angular integration by parts; alpha=e3 has the same vanishing curl.
    # For alpha=e1,e2 the remaining q_j conjugate(q_3) is already trace free.
    # Hence the two entries above exhaust that coefficient for this ordering.
    print("Root and normalized remote unit-sphere moments:", values)


def harmonic_sum():
    x = sp.symbols('x1:4')
    rx = sp.symbols('r1:4')
    for first in range(5):
        for second in range(5-first):
            for third in range(5-first-second):
                monomial = x[0]**first*x[1]**second*x[2]**third
                curl = x[2]*(rx[1]*sp.diff(monomial, x[0])-rx[0]*sp.diff(monomial, x[1]))
                assert sp.expand(sum(sp.diff(curl,x[i],rx[i]) for i in range(3))) == 0
    for k in range(1, 7):
        harmonic = x[1]*rx[0]**(k-1)*rx[2]
        direct = -x[2]*x[0]**k*rx[0]*sp.diff(harmonic,x[1])
        moved = -harmonic*k*x[2]*rx[1]*x[0]**(k-1)
        assert sphere_integral(direct-moved, x+rx) == 0
    # Check the U(3) average by integrating the third row of a unitary matrix
    # on the complex unit sphere, including examples with a complex dot product.
    for p, q in [
        (sp.Matrix([1, 2, 1]), sp.Matrix([2, -1, 4])),
        (sp.Matrix([1, sp.I, 1]), sp.Matrix([2, -1, 1+sp.I])),
    ]:
        row = sp.Matrix(x)
        dual = sp.Matrix(rx)
        v = p.cross(q).conjugate()
        d = q-p
        value = (dual.dot(v)*row.dot(v.conjugate())*dual.dot(q)*row.dot(d.conjugate()))
        a, b, cross = p.dot(p.conjugate()), q.dot(q.conjugate()), q.dot(p.conjugate())
        predicted = (a*b-cross*cross.conjugate())*(b-cross)/12
        assert sp.simplify(sphere_integral(value, x+rx)-predicted) == 0
    # Unit radii. At angular expansion order m the only surviving terms
    # in (z+conj(z))**m/(1-z) have matching exponents. |z|² has density 2(1-s).
    result = []
    for m in range(3):
        val = 0
        for j in range(m//2+1):
            power = m-j
            val += sp.binomial(m,j)*sp.Rational(4, (power+1)*(power+2)*(power+3))
        result.append(-val/(24*factorial(m)))
    assert result == [-sp.Rational(1,36), -sp.Rational(1,144), -sp.Rational(1,120)]
    print("Ordered equal-shell angular series, prior to the common radial Gaussian:", result)


def separated_shell_times():
    """After the three time integrals the ordered angular kernel is a ball Poisson integral."""
    t, s = sp.symbols('t s', real=True)
    f = (t*t-t*s)**2/(1-t*s)**3
    assert sp.cancel(f.subs(s,t)) == 0
    # On S^7, the first complex coordinate has density 3(1-|s|²)^2/pi.
    # Check the harmonic extension order by order without using the Poisson identity.
    b = sp.Symbol('b')
    order = 10
    u = sp.series((1-t*t)*f/(1+t*t-t*(s+b))**4, t, 0, order+1).removeO()
    for k in range(order+1):
        poly = sp.Poly(u.expand().coeff(t,k),s,b)
        mean = sum(coeff*sp.Rational(6*factorial(i),factorial(i+3))
                   for (i,j),coeff in poly.terms() if i==j)
        assert mean == 0
    # Equal shells. Integrate exactly across the slice Im(s) at 1-Re(s)=x.
    x, y, z = sp.symbols('x y z', positive=True)
    integrand = x*(z*z-y*y)**2/(x*x+y*y)
    primitive = x*(y**3/3-(2*z*z+x*x)*y)+(z*z+x*x)**2*sp.atan(y/x)
    assert sp.cancel(sp.diff(primitive,y)-integrand) == 0
    bound = sp.sqrt(2*x-x*x)
    section = 2*primitive.subs(y,z).subs(z,bound)
    assert sp.limit(section/x**2,x,0,dir='+') == 4*sp.pi
    # Integrate the one-sided radial scaling limit. Shift v by w²/2 first.
    angle = sp.symbols('angle', real=True)
    radial = sp.integrate(4*x*x*sp.exp(-x/2),(x,0,sp.oo))
    half_plane = (sp.integrate(sp.exp(-x*x/4),(x,0,sp.oo))
                  *sp.integrate(sp.cos(5*angle),(angle,-sp.pi/2,sp.pi/2)))
    assert radial == 64 and half_plane == 2*sp.sqrt(sp.pi)/5
    assert sp.simplify(-2*256*radial*half_plane/(6*sp.pi)
                       +32768/(15*sp.sqrt(sp.pi))) == 0
    print('Disjoint ordered shells: zero time-integrated tensor; equal-shell leading slice:',4*sp.pi)


if __name__ == "__main__":
    odd_and_dolbeault()
    gaussian_contractions()
    angular_moments()
    harmonic_sum()
    print("The ordered 25- and 16-dimensional determinants and both Gaussian contractions agree.")
    separated_shell_times()
