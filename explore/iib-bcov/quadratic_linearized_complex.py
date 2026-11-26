"""Exact checks of the local linearized descendant Hamiltonian on a quadratic."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod

from flint import fmpq, fmpq_mpoly_ctx
from sympy import Matrix


def rational(value):
    return fmpq(value.numerator, value.denominator)


def adjoint_series(through):
    """Coefficients of the logarithmic derivative of the constant-coefficient jet series."""
    series = [Fraction(1, 4**k*factorial(2*k+1)) for k in range(through+1)]
    log_derivative = [Fraction()]
    for k in range(1, through+1):
        log_derivative.append(k*series[k] - sum(log_derivative[i]*series[k-i]
                                                for i in range(1,k)))
    return log_derivative


def check_complex(hessian, through):
    ctx = fmpq_mpoly_ctx.get(('z1', 'z2', 'w1', 'w2', 'w3', 'q'))
    *coordinates, q = ctx.gens()
    zero = ctx.constant(0)
    determinant = int(hessian.det())
    adjugate = hessian.adjugate()
    coefficients = [(i,j,int(adjugate[i,j])) for i in range(5) for j in range(5) if adjugate[i,j]]
    background = sum(int(hessian[i,j])*coordinates[i]*coordinates[j]
                     for i in range(5) for j in range(5))/2 + coordinates[0] + 2
    slope = adjoint_series(through+1)

    def laplacian(poly):
        return sum((coef*poly.derivative(i).derivative(j) for i,j,coef in coefficients), zero)

    def apply_power(poly, k):
        for _ in range(k):
            poly = laplacian(poly)
        return poly

    def transfer(coefficient, poly):
        """Coefficient of x^coefficient in the transposed quadratic Euler series."""
        result = zero
        for j in range(coefficient//3+1):
            for k in range((coefficient-3*j)//2+1):
                m = coefficient-3*j-2*k
                denom = 24**j*factorial(j)*4**k*factorial(2*k+1)*factorial(m)
                result += determinant**j*q**(j+k)*apply_power(background**m*poly,k)/denom
        return result

    def directions(mask):
        return [(i, (-1)**((mask % (1<<i)).bit_count()), mask^(1<<i))
                for i in range(5) if mask & (1<<i)]

    def nonzero(rows):
        return {key: value for key,value in rows.items() if value}

    def differential(field, mode):
        output = defaultdict(lambda: zero)
        for (level,mask), poly in field.items():
            for i, sign, target in directions(mask):
                derivative = sign*poly.derivative(i)
                output[level+1,target] -= derivative
                if mode == 'original':
                    output[0,target] += sign*transfer(level+1,poly).derivative(i)
                else:
                    output[level,target] += sign*background.derivative(i)*poly
                if mode == 'transformed':
                    for k in range(1,level//2+1):
                        output[level-2*k,target] += (rational(2*determinant*slope[k])*q**k
                                                     *apply_power(derivative,k-1))
        return nonzero(output)

    def transform(field):
        output = defaultdict(lambda: zero)
        for (level,mask), poly in field.items():
            for k in range(level+1):
                output[level-k,mask] += transfer(k,poly)
        return nonzero(output)

    # Polynomials representing U(zeta^r). Each tuple is (level, Laplacian power,
    # loop order). The recurrence uses only constant-coefficient operators.
    basis = [{(0,0,0): Fraction(1)}]
    for _ in range(through+1):
        next_column = defaultdict(Fraction)
        for (level,power,h), coeff in basis[-1].items():
            next_column[level+1,power,h] += coeff
            if power:
                next_column[level,power-1,h] += 2*determinant*power*coeff
            for k in range(1,level//2+1):
                next_column[level-2*k,power+k-1,h+k] -= 2*determinant*slope[k]*coeff
        basis.append(nonzero(next_column))

    def triangular(field):
        output = defaultdict(lambda: zero)
        for (level,mask), poly in field.items():
            for (target,power,h), coeff in basis[level].items():
                output[target,mask] += rational(coeff)*q**h*apply_power(poly,power)
        return nonzero(output)

    poly = prod((1+coordinates[i])**(1+(i%3)) for i in range(5))
    for level in range(through+1):
        field = {(level,0b10101): poly}
        first = differential(field,'original')
        assert differential(first,'original') == {}
        assert transform(first) == differential(transform(field),'transformed')
        assert triangular(differential(field,'tree')) == differential(triangular(field),'transformed')
    print('Hessian rank',hessian.rank(),'determinant',determinant,
          'both intertwiners and square verified through level',through,flush=True)


def check_nonquadratic():
    # The level-zero test Y=z1*w1^2/2 times the oriented (w1,w3) bivector.
    # Keep all five variables and both background strengths in the Hessian.
    from sympy import diff, expand, hessian, symbols
    z1,z2,w1,w2,w3,rho,eps = symbols('z1 z2 w1 w2 w3 rho eps')
    coords = (z1,z2,w1,w2,w3)
    profile = w3**3/6 + w3**5
    matrix = hessian(rho*(z1*w1+z2*w2) + eps*profile, coords).adjugate()
    vector = diff(z1*w1**2/2,w1)
    primary = sum(matrix[i,j]*diff(vector,coords[i],coords[j]) for i in range(5) for j in range(5))
    assert expand(-diff(primary,w3)/24 + rho**3*eps*diff(profile,w3,3)/12) == 0


def check_planar_primitives():
    from sympy import diff, exp, factorial as fac, hessian, simplify, symbols
    z1,z2,x,y,z,t = symbols('z1 z2 x y z t')
    background = (x**2+y**2+z**2)/2

    def shifted(level):
        return sum(t**(level-j)*(-background)**j/fac(j) for j in range(level+1))

    # On tangent divergence, D(ζ^r Y) has tree coefficients -ζ^(r+1)
    # and B^(r+1)/(r+1)!. Test both primitives before any change of variables.
    for level in range(1,5):
        tree = sum((-background)**j/fac(j)*background**(level-j)/fac(level-j)
                   for j in range(level))
        assert simplify(t*shifted(level-1)-tree-shifted(level)) == 0
    scalar, vector = exp(z1), exp(z2)
    assert diff(scalar,z1)==scalar and diff(vector,z1)==0
    # In the one-loop action the vector changes the (2,2) Hessian entry.
    coords = (z1,z2,x,y,z)
    strength = symbols('strength')
    cofactor = hessian(background+strength*vector,coords).adjugate().diff(strength).subs(strength,0)
    gradient = sum(diff(cofactor[i,j]*scalar, coords[i], coords[j])
                   for i in range(5) for j in range(5))
    assert simplify(diff(gradient,z1)-exp(z1+z2)) == 0


def check_primary_threshold():
    from itertools import product
    from bulk_eight import SHAPES, potentials
    from longitudinal_theta import primary_coefficient

    unit = [tuple(int(i==j) for i in range(4)) for j in range(4)]
    scalars = unit[:3]
    vectors = (unit[3], (1,-1,0,0))
    polarization = (1,1,0,0)

    def evaluate(momentum):
        return sum(x*y for x,y in zip(momentum,polarization))

    def quartic(waves, vector):
        amplitudes = []
        for first, second in product((-1,1), repeat=2):
            normal = [(0,0,0,0,first),(0,0,0,0,second)]
            fields = [(*wave,0) for wave in waves]+normal
            amplitudes.append(primary_coefficient(fields,(*vector,0),(*polarization,0)))
        # Each cosh has two coefficients 1/2. The repeated scalar has 1/2!.
        actual = sum(amplitudes)/8
        expected = Fraction(int(Matrix([*waves,vector]).det())**4
                            *sum(map(evaluate,waves)),62208)
        assert actual == expected
        return actual

    # The first vector contracted has a left odd parameter. The second
    # contraction passes two such parameters. Hold their order fixed.
    inner_four = []
    inner_two = []
    total = sum(map(evaluate,scalars))
    for sign, first, second in ((-1,*vectors),(1,*reversed(vectors))):
        inner_four.append(sign*total*quartic(scalars,first))
        for i,scalar in enumerate(scalars):
            fused = list(scalars)
            fused[i] = tuple(x+y for x,y in zip(first,scalar))
            inner_two.append(sign*evaluate(scalar)*quartic(fused,second))
    assert sum(inner_four) == -Fraction(4,62208)
    assert sum(inner_two) == Fraction(28,62208)
    # Integrate the four theta sectors on the original eight conserved waves.
    # The antisymmetrized vector slots remove the common grounded potential.
    theta = []
    for shape in SHAPES:
        normal_signs = []
        for first, second in product((-1,1), repeat=2):
            momenta = [(*wave,0) for wave in (*scalars,*vectors)]
            momenta += [(0,0,0,0,first),(0,0,0,0,second)]
            momenta.append(tuple(-sum(col) for col in zip(*momenta)))
            slots = potentials(momenta,check=True,shapes=(shape,))
            difference = tuple(x-y for x,y in zip(slots[3],slots[4]))
            normal_signs.append(evaluate(difference[:4]))
        assert len(set(normal_signs)) == 1
        theta.append(total*sum(normal_signs)/8)
    assert theta == [Fraction(131,1632960),Fraction(67,408240),Fraction(1,12096),Fraction(1,17010)]
    assert sum(inner_four+inner_two) == sum(theta)
    print('Four-dimensional primary Ward compositions:',sum(inner_four),sum(inner_two),
          'combined:',sum(inner_four+inner_two),flush=True)
    print('Compensating theta shapes:',*theta,flush=True)


if __name__ == '__main__':
    for axial in (0,2):
        matrix = Matrix([[0,0,1,0,0],[0,0,0,1,0],[1,0,0,0,1],
                         [0,1,0,0,0],[0,0,1,0,axial]])
        check_complex(matrix, 6)
    check_nonquadratic()
    check_planar_primitives()
    check_primary_threshold()
