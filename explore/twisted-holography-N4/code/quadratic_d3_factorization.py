"""D3 Koszul matrix factorization and signed Kapustin--Li residue on the regular quadratic."""

from itertools import permutations
from math import factorial

from sympy import diag, diff, eye, hessian, simplify, symbols, zeros


def exterior_matrices():
    wedge = []
    for i in range(3):
        matrix = zeros(8)
        for mask in range(8):
            if mask & (1 << i) == 0:
                matrix[mask | (1 << i), mask] = (-1)**((mask % (1 << i)).bit_count())
        wedge.append(matrix)
    return wedge, [matrix.T for matrix in wedge]


def check_punctured_flux():
    z1, z2, w1, w2, w3, c, rho = symbols("z1 z2 w1 w2 w3 c rho")
    background = rho*(z1*w1 + z2*w2 + c*w3**2/2)
    # In order: scalar, two tangent vectors, tangent bivector. Derive the
    # internal differential from the primary Hamiltonian before restricting.
    def internal(field):
        _, x, y, area = field
        return (diff(background*x,z1) + diff(background*y,z2),
                -diff(background*area,z2), diff(background*area,z1), 0)

    def total(cochain):
        result = {}
        for chart in ((1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)):
            parity = (-1)**(len(chart)-1)
            field = cochain.get(chart, (0,)*4)
            output = [parity*entry for entry in internal(field)]
            for i in range(len(chart)):
                boundary = cochain.get(chart[:i]+chart[i+1:], (0,)*4)
                for component, entry in enumerate(boundary):
                    output[component] += (-1)**i * entry
            result[chart] = tuple(simplify(entry) for entry in output)
        return result

    for pole in (1, 2, 4):
        axial = 1/w3**pole
        flux = {(1,2,3): (0,0,0,axial/(w1*w2)),
                (1,3): (0,-rho*axial/w1,0,0),
                (2,3): (0,0,-rho*axial/w2,0),
                (3,): (-rho**2*axial,0,0,0)}
        assert all(value == (0,)*4 for value in total(flux).values())

        # Multiply by a normal profile which removes the pole on the axis.
        # Subtract a coboundary and recover the same scalar on EVERY chart.
        profile = w3**(pole+1)/2
        regular = simplify(profile*axial)
        primitive = {(1,2): (0,0,0,regular/(w1*w2)),
                     (1,): (0,rho*regular/w1,0,0),
                     (2,): (0,0,rho*regular/w2,0)}
        for chart, coboundary in total(primitive).items():
            field = flux.get(chart, (0,)*4)
            difference = tuple(simplify(profile*x-y) for x,y in zip(field,coboundary))
            expected = (-rho**2*regular,0,0,0) if len(chart)==1 else (0,)*4
            assert difference == expected

    laurent = 1/(w1*w2*w3)
    for normal in (w1,w2,w3):
        cubic = -diff(laurent, normal)
        assert simplify(internal((0,cubic/w1,0,0))[0]-rho*cubic) == 0
    assert simplify(diff(background*(-rho/(2*c)), w3) + rho**2*w3/2) == 0


def check():
    wedge, contraction = exterior_matrices()
    parity = diag(*[(-1)**j.bit_count() for j in range(8)])
    z1, z2, w1, w2, w3, c, rho = symbols("z1 z2 w1 w2 w3 c rho")
    coordinates = z1, z2, w1, w2, w3
    potential = rho*(z1*w1+z2*w2+c*w3**2/2)
    D = sum((w*op for w, op in zip(coordinates[2:], contraction)), zeros(8))
    D += rho*(z1*wedge[0]+z2*wedge[1]+c*w3*wedge[2]/2)
    assert simplify(D**2 - potential*eye(8)) == zeros(8)

    def differential(op, odd):
        return (D*op - (-1)**odd * op*D).applyfunc(simplify)

    eta = contraction[2] - rho*c*wedge[2]/2
    assert differential(eta, True) == zeros(8)
    assert eta**2 == -rho*c*eye(8)/2
    assert all(differential(op, True) == scalar*eye(8)
               for op,scalar in zip(wedge+contraction, (w1,w2,w3,rho*z1,rho*z2,rho*c*w3/2)))
    torsion = z1*contraction[1] - z2*contraction[0]
    assert differential(torsion, True) == torsion**2 == zeros(8)
    assert differential(contraction[0]*contraction[1], False) == rho*torsion

    jacobian = hessian(potential, coordinates).det()
    assert jacobian == rho**5*c
    # Equation (4) puts the boundary operator on the right, in this variable order.
    derivatives = [diff(D, x) for x in coordinates]
    exterior_power = zeros(8)
    for order in permutations(range(5)):
        inversions = sum(a>b for i,a in enumerate(order) for b in order[i+1:])
        term = eye(8)
        for i in order:
            term *= derivatives[i]
        exterior_power += (-1)**inversions * term
    boundary_state = simplify((parity*exterior_power*eta).trace()/factorial(5))
    assert boundary_state == rho**3*c
    assert (parity*exterior_power).trace() == 0
    assert simplify(boundary_state/jacobian) == 1/rho**2
    three_odd = simplify((parity*exterior_power*eta**3).trace()/(factorial(5)*jacobian))
    assert three_odd == -c/(2*rho)
    h = diag(1,1,-2)
    assert simplify(three_odd*2*(h**3).trace()) == 6*c/rho
    print("Hessian determinant:", jacobian, "odd boundary state:", boundary_state)
    print("Disk on one and three ordered odd unit insertions:", 1/rho**2, three_odd)


if __name__ == "__main__":
    check_punctured_flux()
    check()
