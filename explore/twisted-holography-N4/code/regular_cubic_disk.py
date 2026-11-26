"""Regular D3 source, in units of its induced tangent-coordinate commutator."""

from itertools import combinations_with_replacement, product

from sympy import I, Matrix, Rational, diag, eye, simplify, symbols, zeros

from boundary_odd_sectors import wedge


# Forms and odd amplitudes are ordered dz1, dz2, theta1, theta2, theta3, eta1, eta2.
def multiply(left, right):
    result = {}
    for a, x in left.items():
        for b, y in right.items():
            if a & b:
                continue
            inversions = sum((a >> (i + 1)).bit_count() for i in range(7) if b & (1 << i))
            key = a | b
            result[key] = result.get(key, zeros(x.rows)) + (-1)**inversions * x*y
    return result


def word(*factors):
    result = factors[0]
    for factor in factors[1:]:
        result = multiply(result, factor)
    return result


def combine(*polys):
    result = {}
    for poly in polys:
        for key, value in poly.items():
            result[key] = result.get(key, zeros(value.rows)) + value
    return result


def open_product(a, b, coordinates):
    z1, z2 = coordinates
    def derivative(poly, coordinate):
        return {key: entry.diff(coordinate) for key, entry in poly.items()}
    ordinary = multiply(a, b)
    first = multiply(derivative(a, z1), derivative(b, z2))
    second = multiply(derivative(a, z2), derivative(b, z1))
    unit = eye(next(iter(a.values())).rows)
    correction = combine(first, {key: -value for key, value in second.items()})
    return combine(ordinary, multiply({24: unit/2}, correction))


def descendant_graphs():
    cases = 0
    for vertices in range(1, 5):
        possible_edges = list(combinations_with_replacement(range(vertices), 2))
        for count in range(vertices-1, vertices+2):
            for edges in combinations_with_replacement(possible_edges, count):
                connected = {0}
                for _ in range(vertices):
                    for i, j in edges:
                        if i in connected or j in connected:
                            connected.update((i, j))
                if len(connected) < vertices:
                    continue
                value = {0: 1}
                for i, j in edges:
                    delta = {1 << i: 1, 1 << j: -1} if i != j else {}
                    value = wedge(value, delta)
                assert bool(value) == (count == vertices-1)
                cases += 1
    for genus in range(1, 15):
        for shift in range(genus+1):
            for axis in range(3):
                powers = tuple(int(i==axis) + (2*shift if i==2 else 0) for i in range(3))
                for differentiated in range(3):
                    if not powers[differentiated]:
                        continue
                    odd_derivatives = tuple(p - int(i==differentiated) for i, p in enumerate(powers))
                    assert min(odd_derivatives) == 0
                    arity = 2 + sum(odd_derivatives) + 2*(3*genus-shift) + 1
                    assert arity == 6*genus+3
    print("Undifferentiated odd delta on connected open graphs:", cases)


def mixed_edge_rank():
    # Two bulk vertices joined by an ordinary line. Their other closed edges
    # end on two of three D3 vertices; two open lines join the three D3 vertices.
    sizes = [5, 5, 2, 2, 2]
    positions = [(1, 2, 1, 3, 5), (3, 1, 4, 1, -2), (4, -2), (-1, 3), (7, 2)]
    edges = [(0, 1, 5), (0, 2, 5), (1, 3, 5), (2, 4, 2), (3, 4, 2)]
    euler = Matrix([component for point in positions for component in point])
    rows = []
    for i, j, dimension in edges:
        incidence = zeros(dimension, sum(sizes))
        for vertex, sign in [(i, 1), (j, -1)]:
            for slot in range(min(sizes[vertex], dimension)):
                incidence[slot, sum(sizes[:vertex])+slot] = sign
        difference = incidence*euler
        rows.extend(v.T*incidence for v in difference.T.nullspace())
    forms = Matrix.vstack(*rows)
    assert forms.rows == sum(sizes)-2 and forms.rank() == forms.rows-1
    translations = [Matrix([int(j==axis) for size in sizes for j in range(size)]) for axis in (0,1)]
    assert all(forms*v == zeros(forms.rows,1) for v in [*translations, euler])
    assert Matrix.hstack(*translations, euler).rank() == 3
    print("Ordinary separated mixed graph: required form degree", forms.rows, "rank", forms.rank())


def coefficients(colors, momenta):
    d, e, f = colors
    unit = eye(d.rows)
    c = {32: I*d, 64: I*e, 4: f}
    dz1, dz2 = {1: unit}, {2: unit}
    # The derivative slots are ordinary coordinate derivatives in both symbols.
    derivatives = [{key: value*p[j] for (key, value), p in zip(c.items(), momenta)}
                   for j in range(2)]
    # Unit coordinate commutator: antisymmetric bidifferential coefficient 1/2.
    first = word({24: unit}, c, *derivatives)
    second = word({24: unit}, c, *reversed(derivatives))
    cubic = (first[124] - second[124]).trace()/6
    # dtheta1 is even. Extract its coefficient by the three positions of the scalar.
    odd_derivatives = [{key: value for key, value in deriv.items() if key in (32, 64)}
                       for deriv in derivatives]
    tangent = combine(word(dz1, odd_derivatives[0]), word(dz2, odd_derivatives[1]))
    scalar = {0: f}
    graviton = sum(word(*args)[99].trace() for args in
                   [(scalar, tangent, tangent), (tangent, scalar, tangent), (tangent, tangent, scalar)])
    return cubic, graviton


def check():
    h = diag(1, 1, -2)
    t1, t2, t3 = Matrix([[0, 1], [1, 0]])/2, Matrix([[0, -I], [I, 0]])/2, diag(1, -1)/2
    momenta = ((1, 0), (0, 1), (-1, -1))
    assert coefficients((h, h, h), momenta) == (6, -36)
    for colors in product((t1, t2, t3), repeat=3):
        assert coefficients(colors, momenta) == (0, 0)
    # A central external leg detects the rank-two matrix identity.
    assert coefficients((t1, t1, eye(2)), momenta) == (-Rational(1, 2), 3)
    e12, e21 = zeros(3), zeros(3)
    e12[0, 1] = e21[1, 0] = 1
    for p, q in [((2, -1), (1, 2)), ((-1, 2), (3, 1)), ((1, 2), (-2, -4))]:
        recoil = tuple(-x-y for x, y in zip(p, q))
        for colors in [(h, h, h), (e12, e21, h), (e21, e12, h)]:
            d, e, f = colors
            invariant = (f*(d*e+e*d)).trace()
            area = p[0]*q[1]-p[1]*q[0]
            assert coefficients(colors, (p, q, recoil)) == (-area*invariant/2, 3*area*invariant)
    print("Unit coordinate commutator and trace HCS cubic versus protected graviton:",
          coefficients((h, h, h), momenta))
    print("All 27 traceless SU(2) polarizations vanish; one central leg:",
          coefficients((t1, t1, eye(2)), momenta))

    z1, z2 = coordinates = symbols("z1 z2")
    for color in (t3, h):
        a, b = {0: z1*color}, {0: z2*color}
        assert open_product(a, b, coordinates)[24] - open_product(b, a, coordinates)[24] == color**2
    assert t3**2 == eye(2)/4 and h**2 == 2*eye(3)-h
    a = {0: z1*z2*h, 4: z1**2*e12}
    b = {0: z1*e21, 8: z2*h}
    c = {0: z2**2*e12, 16: z1*h}
    left = open_product(open_product(a, b, coordinates), c, coordinates)
    right = open_product(a, open_product(b, c, coordinates), coordinates)
    for key in set(left) | set(right):
        assert simplify(left.get(key, zeros(3)) - right.get(key, zeros(3))) == zeros(3)
    print("Cartan commutators at ranks two and three: I/4 and 2I-H.")
    descendant_graphs()
    mixed_edge_rank()


if __name__ == "__main__":
    check()
