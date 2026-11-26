"""Exact BRST tests on two ghost jets and one scalar in the D3 boundary theory."""

from itertools import product

from sympy import I, Matrix, Rational, diag, eye, integrate, symbols, zeros


def add(p, q):
    return tuple(x + y for x, y in zip(p, q))


def bracket(p, q):
    return p[0] * q[1] - p[1] * q[0]


def handle(p, q):
    return bracket(p, q)**2 * (
        Rational(1, 2) * (p[0]**2 + q[0]**2) - Rational(4, 3) * p[0] * q[0]
    )


def completed_handle(p, q, scalar):
    # Polarized coefficient divided by two; the primitive symbol is handle/g.
    return (handle(p, add(q, scalar)) + handle(q, add(p, scalar))
            - handle(add(p, q), scalar)) / 2


def commutator(x, y):
    return x*y - y*x


def jordan(x, y):
    value = x*y + y*x
    return value - eye(x.rows)*value.trace()/x.rows


def current_variation(symbol, color, ghosts, scalar):
    """Matrix in the trace against the scalar, before the common BRST factor."""
    result = zeros(ghosts[0][1].rows)
    for k in range(3):
        (p, x), (q, y), (r, z) = ghosts[k:] + ghosts[:k]
        result += (symbol(add(p, q), r, scalar) * color(commutator(x, y), z)
                   + symbol(p, q, add(r, scalar)) * commutator(color(x, y), z))
    return result


def check():
    zero = (0, 0)
    p, q, r = (1, 0), (2, 0), (0, 1)
    # Source formula Q0 c^3 = -g c^1 c^2. Relative waves have value zero
    # at the point; polarizing the original two ghost slots supplies a factor 2.
    def relative_product(left, right, third):
        return handle(add(left, right), third) - handle(left, third) - handle(right, third)

    factors = relative_product(p, r, q), relative_product(q, r, p)
    assert factors == (Rational(-26, 3), Rational(-2, 3))
    assert 2 * (factors[0] - factors[1]) == -16
    print("Relative two-loop coefficient in units of the component prefactor times g:", -16)

    # Independently include the variation of the scalar on pure waves of
    # conserved total momentum. The Pauli normalization is tr(t_A t_B)=δ_AB/2.
    t1, t2 = Matrix([[0, 1], [1, 0]])/2, Matrix([[0, -I], [I, 0]])/2
    ghosts = [(p, t1), (q, t1), (r, t2)]
    conserved = (-3, -1)
    bare = current_variation(lambda p, q, s: handle(p, q), commutator, ghosts, conserved)
    assert (bare*t2).trace() == -4
    assert current_variation(completed_handle, commutator, ghosts, conserved) == zeros(2)

    p_sym = symbols("p1 p2")
    q_sym = symbols("q1 q2")
    assert (completed_handle(p_sym, q_sym, zero) - handle(p_sym, q_sym)).expand() == 0

    # Two independent tests of the cubic color identities in sl_3.
    h = diag(1, -1, 0)
    e12, e21, e13 = zeros(3), zeros(3), zeros(3)
    e12[0, 1] = e21[1, 0] = e13[0, 2] = 1

    def tensors(colors):
        return [jordan(commutator(colors[i], colors[(i+1) % 3]), colors[(i+2) % 3])
                for i in range(3)]

    first, second, third = tensors([h, h, e13])
    assert first == zeros(3) and second == -third and second != zeros(3)
    first, second, third = tensors([h, e12, e21])
    assert first == second == third and (first*diag(1, 1, -2)).trace() == 4

    def disk(p, q, s):
        total = add(add(p, q), s)
        return bracket(p, q)*(1 + total[0]*total[1]**2)

    def quintic(p, q, s):
        return bracket(p, q)**5

    ghosts = [((1, 0), h), ((1, 0), e12), ((0, 1), e21)]
    defect = current_variation(quintic, jordan, ghosts, zero)
    assert (defect*diag(1, 1, -2)).trace() == 120
    for scalar in [zero, (-2, -1), (2, 3)]:
        assert current_variation(disk, jordan, ghosts, scalar) == zeros(3)
        assert current_variation(completed_handle, commutator, ghosts, scalar) == zeros(3)

    for p, q in product([(1, 0), (2, -1), (-1, 3)], repeat=2):
        assert completed_handle(p, q, zero) == handle(p, q)
    # Same total momentum, two areas. No polynomial of the total times the
    # linear area can reproduce the quintic even on the constant scalar.
    pairs = [((1, 0), (0, 1)), ((2, 0), (-1, 1))]
    assert [add(*pair) for pair in pairs] == [(1, 1)]*2
    assert [quintic(*pair, zero)//bracket(*pair) for pair in pairs] == [1, 16]
    print("Traceless symmetric-color quintic variation (unit matrix convention):", 120)
    print("Constant scalar quintic/linear ratios at the same total momentum:", (1, 16))

    # One-loop Konishi primitive with two derivatives on only the scalar.
    # Sum both placements of the differentiated letter on the original triangle.
    a, b, p_u, q_u = symbols("a b p_u q_u")
    displacements = ((1-a)*p_u + b*q_u, a*p_u + (1-b)*q_u)
    weight = integrate(sum(x**2 for x in displacements), (a, 0, 1-b), (b, 0, 1))
    assert (weight - (p_u**2 + q_u**2)/3 - p_u*q_u/2).expand() == 0
    assert [weight.subs({p_u: x, q_u: 1-x}) for x in (0, 2)] == [Rational(1, 3), Rational(2, 3)]
    print("One-loop split-derivative weights at fixed total:", (Rational(1, 3), Rational(2, 3)))


if __name__ == "__main__":
    check()
