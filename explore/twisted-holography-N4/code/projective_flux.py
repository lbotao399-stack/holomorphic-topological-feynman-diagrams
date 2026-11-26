"""Special conformal action in the punctured BCOV complex, forms ordered first."""
from itertools import product
from sympy import QQ
from sympy.polys.rings import ring


R, z1, z2, w1, w2, w3, v1, v2, v3, u, N = ring("z1,z2,w1,w2,w3,v1,v2,v3,u,N", QQ)
z = (z1, z2)
w = (w1, w2, w3)
v = (v1, v2, v3)  # complex conjugate normal coordinates
relation = u * sum(x * y for x, y in zip(w, v)) - 1
x = z + w


def add(*terms):
    out = {}
    for term in terms:
        for mask, value in term.items():
            out[mask] = out.get(mask, R.zero) + value
    return {mask: value for mask, value in out.items() if value}


def scale(a, value):
    return {mask: value * coeff for mask, coeff in a.items()}


def wedge(a, b):
    terms = []
    for i, ai in a.items():
        for j, bj in b.items():
            if not i & j:
                crossings = sum((j & ((1 << n) - 1)).bit_count() for n in range(8) if i & (1 << n))
                terms.append({i | j: (-1)**crossings * ai * bj})
    return add(*terms)


def contract(a, n):
    bit = 1 << n
    return {m ^ bit: (-1)**((m & (bit - 1)).bit_count()) * c for m, c in a.items() if m & bit}


def derivative(a, coordinate):
    conjugate = dict(zip(w + v, v + w)).get(coordinate, R.zero)
    return {m: c.diff(coordinate) - u**2 * conjugate * c.diff(u) for m, c in a.items()}


def bar(a):
    return add(*(wedge(anti[i], derivative(a, v[i])) for i in range(3)))


def divergence(a):
    return add(*(derivative(contract(a, j + 3), coordinate) for j, coordinate in enumerate(x)))


def bracket(a, b):
    # Same suspended product rule as code/holography/bulk_x.py::schouten.
    signed = {m: (-1)**m.bit_count() * c for m, c in a.items()}
    return add(*(add(wedge(contract(a, j + 3), derivative(b, coordinate)),
                      wedge(derivative(signed, coordinate), contract(b, j + 3)))
                 for j, coordinate in enumerate(x)))


def check_zero(a):
    assert all(not c.rem([relation]) for c in a.values())


def epsilon(i, j, k):
    return (j - i) * (k - i) * (k - j) // 2


anti = [{1 << i: R.one} for i in range(3)]
xi = [{1 << (i + 3): R.one} for i in range(5)]
b = []
U = {}
for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    b.append(add(scale(anti[k], u**2 * v[j]), scale(anti[j], -u**2 * v[k])))
    U = add(U, scale(wedge(anti[j], anti[k]), 2 * u**3 * v[i]))
P = wedge(xi[0], xi[1])
F = wedge(U, P)


def differential(a):
    return add(bar(a), scale(bracket(F, a), N))


def hamiltonian(polynomial):
    return bracket(P, {0: polynomial})


def lift(s, polynomial):
    return add({0: w[s] * polynomial}, scale(wedge(b[s], hamiltonian(polynomial)), -N))


def check():
    euler = add(*(scale(e, c if j < 2 else -c) for j, (e, c) in enumerate(zip(xi, x))))
    normal = add(*(wedge(b[i], xi[2 + i]) for i in range(3)))
    generators = []
    for coordinate in z:
        corrected = add(scale(euler, coordinate), scale(wedge(normal, hamiltonian(coordinate)), -N))
        check_zero(differential(corrected))
        check_zero(divergence(corrected))
        generators.append(corrected)

    primitive = add(*(scale(wedge(xi[2 + (k + 1) % 3], xi[2 + (k + 2) % 3]), 2 * N * u * v[k])
                      for k in range(3)))
    check_zero(divergence(primitive))
    check_zero(add(bracket(*generators), scale(differential(primitive), -1)))

    check_zero(bracket(F, primitive))
    # Psi sends +z_b d/dtheta_s to -z_b w_s.  The suspended BV
    # bracket is ordered K_a first; reversing it gives the shifted Lie convention.
    for a, s in product(range(2), range(3)):
        other = 1-a
        odd = scale(lift(s, z[other]), -1)
        check_zero(differential(odd))
        check_zero(divergence(odd))
        mixed_primitive = add(*(scale(xi[2+i], (1-2*a)*N*epsilon(i, s, k)*u*v[k])
                                for i, k in product(range(3), repeat=2)))
        check_zero(bracket(F, mixed_primitive))
        check_zero(divergence(mixed_primitive))
        check_zero(add(bracket(generators[a], odd), scale(bar(mixed_primitive), -1)))

    for s, degree in product(range(3), range(5)):
        H = z2**degree
        source = lift(s, H)
        check_zero(differential(source))
        check_zero(divergence(source))
        image = bracket(generators[0], source)
        scalar_image = lift(s, z1 * (z2 * H.diff(z2) - H))
        primitive = {}
        mixing = {}
        for i, k in product(range(3), repeat=2):
            sign = epsilon(i, s, k)
            if sign:
                primitive = add(primitive, scale(xi[2 + i], -N * sign * u * v[k] * H.diff(z2)))
                polarized = wedge(derivative(U, w[k]), wedge(xi[2 + i], xi[0]))
                mixing = add(mixing, scale(polarized, -N**2 * sign * H.diff(z2).diff(z2) / 2))
        check_zero(divergence(primitive))
        check_zero(differential(mixing))
        check_zero(divergence(mixing))
        check_zero(add(image, scale(scalar_image, -1), scale(differential(primitive), -1), scale(mixing, -1)))
    print("All seven projective symmetry brackets have flux-annihilated primitives; 15 scalar actions reproduce the signed Laurent mixing.")


if __name__ == "__main__":
    check()
