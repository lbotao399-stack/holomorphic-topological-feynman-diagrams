"""Four field orientations at fixed, divergence-free external vector polarizations."""
from fractions import Fraction

from sympy import QQ
from sympy.polys.rings import ring

from check_theta166_flow import SWAPPED, mean, reference_setup


def transverse_polynomial(b_axes):
    """B2 through B5 carry the specified vector basis labels, counted from zero."""
    nv, edges, paths, lam, external, flows = reference_setup(SWAPPED)
    from v_integrand import grouped_integrand
    from vcommon import D_theta, kirchhoff

    basis = [[int(i == j) for j in range(5)] for i in range(5)]
    for vertex, axis in zip((2, 3, 4, 5), b_axes):
        external[vertex]["u"] = basis[axis]
    external[6]["u"] = basis[1]
    for vertex in range(2, 11):
        assert sum(v * p for v, p in zip(external[vertex]["u"], lam[vertex])) == 0
    assert all(sum(coords) == 0 for coords in zip(*lam.values()))
    assert D_theta(lam, paths) == 1

    current_ring, *coords = ring("b1,b2,b3,b4,b5,c1,c2,c3,c4,c5", QQ)
    b, c = coords[:5], coords[5:]
    currents = [[-u-v for u, v in zip(b, c)]]
    for initial, (vertices, _) in zip((b, c), paths[1:]):
        current = initial
        currents.append(current)
        for vertex in vertices:
            current = [u+p for u, p in zip(current, lam[vertex])]
            currents.append(current)

    # At B1 the external gamma sits either at a c sink or a b source.
    compatible = [(flow, types) for flow, types in flows
                  if types[1] == "c" or all(head != 1 for _, head in flow)]
    assert len(compatible) == 4
    total = sum((grouped_integrand(nv, edges, paths, f, t, currents, external)
                 for f, t in compatible), start=current_ring.zero)
    # Verify the polynomial against actual Kirchhoff currents for two positive edge metrics.
    for shift in (7, 13):
        times = [1 + (3*e + shift) % 17 for e in range(13)]
        raw = kirchhoff(nv, edges, times, lam)
        point = dict(zip(current_ring.symbols, raw[1] + raw[7]))
        literal = sum(grouped_integrand(nv, edges, paths, f, t, raw, external)
                      for f, t in compatible)
        assert total.as_expr().subs(point) == literal != 0
    # The monomial mean routine expects integer coefficients.
    assert all(coefficient.denominator == 1 for coefficient in total.values())
    return total


if __name__ == "__main__":
    for axes, expected in (((2, 0, 0, 0), Fraction(-25457, 278891278848000)),
                           ((2, 3, 4, 1), Fraction(11120243, 48381574026240000))):
        polynomial = transverse_polynomial(axes)
        value = mean(polynomial, SWAPPED) / 158400
        assert value == expected
        print("B2-B5 vector axes (1-based)", tuple(i + 1 for i in axes),
              "monomials", len(polynomial), "ordered kernel / pi", value)
