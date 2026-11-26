"""Four field flows for one fixed lam, gamma and mu assignment on the 166 graph."""
from fractions import Fraction

from check_theta166_flow import B, C, CURRENT_RING, IDENTITY, SWAPPED, mean, reference_setup


def polynomials():
    x, z = B[0], C[0]
    area = (B[2] - B[4]) * (C[3] - C[4]) - (B[3] - B[4]) * (C[2] - C[4])
    common = (1 + x)**3 * (1 - z)**4 * area
    factors = {
        (0, 1): -2*x*x - 4*x*z + x - 4*z*z + 3*z,
        (0, 11): -(x + 2*z) * (2*x + 1),
        (11, 1): 2*x*x + 4*x*z + x + 4*z*z - z,
        (11, 0): (2*x + 1) * (x + 2*z - 1),
    }
    assert sum(factors.values()) == 2*z - 1
    return {key: common * factor for key, factor in factors.items()}


def main():
    signed = polynomials()
    expected = {(0, 1): -65558839, (0, 11): -11982421,
                (11, 1): 85509601, (11, 0): -27353656}
    for order in (SWAPPED, IDENTITY):
        nv, edges, paths, lam, external, flows = reference_setup(order)
        from vcommon import D_theta, kirchhoff
        from v_integrand import grouped_integrand
        determinant = D_theta(lam, paths)
        seen = set()
        for flow, types in flows:
            incoming_b1 = sum(head == 1 for _, head in flow)
            if (3 if types[1] == 'c' else 1) == incoming_b1:
                continue
            c_vertex, = [v for v in range(nv) if types[v] == 'c']
            key = flow[0][0], c_vertex
            assert key not in seen
            seen.add(key)
            polynomial = signed[key]
            for shift in (4, 6, 11):
                times = [1 + (3*e + shift) % 17 for e in range(13)]
                y = kirchhoff(nv, edges, times, lam)
                currents = dict(zip(polynomial.ring.symbols, y[1] + y[7]))
                reference = grouped_integrand(nv, edges, paths, flow, types, y, external)
                assert polynomial.as_expr().subs(currents) == reference != 0
            averaged = mean(polynomial, order)
            assert averaged == (Fraction(expected[key], 44464437108000) if order == SWAPPED else 0)
            print(order, key, 'mean', averaged)
        assert seen == set(signed)
        total = determinant * mean(sum(signed.values(), CURRENT_RING.zero), order) / 158400
        assert total == (Fraction(-3877063, 1408633367581440000) if order == SWAPPED else 0)
        print(order, 'four-flow kernel/pi', total)


if __name__ == '__main__':
    main()
