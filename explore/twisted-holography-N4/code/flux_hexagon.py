"""Polarized Hessian determinant with one normal Laurent flux on the fivefold."""
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import factorial, prod


def signed_orders(n):
    for order in permutations(range(n)):
        inversions = sum(x > y for i, x in enumerate(order) for y in order[i + 1:])
        yield order, (-1)**inversions


def normal_jet(p, q):
    tangents = ((1, 0, 0, 0, 0), (0, 1, 0, 0, 0))
    waves = (*tangents, (0, 0, *p), (0, 0, *q))
    result = {ij: 0 for ij in combinations_with_replacement(range(3), 2)}
    # Expand the determinant coefficient with four labeled quadratic fields.
    # The last row assignment puts the remaining Hessian on the flux.
    for columns, sign in signed_orders(5):
        for rows in permutations(range(5)):
            i = rows[-1]
            j = columns[i]
            if min(i, j) < 2:
                continue
            coefficient = prod(wave[r] * wave[columns[r]] for wave, r in zip(waves, rows))
            result[tuple(sorted((i - 2, j - 2)))] += sign * coefficient
    return result


def check():
    probes = (((0,1,0),(0,0,1)), ((0,0,1),(1,0,0)), ((1,0,0),(0,1,0)),
              ((0,0,1),(1,-1,0)), ((0,1,0),(-1,0,1)), ((1,0,0),(0,1,-1)))
    coordinates = []
    for p, q in probes:
        actual = normal_jet(p, q)
        v = tuple(p[(i+1)%3]*q[(i+2)%3] - p[(i+2)%3]*q[(i+1)%3] for i in range(3))
        assert actual == {(i,j): (1 if i == j else 2)*v[i]*v[j] for i,j in actual}
        forms = []
        for (i,j), weight in actual.items():
            exponents = [-1, -1, -1]
            laurent = weight
            for k in (i, j):
                laurent *= exponents[k]
                exponents[k] -= 1
            forms.append(Fraction(laurent, prod(factorial(-1-e) for e in exponents)))
        assert forms == list(actual.values())
        coordinates.append(forms)
    determinant = sum(sign * prod(row[j] for row,j in zip(coordinates,order))
                      for order,sign in signed_orders(6))
    assert abs(determinant) == 8
    assert Fraction(normal_jet(*probes[0])[(0,0)] * 2, 24) == Fraction(1,12)
    # After selecting the normal Hessian slot the off-diagonal z-w block
    # has four labeled entries. Extract their product from the full determinant.
    labels = ((-1,-1,0,1), (-1,-1,2,3), (0,2,-1,-1), (1,3,-1,-1))
    mixed = sum(sign for order,sign in signed_orders(4)
                if sorted(row[j] for row,j in zip(labels,order)) == list(range(4)))
    assert Fraction(mixed * 2, 24) == -Fraction(1,6)
    print("Six independent responses; the two normal residues are N/12 and -N/6.")


if __name__ == "__main__":
    check()
