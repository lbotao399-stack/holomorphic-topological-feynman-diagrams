"""Separated one-loop cubic pentagons with one normal exact-pair propagator block."""
from itertools import permutations
from math import prod

from sympy import Matrix, eye, zeros


# Only vertices 0 and 4 lie on the same tangent diagonal.
POINTS = tuple(map(Matrix, [(0, 0, 1, 1, 1), (1, 0, 1, 1, 2),
                           (0, 1, 1, 2, 1), (-1, 0, 1, 1, 2), (0, 0, 0, 1, 1)]))
UNIT = eye(5)


def row(vertex, covector):
    result = zeros(1, 25)
    result[:, 5*vertex:5*vertex+5] = covector.T
    return result


def angular_rows(target, source, direction):
    # Four linear factors of i_direction(volume) in fixed coordinates.
    k = next(i for i, x in enumerate(direction) if x)
    forms = []
    for j in range(5):
        if j != k:
            covector = direction[k]*UNIT[:, j] - direction[j]*UNIT[:, k]
            forms.append(row(target, covector) - row(source, covector))
    return forms, (-1)**k * direction[k]**3


def coefficient(order, derivative):
    target, source, flux, test = 0, 4, 1, 2
    path = [target, *order, source]
    ordinary = []
    scales, distances = [], []
    for i, j in zip(path, path[1:]):
        delta = POINTS[i] - POINTS[j]
        forms, scale = angular_rows(i, j, delta)
        ordinary.extend(forms)
        scales.append(scale)
        distances.append(delta.dot(delta))

    x, y = POINTS[target][2:, :], POINTS[source][2:, :]
    tangent_delta = [row(target, UNIT[:, k]) - row(source, UNIT[:, k]) for k in (0, 1)]
    special = tangent_delta + [row(target, Matrix([0, 0, *y.cross(x)])),
                               row(source, Matrix([0, 0, *y]))]
    flux_normal = POINTS[flux][2:, :]
    flux_forms = [row(flux, flux_normal[0]*UNIT[:, k+2] - flux_normal[k]*UNIT[:, 2])
                  for k in (1, 2)]
    test_forms = [row(test, UNIT[:, k]) for k in (0, 1, 2)]
    anti = Matrix.vstack(*ordinary, *special, *flux_forms, *test_forms).det(method='domain-ge')
    anti /= prod(scales)*flux_normal[0]

    special_odd, special_scale = angular_rows(target, source, UNIT[:, derivative])
    external = [row(flux, UNIT[:, k]) for k in (0, 1)] + [row(test, UNIT[:, k]) for k in (2, 3, 4)]
    odd = Matrix.vstack(*ordinary, *special_odd, *external).det(method='domain-ge')
    odd /= prod(scales)*special_scale
    return anti, odd, tuple(distances)


def check():
    assert len({tuple(POINTS[j][:2]) for j in range(4)}) == 4
    assert POINTS[4][:2] == POINTS[0][:2]
    assert [POINTS[j][2:, :].dot(POINTS[j][2:, :]) for j in (0, 4)] == [3, 2]
    for derivative in (0, 1):
        terms = [coefficient(order, derivative) for order in permutations((1, 2, 3))]
        assert all(distances == (2, 4, 4, 3) for _, _, distances in terms)
        assert all(anti == -2 and odd == 2*derivative for anti, odd, _ in terms)
        total = sum(anti*odd/prod(distances)**6 for anti, odd, distances in terms)
        assert total*32614907904 == -derivative
        print(f"Tangent delta derivative {derivative+1}, six paths, reduced coefficient {total}.")


if __name__ == '__main__':
    check()
