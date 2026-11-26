"""Test the three S5 summands of the genus-four two-soft polarized Betti determinant."""

from fractions import Fraction
from itertools import product
from pathlib import Path
import sys

from flint import nmod_mat

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "iib-bcov"))
from worldsheet_fourfold_positivity import HOLOMORPHIC_SHIFT, POINTS, PRIME


SIMPLEX = [tuple(int(i == j) for j in range(4)) for i in range(4)] + [(-1,) * 4]
SLOTS = (4, 1, 0, 8, 2, 3, 7)


def evaluation_rows():
    p = PRIME
    roots_x = [a for a in range(p) if pow(a, 5, p) == 1]
    roots_y = [b for b in range(p) if pow(b, 3, p) == p - 1]
    for x, y in POINTS:
        assert (y**3 - x**5 + 1) % p == 0
        dy = 5 * x**4 * pow(3 * y**2, -1, p)
        holomorphic = (
            pow(y, -1, p),
            pow(y * y, -1, p),
            x * pow(y * y, -1, p),
            x * x * pow(y * y, -1, p),
        )
        poles = [pow(3 * (x - a), -1, p) for a in roots_x]
        poles += [dy * pow(5 * (y - b), -1, p) for b in roots_y] + [0]
        logarithmic = [poles[i] for i in SLOTS]
        fixed = [
            logarithmic[i]
            - logarithmic[4]
            + sum(a * h for a, h in zip(shift, holomorphic))
            for i, shift in enumerate(HOLOMORPHIC_SHIFT)
        ]
        soft = logarithmic[5] + logarithmic[6] - 2 * logarithmic[4]
        yield holomorphic, fixed, soft


def quadratic_coefficient(direction, rows):
    """Each target-space direction changes four columns, giving degree at most four."""

    def determinant(t):
        matrix = []
        for forms, fixed, soft in rows:
            meromorphic = [a + t * b * soft for a, b in zip(fixed, direction)]
            matrix.append([a * b % PRIME for a, b in product(forms, meromorphic)])
        return int(nmod_mat(matrix, PRIME).det())

    assert determinant(0) == 0
    # Two soft-mark residues imply divisibility by t^2.
    assert (
        8 * (determinant(1) - determinant(-1)) - (determinant(2) - determinant(-2))
    ) % PRIME == 0
    return (
        (16 * (determinant(1) + determinant(-1)) - (determinant(2) + determinant(-2)))
        * pow(24, -1, PRIME)
    ) % PRIME


def projection_norms(v):
    centered = [a - Fraction(sum(v), 5) for a in v]
    second = sum(a * a for a in centered)
    fourth = sum(a**4 for a in centered)
    trivial = second**2 / 4
    standard = Fraction(5, 3) * (fourth - second**2 / 5)
    return trivial, standard, second**2 - trivial - standard


def check():
    assert projection_norms([1, 0, 0, 0, 0]) == (Fraction(4, 25), Fraction(12, 25), 0)
    assert projection_norms([1, 1, 0, 0, 0]) == (Fraction(9, 25), Fraction(2, 25), 1)
    rows = list(evaluation_rows())
    singles = [quadratic_coefficient(v, rows) for v in SIMPLEX]
    assert singles == [21, 91, 144, 205, 120]
    singlet = sum(singles) % PRIME
    standard = (singles[0] - singles[1]) % PRIME
    pairs = [(0, 1), (2, 3), (0, 2), (1, 3)]
    values = [
        quadratic_coefficient(
            tuple(a + b for a, b in zip(SIMPLEX[i], SIMPLEX[j])), rows
        )
        for i, j in pairs
    ]
    five = (values[0] + values[1] - values[2] - values[3]) % PRIME
    print(
        "Singles",
        singles,
        "pairs",
        values,
        "singlet, standard, five:",
        singlet,
        standard,
        five,
    )
    assert singlet and standard and five


if __name__ == "__main__":
    check()
