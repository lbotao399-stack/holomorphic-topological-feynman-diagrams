"""Exact five-probe certificate for the three-loop Ward coefficients."""
import random
from fractions import Fraction

import numpy as np
from sympy import Matrix

from vops import Op, pair, so_n


def verify():
    dim, structure = so_n(5)
    rng = random.Random(3)
    operator = [(('O', k), 'b') for k in range(4)]
    colours = {tag: [rng.randint(-2, 2) for _ in range(dim)] for tag, _ in operator}
    cache = {}
    rows = []
    for probe in range(5):
        types = ['c'] * 4 + ['g1', 'g2'] if probe < 3 else ['c'] * 5 + ['b']
        vectors = np.array([[rng.randint(-2, 2) for _ in range(dim)] for _ in types])
        momenta = [[rng.randint(-4, 4) for _ in range(3)] for _ in types]
        def anticommutator(first, second):
            a, b = Op(first, 'u'), Op(second, 'd')
            args = operator, colours, {}, (types, vectors, momenta), dim, structure, cache
            return pair([a, b], *args) + pair([b, a], *args)

        row = [anticommutator('q0', g) for g in ('G1', 'G2', 'G3', 'G4', 'G5')]
        rhs = -(anticommutator('box', 'A') + anticommutator('box', 'B')) / 4
        rows.append(row + [rhs])
        print('probe', probe + 1, *rows[-1], flush=True)

    matrix = Matrix(rows)
    weights = Matrix([Fraction(-1, 2), Fraction(1, 2), -1, Fraction(1, 8), Fraction(-1, 2)])
    assert matrix[:, :5] * weights == matrix[:, 5]
    residue = matrix[:, :5].applyfunc(lambda q: int(q.p) * pow(int(q.q), -1, 101) % 101)
    determinant = int(residue.det()) % 101
    assert determinant
    print('rank minor modulo 101:', determinant, flush=True)
    print('weights:', tuple(weights))


if __name__ == '__main__':
    verify()
