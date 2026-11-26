"""Exact constant-mode multiplicities, including the project's finite vertex cutoffs."""
import ast
from collections import defaultdict
from math import comb
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def vertex_product(vertices, cutoff):
    coefficients = [1]
    for _ in range(vertices):
        product = [0] * (len(coefficients) + cutoff)
        for power, coefficient in enumerate(coefficients):
            for extra in range(cutoff + 1):
                product[power + extra] += coefficient
        coefficients = product
    return coefficients


def main():
    for m in range(1, 6):
        vertices = 12 - m
        # Truncation at each vertex to 20 is exact through total degree 20.
        coefficients = vertex_product(vertices, 20)
        assert coefficients[:21] == [comb(k + vertices - 1, k) for k in range(21)]
    assert vertex_product(11, 1)[:3] == [1, 11, 55]
    assert vertex_product(11, 2)[:3] == [1, 11, 66]

    rows = defaultdict(dict)
    pattern = r'(\([0-9, ]+\)) m = (\d+) n_a = 0 external (\{[^}]*\}) maps: (\d+)'
    for graph, m, fields, count in re.findall(pattern, (HERE.parent / 'structure.log').read_text()):
        rows[graph, int(m)][ast.literal_eval(fields).get('nu', 0)] = int(count)
    assert len(rows) == 10
    for (_, m), observed in rows.items():
        predicted = vertex_product(12 - m, 1)
        assert observed == {k: observed[0] * n for k, n in enumerate(predicted)}

    path = HERE.parent / 'verify' / 'v_structure.log'
    pattern = r'(\([0-9, ]+\)) Dtheta rule, kmax=(\d+):.*?max K by m: (\{[^}]*\})'
    limits = re.findall(pattern, path.read_text())
    assert len({(graph, cutoff) for graph, cutoff, _ in limits}) == 10
    for _, cutoff, limits_by_m in limits:
        assert ast.literal_eval(limits_by_m) == {m: int(cutoff) * (12 - m) for m in range(1, 6)}
    print('Passed: unrestricted coefficients through K=20 for m=1,...,5; '
          '100 recorded map counts; kmax=1,2 endpoints for five graphs in five sectors.')


if __name__ == '__main__':
    main()
