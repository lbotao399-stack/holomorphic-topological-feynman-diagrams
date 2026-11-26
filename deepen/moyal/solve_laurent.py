"""Determine the two-loop Moyal counterterms by equality of formal trace cochains."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import factorial

import sympy as sp

from check_laurent import ce, cubic, laurent_terms, mark
from theta_moyal import Theta, calculate, ribbon, serialize, to_ce


@lru_cache(maxsize=None)
def classical(degree):
    tensor = "d" if degree % 2 else "f"
    return ce.scale(cubic(degree, tensor), F(1, 2**degree * factorial(degree)))


def verify_disk_phases():
    counts = []
    for shape in ((1, 3, 3), (2, 2, 3)):
        theta = Theta(shape)
        disk_count = 0
        for signs in product((1, -1), repeat=6):
            vertices = {v: [(sign, [[0, 1, 2] if sign == 1 else [0, 2, 1]])]
                        for v, sign in enumerate(signs)}
            colour = ribbon.colour_weight(6, theta.edges, theta.ext, vertices)
            [(words, ranks)] = colour.items()
            [(rank, _)] = ranks.items()
            if rank != 2 or len(words) != 1:
                continue
            disk_count += 1
            phase = sum(sign * vertex for sign, vertex in zip(signs, theta.vertices))
            boundary = 0
            for i, j in combinations(words[0], 2):
                p, q = theta.momenta[i], theta.momenta[j]
                for (k, l), letter in zip(((0, 1), (0, 2), (1, 2)), theta.brackets):
                    boundary += (p[k]*q[l] - p[l]*q[k]) * letter
            assert phase == boundary
        counts.append(disk_count)
    assert counts == [2, 0]


def solve_at_order(n, terms, bare):
    # D and S1 have independent trace types; d S0 vanishes.  At the next
    # even order {j_M, S0} fixes its previously undetermined coefficient.
    residual = bare
    for j in range(-1, n):
        residual = ce.lin((1, residual), (1, mark(classical(n-j), terms[j])))
    columns = {name: ce.d_ce(cubic(n+2, name), 3)
               for name in (("d", "s1") if n % 2 else ("f",))}
    if n % 2 == 0:
        columns["previous s0"] = ce.scale(mark(classical(1), cubic(n+1, "s0")), -1)
    coordinates = set(residual).union(*columns.values())
    assert all(set(value) == {0} for column in columns.values() for value in column.values())
    powers = sorted({rank for value in residual.values() for rank in value})
    coefficients = {name: {} for name in columns}
    for rank in powers:
        rows = {tuple(column.get(key, {}).get(0, 0) for column in columns.values())
                + (residual.get(key, {}).get(rank, 0),) for key in coordinates}
        matrix = sp.Matrix(sorted(rows))
        solution = sp.linsolve((matrix[:, :-1], matrix[:, -1]))
        assert len(solution) == 1, (n, rank, solution)
        for name, value in zip(columns, next(iter(solution))):
            assert value.is_Rational, (n, name, value)
            if value:
                coefficients[name][rank] = F(value)

    single = "d" if n % 2 else "f"
    degree = n + 2
    assert coefficients[single][2] == -F(degree*(degree+4), 10 * 2**degree * factorial(degree))
    terms[n] = {}
    for name, polynomial in coefficients.items():
        if name == "previous s0":
            terms[n-1] = ce.lin((1, terms[n-1]), (polynomial, cubic(n+1, "s0")))
        else:
            terms[n] = ce.lin((1, terms[n]), (polynomial, cubic(n+2, name)))
    check = ce.lin((-1, ce.d_ce(terms[n], 3)), (1, bare))
    for j in range(-1, n):
        check = ce.lin((1, check), (1, mark(classical(n-j), terms[j])))
    assert not check
    return coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("max_order", type=int, default=2, nargs="?")
    args = parser.parse_args()
    if args.max_order < 2:
        parser.error("max_order must be at least 2")
    verify_disk_phases()
    print("Disk phases verified on both theta shapes.", flush=True)
    terms = dict(zip((-1, 0, 1), laurent_terms()))
    for n in range(2, args.max_order+1):
        bare = to_ce(serialize(calculate(n)))
        coefficients = solve_at_order(n, terms, bare)
        for name, polynomial in coefficients.items():
            print(f"e^{n}, {name}: " + ", ".join(f"N^{rank} * {c}" for rank,c in polynomial.items()),
                  flush=True)


if __name__ == "__main__":
    main()
