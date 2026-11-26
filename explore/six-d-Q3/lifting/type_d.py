"""First extra D4 polarization, its Hamiltonian, and native superspace shears."""

from fractions import Fraction
from itertools import product
from math import comb, prod
from random import Random
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from ops import L, Theory
from salg import EPS, add, gen, lderiv, letters_of, mul, scale
from quadratic import exact_rank


def pf_coeff(x, y, k):
    return sum(prod(y[j] if bit else x[j] for j, bit in enumerate(bits))
               for bits in product((0, 1), repeat=len(x)) if sum(bits) == k)


def pfaffian_test():
    rng = Random(8451)
    data = [[tuple(rng.randint(-5, 5) for _ in range(4)) for _ in range(2)]
            for _ in range(34)]
    generators = []
    for d, tag in [(2, "trace"), (4, "trace"), (6, "trace"), (4, "pf")]:
        for k in range(d + 1):
            if d - k > 3 or k > 3:
                continue
            row = [pf_coeff(x, y, k) if tag == "pf" else
                   comb(d, k) * sum(a ** (d - k) * b ** k for a, b in zip(x, y))
                   for x, y in data]
            generators.append((d-k, k, row))
    products = {(0, 0): [[1] * len(data)]}
    for dp, dq, values in generators:
        for p in range(4):
            for q in range(4):
                if p + dp <= 3 and q + dq <= 3:
                    for row in products.get((p, q), ()):
                        products.setdefault((p + dp, q + dq), []).append([a*b for a,b in zip(row, values)])
    rows = [dict(enumerate(v)) for v in products[3, 3]]
    new = [sum(b**3 * prod(x[j] for j in range(4) if i != j) for i, b in enumerate(y))
           for x, y in data]
    rank = exact_rank(rows)
    assert exact_rank([*rows, dict(enumerate(new))]) == rank + 1
    # On x,y multiples of a common Cartan element the formula is also the
    # linear coefficient of Pf(x+t y^3), term by term.
    for k, (x, y) in enumerate(data):
        assert new[k] == pf_coeff(x, tuple(b**3 for b in y), 1)
    return rank


def shear(poly, flavor, axis):
    result = {}
    for letter in letters_of(poly):
        kind, color, *powers = letter
        if kind == 2:
            target, sign = (4, 1) if flavor == 1 else (3, -1)
        elif kind == 2 + flavor:
            target, sign = 1, -1
        else:
            continue
        powers[axis] += 1
        result = add(result, mul(gen((target, color, *powers)), lderiv(poly, letter)), sign)
    return result


def shear_test():
    theory = Theory(3, {key: Fraction(value) for key, value in EPS.items()})
    probes = (mul(L("b", 0), L("b", 1)), mul(L("b", 0), L("g1", 1)))
    count = 0
    for probe, flavor, axis, Q in product(probes, (1, 2), (0, 1, 2), (theory.Q0, theory.Q1)):
        lhs = Q(shear(probe, flavor, axis))
        rhs = shear(Q(probe), flavor, axis)
        assert not add(lhs, rhs)
        count += bool(lhs)
    assert count
    return count


if __name__ == "__main__":
    base_rank = pfaffian_test()
    count = shear_test()
    print(f"D4 bidegree (3,3) polarized-trace/Pfaffian test rank {base_rank} -> {base_rank+1}")
    print(f"24 native shear anticommutators vanish, including {count} nonzero cancellations")
