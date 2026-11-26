"""Check the project component calculus, invariant-jet Jacobian and integrated sample."""

from collections import Counter
from itertools import product
from math import factorial, prod
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verify"))
from vops import GRAPH_SPECS, apply_graph


def check_theta():
    u, v = 2, -3
    count = 0
    for graph, spec in GRAPH_SPECS.items():
        k = len(spec["spokes"])
        bindings = tuple(range(k))
        coefficients = Counter()
        for species in product(("g1", "g2"), repeat=k):
            weight = prod(u if field == "g1" else v for field in species)
            for sign, words in apply_graph(
                graph, "output", tuple(enumerate(species)), bindings
            ):
                coefficients[tuple(words)] += weight * sign
        assert not any(coefficients.values())
        for species in product(("c", "g1", "g2", "b"), repeat=k):
            for _, words in apply_graph(
                graph, "output", tuple(enumerate(species)), bindings
            ):
                assert any(field != "g1" for _, field in words)
            count += 1
    return count


def check_jets():
    for n in range(2, 10):
        h = list(range(1, n))
        h.append(-sum(h))
        jac = sp.Matrix(
            [
                [j * (h[i] ** (j - 1) - h[-1] ** (j - 1)) for i in range(n - 1)]
                for j in range(2, n + 1)
            ]
        )
        assert abs(jac.det()) == factorial(n) * abs(
            prod(h[i] - h[j] for i in range(n) for j in range(i + 1, n))
        )
    # Two distinct derivatives, including the nonlinear mixed term, in sl(3).
    t, s = sp.symbols("t s")
    entries = sp.symbols("u:8")
    h = [
        1 + entries[0] + entries[1] * t + entries[2] * s + entries[3] * s * t,
        2 + entries[4] + entries[5] * t + entries[6] * s + entries[7] * s * t,
    ]
    h.append(-sum(h))
    jets = []
    for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
        for j in (2, 3):
            jets.append(sum(sp.diff(e**j, t, a, s, b).subs({t: 0, s: 0}) for e in h))
    at = {var: 0 for var in entries}
    assert abs(sp.Matrix(jets).jacobian(entries).subs(at).det()) == 120**4
    mode = sp.symbols("w")
    density = sp.expand((mode * sp.diff(2 * (mode + 1 / mode) ** 2, mode)) ** 2)
    assert density.coeff(mode, 0) == -32


if __name__ == "__main__":
    count = check_theta()
    check_jets()
    print(
        count, "component bindings; Casimir jet Jacobian and integrated density passed"
    )
