"""Mixed one-jet cochain, regular coordinates and translation-invariant evaluation."""

from fractions import Fraction
from math import prod
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mixed import GhostBox, EPS, L, jet, word, add


def check_mixed():
    theory = GhostBox(3, {p: Fraction(s) for p, s in EPS.items()})
    for axis in range(3):
        terms = [{}, {}, {}]
        for i in range(3):
            for j, (left, right) in enumerate((("g1", "g2"), ("g2", "g1"), ("b", "c"))):
                terms[j] = add(terms[j], word(L(left, i), L(right, i, jet(axis))))
        assert theory.Q0(terms[2])
        assert not theory.Q0(add(add(terms[0], terms[1], -1), terms[2], 2))
        assert all(not theory.Q1(term) for term in terms)
    w = sp.symbols("w")
    D = lambda expression: w * sp.diff(expression, w)
    x, y = w, 1 / w
    A = 2 * x * x
    C = 2 * x * y
    assert sp.cancel(D(C) - C * D(A) / A) == 2 * (x * D(y) - y * D(x)) == -4


def check_jacobian():
    for n in range(2, 7):
        x = list(range(1, n))
        x.append(-sum(x))
        y = [i**2 for i in range(1, n)]
        y.append(-sum(y))
        size = n - 1
        jac = sp.zeros(2 * size)
        for a in range(1, n):
            for i in range(size):
                c = x[i] ** a - x[-1] ** a
                jac[a - 1, i] = (a + 1) * c
                jac[a - 1 + size, i] = a * (
                    x[i] ** (a - 1) * y[i] - x[-1] ** (a - 1) * y[-1]
                )
                jac[a - 1 + size, i + size] = c
        vand = prod(x[i] - x[j] for i in range(n) for j in range(i + 1, n))
        assert abs(jac.det()) == sp.factorial(n) * vand**2


if __name__ == "__main__":
    check_mixed()
    check_jacobian()
    print(
        "three mixed jet directions and sl2 through sl6 regular coordinate Jacobians passed"
    )
