"""Exact checks of holography.cyclic (Connes complex, cyclic-word CE model) against HKR.

Run:  python3 code/holography/tests/test_cyclic.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from fractions import Fraction
from holography import cyclic as cy, derham as dr
from holography.linalg import vadd

BETAS = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 0, 1, 0, 0), (1, 1, 0, 0, 0), (0, 0, 1, 1, 0), (0, 0, 2, 0, 0),
         (1, 1, 1, 0, 0), (0, 0, 1, 1, 1), (2, 1, 0, 0, 0), (1, 1, 1, 1, 0), (0, 0, 2, 1, 0)]


def b_lambda(x):
    out = {}
    for y, c in cy.hochschild_b(x).items():
        r = cy.canon(y)
        if r:
            vadd(out, {r[0]: Fraction(r[1])}, c)
    return out


def test_b_squared_and_descent():
    for beta in BETAS:
        for n in range(4):
            for x in cy.chains(beta, n):
                acc = {}
                for y, c in cy.hochschild_b(x).items():
                    vadd(acc, cy.hochschild_b(y), c)
                assert not acc
                r = cy.canon(x)
                if r is None:
                    assert not b_lambda(x)
                else:
                    assert b_lambda(x) == {k: r[1] * v for k, v in b_lambda(r[0]).items()}


def test_HC_vs_HKR():
    for beta in BETAS:
        B = dr.B_list(beta)
        for n in range(4):
            Bn = B[n + 1] if n + 1 < len(B) else 0
            assert cy.HC(beta, n)[0] == Bn + (1 if beta == (0,) * 5 and n % 2 == 0 else 0)
            assert cy.HC(beta, n, reduced=True)[0] == Bn


def test_CE_words():
    for beta in BETAS:
        B = dr.B_list(beta)
        for Y in range(1, 5):
            for rel in (False, True):
                for x in cy.words(beta, Y, rel):
                    acc = {}
                    for y, c in cy.Q_word(x, rel).items():
                        vadd(acc, cy.Q_word(y, rel), c)
                    assert not acc
            BY = B[Y] if Y < len(B) else 0
            assert cy.CE_single_trace(beta, Y, relative=True)[0] == BY
            assert cy.CE_single_trace(beta, Y)[0] == BY + (1 if beta == (0,) * 5 and Y % 2 else 0)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
