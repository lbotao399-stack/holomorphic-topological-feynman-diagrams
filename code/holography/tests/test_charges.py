"""Exact checks of holography.charges: LLC letter charges and the KMMR index.

Run:  python3 code/holography/tests/test_charges.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from fractions import Fraction
from holography import charges as ch, derham as dr
from holography.linalg import vadd


def test_letter_charges():
    # LLC Table appE2-0020: gamma^1, beta_1, b, d_{+dot+} c  (Y = 1)
    assert ch.llc_charges((1, 0, 0), (0, 0), 1) == (1, (1, 0, 0), 0, 0)
    E, R, JL, JR = ch.llc_charges((0, 1, 1), (0, 0), 1)
    assert (E, R, JL, JR) == (Fraction(3, 2), (Fraction(-1, 2), Fraction(1, 2), Fraction(1, 2)), Fraction(1, 2), 0)
    assert ch.llc_charges((1, 1, 1), (0, 0), 1) == (2, (0, 0, 0), 1, 0)
    E, R, JL, JR = ch.llc_charges((0, 0, 0), (1, 0), 1)
    assert (E, JL, JR) == (Fraction(3, 2), 0, Fraction(1, 2)) and R == (Fraction(1, 2),) * 3


def series(tmax, pmin):
    out = {}
    for beta in ch.sector_list(tmax):
        for p, dim in enumerate(dr.B_list(beta)):
            if p >= pmin and dim:
                s, m = ch.kmmr_monomial(beta[2:], beta[:2], p)
                vadd(out, {m: Fraction(s * dim)})
    return out


def test_kmmr_low_order():
    tmax = 12
    f = ch.kmmr_f(tmax)
    Isp = ch.kmmr_Isp(tmax)
    B1 = vadd(series(tmax, 1), series(tmax, 2), -1)
    assert B1 == f
    assert series(tmax, 1) == Isp
    assert series(tmax, 2) == vadd(dict(Isp), f, -1)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
