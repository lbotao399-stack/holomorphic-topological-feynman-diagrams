"""Brackets of regular scalar sources with primary Laurent states at D3 flux."""
from collections import defaultdict
from fractions import Fraction
from math import factorial

from single_particle import contract, div_A2


def divergence(terms):
    out = defaultdict(Fraction)
    for term, coeff in terms.items():
        for key, value in div_A2(term).items():
            out[key] += coeff * value
    return {term: coeff for term, coeff in out.items() if coeff}


def source_bracket(terms, normal, degree):
    """Contract w_normal^(degree-1) dw_normal/(degree-1)! into each state."""
    out = defaultdict(Fraction)
    for (tangent, powers, wedge), coeff in terms.items():
        entry = contract(wedge, normal + 2)
        if entry is None or powers[normal] < degree - 1:
            continue
        sign, remaining = entry
        lowered = list(powers)
        lowered[normal] -= degree - 1
        out[tangent, tuple(lowered), remaining] += sign * coeff / factorial(degree - 1)
    return {term: coeff for term, coeff in out.items() if coeff}


def check():
    for tangent in ((1,0), (0,3), (3,2)):
        for pole in (1,3,5,7):
            powers = (1,0,pole-1)
            target = div_A2((tangent, powers, (0,1)))
            for normal in range(3):
                for degree in range(1,6):
                    raised = list(powers)
                    raised[normal] += degree - 1
                    primitive = (tangent, tuple(raised), (0,1,normal+2))
                    state = divergence({primitive: -Fraction(factorial(degree-1))})
                    assert state and not divergence(state)
                    assert source_bracket(state,normal,degree) == target
    # The third-order pole in the two-loop five-source fixture.
    primitive = ((1,0),(0,3,2),(0,1,3))
    state = divergence({primitive: -Fraction(1,216)})
    output = source_bracket(state,1,4)
    assert output == {((0,0),(0,0,2),(1,)): Fraction(1,1296)}
    assert not divergence(state)
    print('The fourth-order source change contributes +1/1296 to the selected five-source term.')


if __name__ == '__main__':
    check()
