"""Four-corner ladder density and projected current moments in the S4 edge order."""
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
from htk import brute
from htk.graph import ladder


def period(n):
    return Fraction((-1)**(n-1) * 2**(n+1),
                    factorial(n+2)*factorial(n)**2*factorial(n+1)**2)


def beta(i, j):
    return Fraction(factorial(i-1)*factorial(j-1), factorial(i+j-1))


def projected_moment(powers):
    """Insert powers of downward rung currents, projected with det(A,B,-)/det(A,B,C)."""
    moment = period(len(powers)-1) * (-1)**sum(powers)
    total = powers[0]
    for j, power in enumerate(powers[1:], 1):
        moment *= beta(3+power,j+total) / beta(3,j)
        total += power
    return moment


def times(ratios):
    # Fix the final effective resistance to 1 and undo the parallel reductions.
    rungs, upper, lower = {}, {}, {}
    jacobian = effective = 1
    for j in range(len(ratios), 0, -1):
        rho, s1, s2 = ratios[j-1]
        jacobian *= effective**3 / (rho**2*(1-rho)**4)
        series = effective/(1-rho)
        rungs[j] = effective/rho
        upper[j], lower[j] = series*(1-s1), series*s2
        effective = series*(s1-s2)
    rungs[0] = effective
    original = []
    n = len(ratios)
    for i in range(n-1):
        original.extend((upper[i+1],rungs[i],lower[i+1]))
    return original + [rungs[n-1],lower[n],rungs[n],upper[n]], jacobian


def verify():
    rng = Random(94)
    for n in range(1,9):
        graph, cotree, a, b = ladder(n)
        # Unit determinant at the four corners. The source at b(n) balances.
        sources = {a(0): (1,0,0), b(0): (0,1,0), a(n): (-1,-1,-1)}
        ratios = []
        for _ in range(n):
            s1, s2 = sorted(rng.sample(range(1,16),2), reverse=True)
            ratios.append((Fraction(rng.randrange(1,16),16),Fraction(s1,16),Fraction(s2,16)))
        t, jacobian = times(ratios)
        raw = brute.jacobian_density(graph,t,sources,cotree)
        expected = (-1)**(n-1)*prod(rho**2*(1-rho)**(j-1)*(s1-s2)**(n-j)
                                    for j,(rho,s1,s2) in enumerate(ratios,1))
        assert raw*jacobian == expected
        currents = brute.currents(graph,t,sources)
        for i in range(n+1):
            edge = next(k for k,(u,v) in enumerate(graph.edges) if {u,v} == {a(i),b(i)})
            sign = 1 if graph.edges[edge][0] == a(i) else -1
            share = ratios[i-1][0] if i else 1
            assert currents[edge][2]*sign == -share*prod(1-r[0] for r in ratios[i:])
        integrated = (-1)**(n-1)*prod(beta(3,j)/((n-j+1)*(n-j+2)) for j in range(1,n+1))
        assert integrated == period(n)
        first = projected_moment([1]+[0]*n)
        assert first == -period(n)*Fraction(6*factorial(n),factorial(n+3))
        print(n, 'period', period(n), 'left rung first moment', first)


if __name__ == '__main__':
    verify()
