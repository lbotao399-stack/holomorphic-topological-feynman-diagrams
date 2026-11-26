"""Nonrecursive two-loop counterterms in all four trace/rank sectors."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, prod

from sympy import bernoulli

from probe_laurent import RibbonProbe, SINGLE


COLLINEAR = ((1,0),(0,1),(0,1),(-1,-2))


@lru_cache(maxsize=None)
def bare_double(n):
    if n < 1 or not n%2:
        return F(0)
    m = n+3
    # Beta integral on the (0)(123) probe; Gamma(1/2) cancels.
    def rising_half(length):
        return prod((F(2*i+1,2) for i in range(length)), start=F(1))
    return F(3,8)/rising_half(m-1)*sum(
        rising_half(i-2)/factorial(i) for i in range(2,m-1,2))


@lru_cache(maxsize=None)
def double_quotient(j):
    # [x^j] Theta(x)/sinh(x), via the expansion of x/sinh(x).
    value = F(0)
    for i in range(1,j+2,2):
        m = j+1-i
        inverse = 2*(1-F(2)**(m-1))*F(bernoulli(m))/factorial(m)
        value += bare_double(i)*inverse
    return value


def double(n):
    if not n%2:
        return F(0)
    homogeneous = F(1,2**(n+2)*factorial(n+2))
    return homogeneous + sum(double_quotient(j) /
                             ((2**(j+2)-1)*2**(n-j)*factorial(n-j))
                             for j in range(0,n,2))


def disk(n):
    k = n+2
    return -F(k*(k+4),10*2**k*factorial(k))


def triple(n):
    k = n+2
    return -F(6,(k+2)*2**k*factorial(k)) if k%2 else F(0)


def genus_one(max_order):
    bare = [F(0)]*(max_order+1)
    for shape in ((1,3,3),(2,2,3)):
        probe = RibbonProbe(shape, COLLINEAR, SINGLE)
        probe.integrate(range(max_order+1))
        for n in range(max_order+1):
            bare[n] += probe.moments[n].get(0,0)
    # First multiply the graph by exp(-e). The dilation 2e produces a
    # diagonal divisor on each Taylor coefficient. Then multiply by exp(e/2).
    reduced = [sum((-1)**(j-r)*bare[r]/factorial(j-r) for r in range(j+1)) /
               (2*(2**(j+1)-1)) for j in range(max_order+1)]
    return {n:sum(reduced[j]/(2**(n-j)*factorial(n-j)) for j in range(n+1))
            for n in range(max_order+1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("max_order", type=int, default=2, nargs="?")
    args = parser.parse_args()
    if args.max_order < 0:
        parser.error("max_order must be nonnegative")
    handle = genus_one(args.max_order)
    print("n   N^2 single       N^0 single       N double       triple")
    for n in range(-1,args.max_order+1):
        print(n, disk(n), handle.get(n,F(0)), double(n), triple(n))


if __name__ == "__main__":
    main()
