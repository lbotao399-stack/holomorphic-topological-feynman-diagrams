"""Convergent projective K4 monomial integrals in every even cut dimension."""
from fractions import Fraction as F
from functools import cache
from itertools import permutations
from pathlib import Path
import sys

from flint import fmpq_poly

sys.path.insert(0,str(Path(__file__).resolve().parents[4]/'deepen/o3/higher_orders/four_loop_reduction'))
from reconstruct_prism.evaluate import Rational
from kite_native import native_kite

EDGES = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
PERMS = [tuple(EDGES.index(tuple(sorted((p[i],p[j])))) for i,j in EDGES)
         for p in permutations(range(4))]


def orbit(powers):
    return min(tuple(powers[i] for i in p) for p in PERMS)


def gamma_ratio(n,slope):
    """Rational factor in Gamma(n+slope*ep) / Gamma(1+slope*ep)."""
    z = fmpq_poly([1])
    for j in range(min(1,n),max(1,n)):
        z *= fmpq_poly([j,slope])
    return Rational(z) if n>=1 else Rational(1,z)


@cache
def prefactors(q):
    x = gamma_ratio(2-q,1)
    y = gamma_ratio(q-1,-1)
    z = gamma_ratio(2*q-2,-2)
    first = x*x*y*y*y*y/(z*z)
    second = y*y*y*gamma_ratio(3-2*q,2)/gamma_ratio(3*q-3,-3)
    return first,second


@cache
def cut(powers):
    """Two-loop propagator cut at 01, with its weights and Q^2=1."""
    q,r = divmod(sum(powers),3)
    assert not r and min(powers)>=1 and max(powers)<q
    for triangle in ((0,1,3),(0,2,4),(1,2,5),(3,4,5)):
        assert sum(powers[j] for j in triangle)>q
    first,second,common = native_kite(*powers[1:])
    dimension = fmpq_poly([q,-1])
    denominator = fmpq_poly(common)(dimension)
    x,y = prefactors(q)
    a = Rational(fmpq_poly(first)(dimension),denominator)*x
    b = Rational(fmpq_poly(second)(dimension),denominator)*y
    assert b.valuation()>=-3
    summed = (a+b).series(0)
    assert all(i>=0 for i in summed)
    return F(summed.get(0,0)), -6*F(b.series(-3).get(-3,0))


@cache
def checked_cut(powers):
    first = cut(powers)
    # Reverse the cut and the central propagator by a K4 vertex permutation.
    other = tuple(powers[j] for j in (5,4,2,3,1,0))
    assert other in {tuple(powers[j] for j in p) for p in PERMS}
    assert first==cut(other)
    return first
