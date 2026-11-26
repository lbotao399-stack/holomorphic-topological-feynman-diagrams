"""Independent exact four-flow integration for two transverse theta_166 fixtures.

Only the project's literal Wick words and Kirchhoff solver are used. Each ordered
path is integrated with five independent beta ratios. The conductance simplex
is integrated using two separate beta ratios with a positive binomial expansion.
"""

import sys
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import comb, prod
from pathlib import Path

from sympy import QQ
from sympy.polys.rings import ring

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verify"))
from vcommon import D_theta, kirchhoff, theta_graph
from v_integrand import flows_m1, grouped_integrand


def rising(x, n):
    return prod((x+i for i in range(n)), start=Fraction(1))


@cache
def beta(a, b, p, q):
    """E[X**p (1-X)**q] when X has normalized Beta(a,b) law."""
    return rising(a,p)*rising(b,q)/rising(a+b,p+q)


@cache
def tail_moment(powers):
    # Five descending uniform order statistics t_i have independent ratios
    # t_i/t_(i-1) ~ Beta(6-i,1), i=1,...,5 and t_0=1.
    return prod((beta(5-j,1,sum(powers[j:]),0) for j in range(5)), start=Fraction(1))


@cache
def weight_moment(bs, bt, cs, ct):
    # b_i=-(1-wB) S_i - wB T_i; c_i=wC S_i + (1-wC) T_i.
    # p=wB ~ Beta(1/2,6), q=wC/(1-wB) ~ Beta(1/2,11/2), independently.
    # Expand (1-wC)^ct = [p+(1-p)(1-q)]^ct with positive coefficients.
    half = Fraction(1,2)
    return (-1)**(bs+bt)*sum(
        comb(ct,k)*beta(half,6,bt+k,bs+cs+ct-k)*beta(half,11*half,cs,ct-k)
        for k in range(ct+1))


def integrate(poly, order):
    total = Fraction(0)
    for exponent, coefficient in poly.items():
        factors = [i for i, power in enumerate(exponent) for _ in range(power)]
        moment = Fraction(0)
        for choices in product((0,1),repeat=len(factors)):
            powers = [[0]*5, [0]*5]
            counts = [0]*4
            for i, path in zip(factors, choices):
                powers[path][i%5] += 1
                counts[2*(i//5)+path] += 1
            bpath = tuple(powers[0])
            cpath = tuple(powers[1][j] for j in order)
            moment += tail_moment(bpath)*tail_moment(cpath)*weight_moment(*counts)
        total += Fraction(int(coefficient.numerator),int(coefficient.denominator))*moment
    return total


def parity(seq):
    return (-1)**sum(x>y for i,x in enumerate(seq) for y in seq[i+1:])


def stabilizer(legs):
    """Signed coordinate permutations preserving momentum and polarization lines.

    Each leg is (signed momentum axis, unsigned polarization axis), one-based.
    The gamma has (+1,2). Its magnitude and all momentum signs are exact;
    polarization reversal is recorded by multilinearity.
    """
    preserved = []
    for permutation in permutations(range(1,6)):
        if permutation[:2] != (1,2):
            continue
        for tail in product((1,-1),repeat=4):
            signs = (1,)+tail
            image = [(signs[abs(k)-1]*permutation[abs(k)-1]*(1 if k>0 else -1),
                      permutation[u-1]) for k,u in legs]
            if set(image) != set(legs):
                continue
            permuted_odd = [legs.index(x) for x in image]
            polarization = signs[1]*prod(signs[u-1] for _,u in legs)
            det = parity(permutation)*prod(signs)
            character = det**3*parity(permuted_odd)*polarization
            preserved.append((permutation,signs,character))
    return preserved


def check_fixture(name, polarizations, expected):
    n,edges,paths = theta_graph(1,6,6)
    order = (0,1,3,2,4)
    basis = [[int(i==j) for j in range(5)] for i in range(5)]
    # Positive and negative momentum sequences in the circuit path order.
    axes = {v+1:v+1 for v in range(5)}
    axes.update({v+6:-j-1 for v,j in enumerate(order)})
    lam = {v:[(1 if i>0 else -1)*x for x in basis[abs(i)-1]] for v,i in axes.items()}
    lam.update({0:[0]*5,11:[0]*5})
    ext = {v:dict(lam=lam[v],u=basis[polarizations.get(v,1)-1],vg=basis[1]) for v in range(n)}
    assert all(sum(t) == 0 for t in zip(*lam.values())) and D_theta(lam, paths) == 1
    legs = [(axes[v], polarizations.get(v,1)) for v in range(2,11)]
    assert all(sum(x*y for x,y in zip(lam[v],ext[v]['u'])) == 0 for v in range(2,11))
    group = stabilizer(legs)
    if name == 'II':
        assert group == [((1,2,3,4,5),(1,1,1,1,1),1)]
    else:
        assert len(group) == 16 and sum(char < 0 for _,_,char in group) == 8

    R,*variables = ring('b1,b2,b3,b4,b5,c1,c2,c3,c4,c5',QQ)
    b,c = variables[:5],variables[5:]
    currents = [[-p-q for p,q in zip(b,c)]]
    for start,path in zip((b,c),paths[1:]):
        cur = start
        currents.append(cur)
        for v in path[0]:
            cur = [x+y for x,y in zip(cur,lam[v])]
            currents.append(cur)
    physical = []
    for r in (13,17):
        times = [Fraction(1+(e*e+3*e+5)%r,1+(5*e+1)%7) for e in range(13)]
        physical.append(kirchhoff(n,edges,times,lam))

    total = R.zero
    count = 0
    for flow,types in flows_m1(n,edges):
        cv = next(v for v in types if types[v]=='c')
        if not (cv==1 or (cv in (0,11) and all(h!=1 for _,h in flow))):
            continue
        poly = grouped_integrand(n,edges,paths,flow,types,currents,ext)
        for y in physical:
            actual = grouped_integrand(n,edges,paths,flow,types,y,ext)
            value = sum(Fraction(int(c.numerator),int(c.denominator))
                        * prod(x**p for x,p in zip(y[1]+y[7],ex)) for ex,c in poly.items())
            assert actual == value
        count += 1
        total += poly
    assert count == 4
    mean = integrate(total,order)
    coefficient = mean/158400
    assert coefficient == expected
    print(name,"current monomials",len(total),"mean",mean,"kernel/pi",coefficient,
          "signed stabilizer order",len(group),"negative characters",
          sum(char<0 for _,_,char in group))


if __name__ == '__main__':
    check_fixture('II',{2:3,3:4,4:5,5:2,6:2},Fraction(11120243,48381574026240000))
    check_fixture('I',{2:3,6:2},Fraction(-25457,278891278848000))
