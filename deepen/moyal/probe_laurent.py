"""All formal-e two-loop counterterms from two four-input trace probes."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, prod

import sympy as sp
from sympy.polys.rings import ring

from check_laurent import ce, laurent_terms
from theta_moyal import Theta, ribbon


SINGLE = ((0,1,2,3),)
DOUBLE = ((0,), (1,2,3))
TRIPLE = ((0,1), (2,), (3,))
CYCLE = ((1,0), (0,1), (1,1), (-2,-2))
PAIRS = ((1,0), (0,1), (-1,0), (0,-1))
COLOURS = {name: [(ce.canon_tau(cycles), sign) for sign, cycles in terms]
           for name, terms in ribbon.TENS.items()}


def det(p, q):
    return p[0] * q[1] - p[1] * q[0]


class DirichletMoments:
    """Integer moment tables sharing one denominator at one Taylor degree."""
    def __init__(self, shape, degree):
        def rising(first, count, step=1):
            return prod(first+step*i for i in range(count))
        self.shares = {}
        self.denominator = rising(7,degree,2)
        for i in range(degree+1):
            for j in range(degree+1-i):
                self.shares[i,j] = (rising(7-2*shape[0],i,2)*rising(7-2*shape[1],j,2)
                                    *rising(7+2*(i+j),degree-i-j,2))
        self.paths = []
        start = 2
        for length in shape:
            if length == 1:
                continue
            # The phase has at most two path ratios and one share per monomial.
            bound = 2*degree
            divisor = rising(length,bound)
            self.denominator *= divisor
            table = {powers:prod(factorial(v) for v in powers)*rising(length+sum(powers),bound-sum(powers))
                     for powers in product(range(bound+1),repeat=length-1) if sum(powers) <= bound}
            self.paths.append((slice(start,start+length-1),table))
            start += length-1

    def integrate(self, polynomial):
        numerator = 0
        for exponents, integer in polynomial.items():
            weight = self.shares[exponents[:2]]
            for block,table in self.paths:
                weight *= table[exponents[block]]
            numerator += int(integer)*weight
        return F(numerator,self.denominator)


class RibbonProbe:
    def __init__(self, shape, momenta, word):
        self.theta = theta = Theta(shape)
        self.ring = ring("A B x1 x2 x3 x4", sp.ZZ)[0]
        self.phases = defaultdict(int)
        for signs in product((1, -1), repeat=6):
            vertices = {i: [(s, [[0,1,2] if s==1 else [0,2,1]])] for i,s in enumerate(signs)}
            [(trace, ranks)] = ribbon.colour_weight(6, theta.edges, theta.ext, vertices).items()
            [(rank, weight)] = ranks.items()
            phase = sum(s * v for s, v in zip(signs, theta.vertices))
            for order in permutations(range(4)):
                if ribbon.canon_tau(tuple(order[i] for i in cycle) for cycle in trace) != word:
                    continue
                p = [momenta[i] for i in order]
                prefactor = det(p[0],p[1]) * det(p[2],p[3])
                if not prefactor:
                    continue
                brackets = [det(p[i],p[j]) for i,j in ((0,1),(0,2),(1,2))]
                terms = defaultdict(int)
                for exps, integer in phase.items():
                    terms[exps[3:]] += int(integer)*prod(x**a for x,a in zip(brackets,exps[:3]))
                key = tuple(sorted((powers,c) for powers,c in terms.items() if c))
                self.phases[rank,key] += int(weight)*ribbon.perm_sign(order)*prefactor
        self.phases = {k:v for k,v in self.phases.items() if v}
        self.moments = {}

    def integrate(self, degrees):
        """Keep only one ribbon polynomial power in memory at a time."""
        degrees = sorted(degrees)
        self.moments = {degree: defaultdict(F) for degree in degrees}
        tables = {degree: DirichletMoments(self.theta.shape, degree) for degree in degrees}
        for (rank,key), weight in self.phases.items():
            phase = self.ring.from_dict(dict(key))
            power = self.ring.one
            last = 0
            for degree in degrees:
                power *= phase ** (degree-last)
                last = degree
                moment = tables[degree].integrate(power)
                scale = self.theta.prefactor / (4 * 2**degree * factorial(degree))
                self.moments[degree][rank] += scale * weight * moment


@lru_cache(maxsize=None)
def tree_terms(word, momenta, left, right):
    result = []
    for i,j in combinations(range(4),2):
        k,l = [s for s in range(4) if s not in (i,j)]
        rename = {s:(t,) for s,t in enumerate((i,j,k,l))}
        weight = 0
        for a,ca in COLOURS[left]:
            for b,cb in COLOURS[right]:
                for trace,col in ce._contract(a,3,b,3,False):
                    # The left vertex is single trace and cannot close an index loop at this edge.
                    assert col == {0:F(1)}
                    if ce.relabel(trace,rename) == word:
                        weight += ca*cb
        if weight:
            result.append(((-1)**(i+j-1)*weight,det(momenta[i],momenta[j]),det(momenta[k],momenta[l])))
    return result


def residual(n, word, probes, terms):
    result = defaultdict(F)
    for probe in probes[word]:
        for rank,q in probe.moments[n].items():
            result[rank] += q
    momenta = PAIRS if word == TRIPLE else CYCLE
    for j, deformations in terms.items():
        left_degree = n-j
        left = "d" if left_degree%2 else "f"
        scale = F(1,2**left_degree*factorial(left_degree))
        for right, ranks in deformations.items():
            value = sum(sign*u**left_degree*v**(j+2)
                        for sign,u,v in tree_terms(word,momenta,left,right))
            for rank,c in ranks.items():
                result[rank] += scale*value*c
    return {rank:c for rank,c in result.items() if c}


def initial_terms():
    terms = {}
    for n,cochain in zip((-1,0,1), laurent_terms()):
        by_word = {((0,1,2),): "d" if n%2 else "f", ((0,),(1,2)):"s1", ((0,),(1,),(2,)):"s0"}
        terms[n] = {}
        # At p0=(1,0), p1=(0,1) each cubic determinant equals one.
        for (word,powers),ranks in cochain.items():
            if powers[1] or powers[2] or word not in by_word:
                continue
            terms[n][by_word[word]] = ranks
    return terms


def coefficients(max_order):
    probes = {word: [RibbonProbe(shape,momenta,word) for shape in ((1,3,3),(2,2,3))]
              for word,momenta in ((SINGLE,CYCLE),(DOUBLE,CYCLE),(TRIPLE,PAIRS))}
    for word, instances in probes.items():
        degrees = [n for n in range(2, max_order+1)
                   if word == SINGLE or (word == DOUBLE) == bool(n%2)]
        for probe in instances:
            probe.integrate(degrees)
    terms = initial_terms()
    for n in range(2,max_order+1):
        degree = n+2
        single = residual(n,SINGLE,probes,terms)
        denom = 2**degree + (-1)**degree - 1
        answer = {"d" if n%2 else "f": {rank:c/denom for rank,c in single.items()}}
        if n%2:
            double = residual(n,DOUBLE,probes,terms)
            answer["s1"] = {rank:c/(2**degree-2) for rank,c in double.items()}
        else:
            triple = residual(n,TRIPLE,probes,terms)
            terms[n-1]["s0"] = {rank:-c for rank,c in triple.items()}
        terms[n] = answer
    return terms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("max_order", type=int)
    args = parser.parse_args()
    if args.max_order < 2:
        parser.error("max_order must be at least 2")
    for n,terms in coefficients(args.max_order).items():
        print(n, " ".join(f"{name} = " + " + ".join(f"({c}) N^{r}" for r,c in sorted(poly.items()))
                           for name,poly in terms.items()))


if __name__ == "__main__":
    main()
