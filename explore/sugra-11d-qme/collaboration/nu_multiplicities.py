"""Exact labeled nu-attachment sums from finitely many polarized graph evaluations."""
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod, perm


def rising(base, degree):
    return prod(range(base, base+degree))


def falling_moment(vertices, multiplicities, orders):
    """Average of a product of falling powers of occupations (species, vertex)."""
    per_mode = [0]*len(multiplicities)
    per_vertex = [0]*vertices
    for (mode, vertex), order in orders.items():
        if not 0<=mode<len(per_mode) or not 0<=vertex<vertices or order<0:
            raise ValueError("invalid species, vertex, or falling power")
        per_mode[mode] += order
        per_vertex[vertex] += order
    return Fraction(prod(perm(k,r) for k,r in zip(multiplicities,per_mode))
                    *prod(factorial(r) for r in per_vertex),
                    rising(vertices,sum(per_vertex)))


def continue_multiplicities(vertices, multiplicities, samples):
    """Evaluate the raw labeled coefficient. samples maps multi-indices to raw values
    on the full total-degree simplex required by the projected kernel's degree bound.
    Newton interpolation uses the source factorials and accepts outputs in common units.
    """
    if vertices<1 or not samples or not multiplicities or min(multiplicities)<0:
        raise ValueError("invalid vertex count, samples, or multiplicities")
    dimension = len(multiplicities)
    if any(len(k)!=dimension or min(k)<0 for k in samples):
        raise ValueError("invalid sample index")
    degree = max(map(sum,samples))
    if len(samples)!=comb(degree+dimension, dimension):
        raise ValueError("incomplete total-degree simplex")
    target = Fraction(0)
    for order in samples:
        weight = prod(comb(k,r) for k,r in zip(multiplicities,order))
        if weight==0:
            continue
        difference = Fraction(0)
        for point in product(*(range(r+1) for r in order)):
            value = Fraction(samples[point],rising(vertices,sum(point)))
            coefficient = prod(comb(r,i) for r,i in zip(order,point))
            difference += (-1)**(sum(order)-sum(point))*coefficient*value
        target += weight*difference
    return target*rising(vertices,sum(multiplicities))
