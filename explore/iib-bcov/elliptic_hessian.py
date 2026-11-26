"""Direct finite-jet polarization of the genus-one Hessian action and its BV divergence."""
from fractions import Fraction
from math import factorial, prod

from genus_one_graph_ward import circle_groups, classical, graph_vertex, plus, square


def add(*terms):
    out = {}
    for term in terms:
        for key,value in term.items():
            out[key] = out.get(key,0)+value
    return {k:v for k,v in out.items() if v}


def scale(series, factor):
    return {k:factor*v for k,v in series.items() if factor*v}


def product(left, right):
    # Each generator records one labeled wave and is square zero.
    out = {}
    for i, x in left.items():
        for j, y in right.items():
            if not i&j:
                out[i|j] = out.get(i|j,0) + x*y
    return {k:v for k,v in out.items() if v}


def determinant(rows, multiply=product, unit=None):
    # Laplace expansion by successive rows, retaining each minor once.
    minors = {0:{0:1} if unit is None else unit}
    for row in rows:
        next_minors = {}
        for columns,value in minors.items():
            for j,entry in enumerate(row):
                if columns & (1<<j):
                    continue
                term = multiply(value,entry)
                key = columns | (1<<j)
                sign = (-1)**((columns >> (j+1)).bit_count())
                next_minors[key] = add(next_minors.get(key,{}),scale(term,sign))
        minors = next_minors
    return minors[(1<<len(rows))-1]


def field_series(levels):
    """The formal solution U=T(U), H=(1-T'(U))^{-1}, log H, and T''(U)."""
    n = len(levels)
    fields = [{1<<i:1 for i,k in enumerate(levels) if k==j}
              for j in range(max(levels)+1)]
    def potential_derivative(u, order):
        result = {}
        for j in reversed(range(order,len(fields))):
            result = add(product(result,u),scale(fields[j],Fraction(1,factorial(j-order))))
        return result
    u = {}
    for _ in range(n):
        new = potential_derivative(u,0)
        if u==new:
            break
        u = new
    assert u==potential_derivative(u,0)
    prime = potential_derivative(u,1)
    h, log_h, power = {0:1}, {}, {0:1}
    for k in range(1,n+1):
        power = product(power,prime)
        if not power:
            break
        h = add(h,power)
        log_h = add(log_h,scale(power,Fraction(1,k)))
    return u,h,log_h,potential_derivative(u,2)


def momenta_by_subset(momenta):
    dim = len(momenta[0])
    n = len(momenta)
    sums = {0:(0,)*dim}
    for mask in range(1,1<<n):
        bit = mask & -mask
        sums[mask] = plus((sums[mask^bit],momenta[bit.bit_length()-1]))
    return sums


def hessian(series, sums):
    dim = len(sums[0])
    return [[{key:value*sums[key][i]*sums[key][j] for key,value in series.items()
              if sums[key][i] and sums[key][j]} for j in range(dim)] for i in range(dim)]


def cofactor_contraction(matrix, direction):
    result = {}
    for i in range(len(matrix)):
        rows = [direction[j] if j==i else row for j,row in enumerate(matrix)]
        result = add(result,determinant(rows))
    return result


def local_symbols(momenta, levels, variation=False, source_level=0):
    """Compute 24 times the indicated Taylor coefficient in the original A_k variables."""
    n = len(levels)
    u,h,log_h,second = field_series(levels)
    sums = momenta_by_subset(momenta)
    matrix = hessian(u,sums)
    det = determinant(matrix)
    if not variation:
        return product(log_h,det).get((1<<n)-1,0)
    first = product(product(h,h),product(second,det))
    cofactor = cofactor_contraction(matrix,hessian(log_h,sums))
    derivative = add(first,product(h,cofactor))
    if source_level:
        preceding = {0:1}
        for _ in range(source_level-1):
            preceding = product(preceding,u)
        indirect = scale(product(product(preceding,u),derivative),Fraction(1,factorial(source_level)))
        direct = scale(product(preceding,product(h,det)),Fraction(1,factorial(source_level-1)))
        derivative = add(indirect,direct)
    return derivative.get((1<<n)-1,0)


def point_square(momenta, levels, vector_slot, polarization):
    """The rooted hexagon on arbitrary descendant waves by two independent calculations."""
    labels = tuple(range(len(levels)))
    all_labels = (1<<len(labels))-1
    u = field_series(levels)[0]
    scalars = {mask:v for mask,v in u.items() if not mask&(1<<vector_slot)}
    vector = {mask:v for mask,v in u.items() if mask&(1<<vector_slot)}
    sums = momenta_by_subset(momenta)
    def pairing(p):
        return sum(x*y for x,y in zip(p,polarization))
    matrix = hessian(scalars,sums)
    total_div = pairing(sums[all_labels])*cofactor_contraction(
        matrix,hessian(vector,sums)).get(all_labels,0)
    div_w = {mask:v*pairing(sums[mask]) for mask,v in vector.items()}
    local = total_div+cofactor_contraction(matrix,hessian(div_w,sums)).get(all_labels,0)

    graph = 0
    for partition in circle_groups(labels,levels,5):
        charges = tuple(plus(momenta[i] for i in group) for group in partition)
        slot = next(i for i,group in enumerate(partition) if vector_slot in group)
        voltage = pairing(plus((*charges,charges[slot])))
        graph += square(charges)*voltage*prod(classical(levels,group) for group in partition)
    assert local==graph
    return Fraction(local,24)


def check():
    from random import Random
    basis = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    directed = plus(basis[:2])
    points = (*basis,directed,tuple(-x for x in plus((*basis,directed))))
    for degrees, expected in [((0,0,0,0,0,2,0),Fraction(1,2)),
                              ((0,0,0,0,0,1,1),Fraction(7,12))]:
        action = local_symbols(points,degrees)/24
        gradient = local_symbols(points[1:],degrees[1:],variation=True)/24
        assert expected==graph_vertex(points,degrees)==action==gradient
        # Q_0 phi=div(phi mu) fixes the divergence sign (+ sum of incoming momenta).
        vector = basis[0]
        contraction = sum(x*y for x,y in zip(plus(points[1:]),vector))
        assert contraction==-1
        one_edge = contraction*gradient
        # For genus one the source trace contributes -1/(2*pi*i)^5.
        assert -one_edge==expected
        print('Seven-wave fixture: action',action,'edge',one_edge,'heat numerator',-one_edge)
    rng = Random(2601)
    for levels in ((0,0,0,0,0,1),(0,0,0,0,1,1,1,0),(0,0,0,0,1,0,2,0),
                   (0,0,0,0,0,0,0,3),(0,2,0,1,0,0,0,0,1)):
        waves = tuple(tuple(rng.randrange(-3,4) for j in range(5)) for i in levels[1:])
        root = tuple(-x for x in plus(waves))
        points = (root,*waves)
        graph = graph_vertex(points,levels)
        assert graph and local_symbols(points,levels)==24*graph
        assert local_symbols(waves,levels[1:],variation=True)==24*graph
    for r in range(1,5):
        # A single level-r source and five active primary Hessian waves.
        waves = (*basis,*((0,)*5 for _ in range(r-1)))
        root = tuple(-x for x in plus(waves))
        action = local_symbols((root,*waves),(r,)+(0,)*len(waves))
        assert action==1==local_symbols(waves,(0,)*len(waves),variation=True,source_level=r)
        # Move a positive-level external wave to retain the indirect Euler derivative.
        degrees = (0,)*len(waves)+(1,)
        waves = (*waves,(1,2,3,4,5))
        root = tuple(-x for x in plus(waves))
        action = local_symbols((root,*waves),(r,*degrees))
        assert action and action==local_symbols(waves,degrees,variation=True,source_level=r)
    print('The complete Euler derivatives at source levels one through four pass direct polarization.')
    assert point_square(basis,(0,)*5,4,basis[0])==Fraction(1,24)
    for levels in ((0,0,0,0,1,0),(1,1,0,0,0,0,0),(0,0,0,0,0,2,0),
                   (0,0,0,0,0,0,0,3),(0,0,0,0,2,0,0,1)):
        waves = tuple(tuple(rng.randrange(-3,4) for j in range(5)) for i in levels)
        for slot in (0,len(levels)-1):
            assert point_square(waves,levels,slot,(1,-1,2,1,3))
    print('Rooted hexagon and local second-order point operator agree on moving descendants.')
    print('Direct compact action and variational Hessians agree with the graph recursion through nine waves.')


if __name__=='__main__':
    check()
