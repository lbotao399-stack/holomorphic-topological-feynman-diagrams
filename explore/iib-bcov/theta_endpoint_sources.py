"""Integrated fivefold theta potential with a source at each trivalent endpoint."""

from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import factorial, lcm, prod

from flint import fmpz_mat

from bulk_eight import SHAPES, SLOTS, coeffs, determinants, grounded_potentials, template
from marked_theta import conductance, green_moments, partitions


@cache
def endpoint_template(shape):
    groups = next(partitions(shape))
    path = {i:h for h,group in enumerate(groups) for i in group}
    regular,_,_,_,splits = template(shape)
    # The extra coordinate s_8=1 changes L to p_8 + sum_i s_i p_i.
    pairs = [*regular,*((i,8) for i in SLOTS),(8,8)]
    products = list(combinations_with_replacement(range(len(pairs)),2))
    degrees = tuple(6-n for n in shape)
    mass = conductance(degrees)/(2*prod(factorial(n) for n in Counter(shape).values()))
    green_pairs = tuple(combinations_with_replacement(range(9),2))
    matrix = []
    denominator = 1
    for i,j in products:
        word = pairs[i]+pairs[j]
        powers = tuple(word.count(k) for k in SLOTS)
        cube = F(1,prod(n+1 for n in powers))
        row = []
        for a,b in green_pairs:
            if a==8:
                value = cube
            elif b==8:
                value = cube*F(powers[a]+1,powers[a]+2)
            else:
                value,bridge = green_moments(powers,a,b)
                if path[a]==path[b]:
                    value += F(9,degrees[path[a]])*bridge
            value *= mass*(1+(i!=j))
            denominator = lcm(denominator,value.denominator)
            row.append(value)
        matrix.append(row)
    integer = fmpz_mat([[int(x*denominator) for x in row] for row in matrix])
    return pairs,products,green_pairs,integer,denominator,splits


def endpoint_minors(waves):
    minors = determinants(waves[:8])
    for group,row in minors.items():
        row.append(int(fmpz_mat([waves[8],*(waves[i] for i in group)]).det()))
    return minors


def quadratic(order,splits,pairs,minors):
    q = coeffs(order,splits,pairs,minors)
    # The symmetric cross-term formula includes both orders; the diagonal
    # p_8 in the first block and p_8 in the second block occurs only once.
    assert pairs[-1]==(8,8) and q[-1]%2==0
    q[-1] //= 2
    return q


def theta_potentials(waves, endpoint, ground, shapes=SHAPES):
    """Return eight bivalent rows, then the ungrounded and grounded rows.

    Endpoint assignments are ordered; the mass uses the same factor 1/2 for
    path reversal as bulk_eight. Sum the assigned classical vertex factors
    over both placements when the two endpoint decorations are distinct.
    """
    charges = (*waves,endpoint,ground)
    if len(waves)!=8 or any(len(p)!=5 for p in charges) or any(map(sum,zip(*charges))):
        raise ValueError('Eight bivalent and two endpoint fivefold charges must be conserved.')
    scale = lcm(*(F(x).denominator for p in charges for x in p))
    p = [tuple(int(F(x)*scale) for x in q) for q in charges]
    minors = endpoint_minors(p)
    result = [[F() for _ in range(5)] for _ in range(10)]
    for shape in shapes:
        pairs,products,green_pairs,moments,den,splits = endpoint_template(shape)
        orders,values = [],[]
        for groups in partitions(shape):
            order = (*sum(groups,()),8)
            q = quadratic(order,splits,pairs,minors)
            if any(q):
                orders.append(order)
                values.append([q[i]*q[j] for i,j in products])
        if not orders:
            continue
        weighted = (fmpz_mat(values)*moments).tolist()
        total = [[0]*9 for _ in range(9)]
        for order,row in zip(orders,weighted):
            for (i,j),coefficient in zip(green_pairs,row):
                a,b = order[i],order[j]
                total[a][b] += int(coefficient)
                if a!=b:
                    total[b][a] += int(coefficient)
        for i,row in enumerate(total):
            for k in range(5):
                result[i][k] += F(sum(v*q[k] for v,q in zip(row,p)),den*scale**21)
    return tuple(map(tuple,result))


def fixed_metric(waves, endpoint, ground, shape):
    """Compare the Laplace solution and the 11-column current determinant."""
    from point_graph import solve_graph
    from scalar_fan import cofactor

    positions = tuple(F(i,11) for i in (2,8,1,5,7,4,9,3))
    resistance = (F(2),F(3),F(6))  # parallel resistance one
    pairs,_,_,_,_,splits = endpoint_template(shape)
    minors = endpoint_minors((*waves,endpoint))
    for groups in partitions(shape):
        order = (*sum(groups,()),8)
        q = quadratic(order,splits,pairs,minors)
        x = (*positions,1)
        value = sum(v*x[order[i]]*x[order[j]] for v,(i,j) in zip(q,pairs))
        if not value:
            continue
        edges,times,parents = [],[],[]
        for h,group in enumerate(groups):
            nodes = [9,*sorted(group,key=positions.__getitem__),8]
            for i,j in zip(nodes,nodes[1:]):
                edges.append((i,j))
                times.append(resistance[h]*(x[j]-(x[i] if i!=9 else 0)))
                parents.append(h)
        flow,phi = solve_graph(edges,[*waves,endpoint,ground],times,9)
        cycles = [[(int(h==j)-int(h==1)) for h in parents] for j in (0,2)]
        tree = sum(resistance[i]*resistance[j] for i,j in combinations(range(3),2))
        assert cofactor(flow,cycles,times)==abs(value)*prod(resistance)/tree**2
        end = [endpoint[k]+sum(t*p[k] for t,p in zip(positions,waves)) for k in range(5)]
        assert end==phi[8]
        for h,group in enumerate(groups):
            for i in group:
                s = positions[i]
                for k in range(5):
                    interior = sum((min(s,positions[j])-s*positions[j])*waves[j][k] for j in group)
                    assert phi[i][k] == s*end[k]+resistance[h]*interior
        return value
    raise AssertionError('The fixed-metric test requires a nonzero current cofactor.')


if __name__=='__main__':
    waves = tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
    waves += ((1,2,-1,3,2),(-2,1,3,1,-2),(3,-1,1,0,2))
    endpoint = (1,1,0,-1,2)
    ground = tuple(-sum(column) for column in zip(*waves,endpoint))
    zero = (0,)*5
    for shape in SHAPES:
        print('Fixed electrical determinant',shape,fixed_metric(waves,endpoint,ground,shape))
        first = theta_potentials(waves,endpoint,ground,(shape,))
        swapped = theta_potentials(waves,ground,endpoint,(shape,))
        assert any(first[8]) and first[9]==zero
        assert all(tuple(a-b for a,b in zip(first[i],first[8]))==swapped[i] for i in SLOTS)
        assert tuple(-x for x in first[8])==swapped[8]
    first = theta_potentials(waves,endpoint,ground)
    ground = tuple(-sum(column) for column in zip(*waves))
    assert theta_potentials(waves,zero,ground)[:8]==tuple(grounded_potentials(waves))
    waves = (*waves[:7],tuple(-sum(column) for column in zip(*waves[:7])))
    rational = [tuple(F(x,2) for x in p) for p in waves]
    ordinary = grounded_potentials(waves)
    assert theta_potentials(rational,zero,zero)[:8]==tuple(tuple(x/2**21 for x in row) for row in ordinary)
    print('All four shapes pass endpoint reversal and the grounded-source reduction.')
    print('Eight marked voltages and both endpoint voltages evaluated. First endpoint:',*first[8])
