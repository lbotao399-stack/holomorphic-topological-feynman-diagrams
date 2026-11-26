"""Highest scalar symbol with 3n labeled momenta in a three-plane."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import factorial, lcm, prod

from flint import fmpq_mat, fmpz_mat

from marked_theta import partitions, square, conductance, green_moments
from scalar_fan import chart, closed, cofactor, fan, solve_graph, trees


def volume(rows):
    return F(str(fmpq_mat([[str(x) for x in row] for row in rows]).det()))


def triple_sum(waves, scalar, vector):
    """Product of squared five-volumes, summed over unordered triple partitions."""
    @cache
    def squared_volume(p, q, r):
        return volume([p, q, r, scalar, vector])**2

    @cache
    def partition(slots):
        if not slots:
            return 1
        first, *rest = slots
        result = 0
        for j, k in combinations(rest, 2):
            value = squared_volume(waves[first], waves[j], waves[k])
            if value:
                remaining = tuple(i for i in rest if i != j and i != k)
                result += value * partition(remaining)
        return result

    return partition(tuple(range(len(waves))))


def reduced(waves, scalar, vector, polarization):
    waves = tuple(map(tuple, waves))
    if not waves or len(waves) % 3:
        raise ValueError("Expected a positive multiple of three plane waves.")
    if any(len(p)!=5 for p in (*waves,scalar,vector,polarization)):
        raise ValueError("Expected five-component covectors and polarization.")
    if fmpq_mat([[str(x) for x in p] for p in waves]).rank() > 3:
        raise ValueError("The selected scalar momenta must span at most three dimensions.")
    if sum(x*y for x, y in zip(vector, polarization)):
        raise ValueError("The vector wave must have zero divergence.")
    n = len(waves)//3
    momentum = tuple(map(sum, zip(scalar, *waves)))
    return (F(1, 2**n * (n+1) * factorial(2*n+1) * factorial(n)) *
            triple_sum(waves, scalar, vector) *
            sum(x*y for x, y in zip(momentum, polarization)))


def direct_two_loop(waves):
    """Enumerate all 5^6 assignments of six mixed waves to the reduced fan."""
    edges, cycles = fan(2)
    # The independent beta means integrate every multiaffine voltage after
    # the path sum. Keep each wave and each edge marginal separate until then.
    times, _ = chart([F(3,4), F(1,2)], [F(2,5), F(4,7)])
    _, phi = solve_graph(edges, [(-1,-1),(1,0),(0,0),(0,1)], times, 0)
    green = [phi[i][1] for i in range(4)]
    cuts = trees(cycles)
    u = sum(prod(times[e] for e in s) for s in cuts)
    weights = {}
    columns = [[[c[e]*x for c in cycles for x in p[:3]] for e in range(len(edges))]
               for p in waves]
    counts = [[F() for e in edges] for p in waves]
    total = F()
    for order in product(range(len(edges)), repeat=len(waves)):
        if max(order.count(e) for e in range(len(edges))) > 3:
            continue
        det = int(fmpz_mat([row[e] for row, e in zip(columns, order)]).det())
        if not det:
            continue
        key = tuple(sorted(order))
        if key not in weights:
            weights[key] = prod(times[e] for e in key)
        mass = det**2 * weights[key]
        total += mass
        for marginal, e in zip(counts, order):
            marginal[e] += mass
    scalar, vector = (0,0,0,1,0), (0,0,0,0,1)
    assert total == u**3 * 2 * triple_sum(waves, scalar, vector)
    for marginal in counts:
        for e in range(len(edges)):
            cotrees = sum(prod(times[j] for j in s) for s in cuts if e in s)
            assert 2*u*marginal[e] == total*cotrees
    mass = F(2**2*factorial(2), factorial(4)*factorial(3)*factorial(5))
    response = []
    for i in range(4):
        voltage = total*scalar[i]*green[1]
        for p, marginal in zip(waves, counts):
            voltage += p[i]*sum(m*(green[a]+green[b])/2
                                for m, (a,b) in zip(marginal, edges))
        actual = mass*voltage/u**3
        v = tuple(int(j == i) for j in range(5))
        assert actual == reduced(waves, scalar, vector, v)
        response.append(actual)
    return response


def restored_graph(waves):
    """Compare compressed and original five-component determinants at unequal times."""
    edges, cycles = fan(2)
    times = [F(x) for x in (2,5,3,7,11)]
    assignment = (0,2,1,4,3,2)
    source = [(0,)*5, (0,0,0,1,0), (0,)*5, (0,0,0,0,1), *waves]
    source[0] = tuple(-sum(c) for c in zip(*source[1:]))
    child_edges, child_times, parents = [], [], []
    for e, (a,b) in enumerate(edges):
        vertices = [a,*[4+i for i,k in enumerate(assignment) if k==e],b]
        parts = list(range(1,len(vertices)))
        child_times.extend(times[e]*k/sum(parts) for k in parts)
        parents.extend([e]*len(parts))
        child_edges.extend(zip(vertices,vertices[1:]))
    full_current, _ = solve_graph(child_edges,source,child_times,0)
    pure_current, _ = solve_graph(edges,[(-1,-1),(1,0),(0,0),(0,1)],times,0)
    columns = [[c[e]*x for c in cycles for x in p[:3]] for p,e in zip(waves,assignment)]
    deleted = abs(int(fmpz_mat(columns).det()))
    assert deleted
    lifted = [[c[e] for e in parents] for c in cycles]
    assert cofactor(full_current,lifted,child_times) == deleted*cofactor(pure_current,cycles,times)


def theta_weights(momenta, polarization):
    """Full two-loop highest-order coefficient from labeled theta partitions."""
    if len(momenta)!=8 or any(len(p)!=5 for p in (*momenta,polarization)):
        raise ValueError("Expected seven fivefold scalar waves and one vector wave.")
    scale = lcm(*(F(x).denominator for p in momenta for x in p))
    momenta = tuple(tuple(int(F(x)*scale) for x in p) for p in momenta)
    vector_slot = 7
    current = [sum(x*y for x,y in zip(p,polarization)) for p in momenta]
    if current[vector_slot]:
        raise ValueError("The vector wave must have zero divergence.")

    @cache
    def global_voltage(powers):
        return sum(q*green_moments(powers,vector_slot,j)[0]
                   for j,q in enumerate(current) if q)

    @cache
    def path_voltage(powers,group):
        return sum(current[j]*green_moments(powers,vector_slot,j)[1]
                   for j in group if current[j])

    result = []
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        degrees = tuple(6-length for length in shape)
        accumulated = F()
        for groups in partitions(shape):
            h = next(i for i,g in enumerate(groups) if vector_slot in g)
            bridge = F(9,degrees[h])
            for powers,c in square(groups,momenta).items():
                accumulated += c*(global_voltage(powers)+bridge*path_voltage(powers,groups[h]))
        symmetry = prod(factorial(m) for m in Counter(shape).values())
        # Two ten-column determinants and one linear voltage.
        result.append(accumulated*conductance(degrees)/(symmetry*scale**21))
    return result


if __name__ == "__main__":
    basis = [tuple(int(j==i) for j in range(5)) for i in range(5)]
    for n in range(1,7):
        waves = [p for p in basis[:3] for _ in range(n)]
        assert reduced(waves,basis[3],basis[4],(1,1,1,1,0)) == closed(n)
    plane = [(1,0,1),(2,-1,0),(0,1,1),(-1,2,1),(1,2,0),(3,0,-1)]
    waves = [(*p,0,0) for p in plane]
    restored_graph(waves)
    print("Two-loop transverse components:", *direct_two_loop(waves))
    v = (1,-1,0,0,0)
    assert sum(theta_weights([*waves,basis[3],basis[4]],v)) == reduced(waves,basis[3],basis[4],v)
    generic = [*basis[:4], (1,1,1,1,1),(1,2,1,0,-1),(2,1,0,1,2),basis[4]]
    assert sum(p[0]-p[1] for p in generic) == 0
    result = theta_weights(generic,v)
    assert sum(result) == -F(91,8640)
    print("Five-dimensional scalar support, vanishing total momentum pairing:", *result, "sum", sum(result))
    generic[-1] = tuple(-sum(c) for c in zip(*generic[:-1]))
    result = theta_weights(generic,v)
    assert result == [F(733,181440),F(197,25920),F(5,1344),F(121,60480)]
    assert sum(result) == F(5,288)
    print("Zero total external momentum:", *result, "sum", sum(result))
