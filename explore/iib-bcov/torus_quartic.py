"""Exact three-loop primary symbol on an oriented H^1(T^4) quartic."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import factorial, lcm, prod

from bulk_three import forest_data
from square_three import G, Graph, plan_generic, sp_monomial
from current_minor import numerator
from scalar_fan import cofactor
from point_graph import solve_graph


CORE = ((0,1), (0,1), (0,2), (0,2), (1,2))
FACTORIAL = tuple(factorial(i) for i in range(16))
# Five path denominators, two parallel reductions and one final triangle.
COMMON_DENOMINATOR = FACTORIAL[9]**8


def shapes():
    for counts in product(range(6), repeat=5):
        left, right = counts[:2], counts[2:4]
        if sum(counts) != 5 or left > right:
            continue
        if left[0] <= left[1] and right[0] <= right[1] and min(sum(left), sum(right)) >= 2:
            yield counts


@cache
def geometry(counts, core=CORE):
    edges, paths = [], []
    vertex = 1 + max(max(edge) for edge in core)
    for (start, end), count in zip(core, counts):
        nodes = (start, *range(vertex, vertex+count), end)
        paths.append(range(len(edges), len(edges)+count+1))
        edges.extend(zip(nodes, nodes[1:]))
        vertex += count
    graph = G(vertex, edges)
    _, parent = graph._tree()
    _, kirchhoff, separating_forests = forest_data(edges, vertex)
    return graph, parent, paths, kirchhoff, separating_forests


@cache
def path_weight(powers):
    a,b,c,d,e = powers
    left, right = a+b-10, c+d-10
    assert left+right+e == 10 and min(left, right, e, 10-max(a,b,c,d)) > 0
    moment = Fraction(FACTORIAL[left-1]*FACTORIAL[right-1]*FACTORIAL[e-1], FACTORIAL[9])
    for x,y in ((a,b), (c,d)):
        moment *= Fraction(FACTORIAL[9-x]*FACTORIAL[9-y], FACTORIAL[19-x-y])
    # Independent reduction on all five retained powers before expanding any path.
    assert moment == sp_monomial(tuple(k-1 for k in powers), 10, core_plan())
    moment /= prod(FACTORIAL[k-1] for k in powers)
    assert COMMON_DENOMINATOR % moment.denominator == 0
    return moment.numerator*(COMMON_DENOMINATOR//moment.denominator)


@cache
def core_plan():
    return plan_generic(Graph(3, CORE))


def integrate(polynomial, paths):
    total = 0
    for exponent, coefficient in polynomial.terms():
        powers = tuple(sum(int(exponent[i])+1 for i in path) for path in paths)
        total += int(coefficient)*prod(FACTORIAL[int(k)] for k in exponent)*path_weight(powers)
    return Fraction(total, COMMON_DENOMINATOR)


def tree_flow(graph, parent, sources):
    flow = [[0,0,0] for _ in graph.edges]
    for vertex, source in enumerate(sources):
        while parent[vertex] is not None:
            neighbor, edge = parent[vertex]
            sign = 1 if graph.edges[edge] == (vertex, neighbor) else -1
            flow[edge] = [x+sign*y for x,y in zip(flow[edge], source)]
            vertex = neighbor
    return flow


def shape_integral(counts, waves, polarization, core=CORE):
    graph, parent, paths, kirchhoff, forests = geometry(counts, core)
    branch_count = graph.nv - len(waves)
    root = tuple(-sum(values) for values in zip(*waves))
    jets = tuple((p, i == len(waves)-1) for i,p in enumerate(waves))
    total = None
    tested = False
    times = tuple(2+(i*i+2*i)%11 for i in range(graph.N))
    for ordered, multiplicity in Counter(permutations(jets)).items():
        sources = [root, *([(0,0,0)]*(branch_count-1)), *(p for p,_ in ordered)]
        at = branch_count + next(i for i,(_,vector) in enumerate(ordered) if vector)
        current = numerator(graph, sources, tree_flow(graph, parent, sources))
        if not current:
            continue
        charges = [sum(x*y for x,y in zip(p, polarization)) for p in sources]
        voltage = sum(q*(forests[at][0]+forests[j][0]-forests[at][j])
                      for j,q in enumerate(charges))/2
        assert all(c.denominator == 1 for _,c in voltage.terms())
        voltage = current.context().from_dict({m:int(c) for m,c in voltage.terms()})
        if not tested:
            flows, potentials = solve_graph(graph.edges, sources, times, 0)
            u = int(kirchhoff(*times))
            assert abs(int(current(*times))) == u**graph.n*cofactor(flows, graph.C, times)
            assert int(voltage(*times)) == u*sum(x*y for x,y in zip(potentials[at], polarization))
            tested = True
        contribution = multiplicity*current**2*voltage
        total = contribution if total is None else total+contribution
    if total is None:
        return Fraction(0)
    if graph.n == 2:
        plan = plan_generic(Graph(graph.nv, graph.edges))
        value = sum(int(c)*sp_monomial(tuple(map(int,m)), 8, plan) for m,c in total.terms())
        return value/prod(factorial(m) for m in Counter(counts).values())
    automorphisms = (1+(counts[0]==counts[1]))*(1+(counts[2]==counts[3]))*(1+(counts[:2]==counts[2:4]))
    return integrate(total, paths)/automorphisms


def quartic(waves, polarization):
    """Four labeled scalar momenta followed by the transverse vector momentum; unit outgoing colors."""
    if len(waves) != 5 or len(polarization) != 3 or any(len(p) != 3 for p in waves):
        raise ValueError('Expected five three-dimensional momenta and one polarization.')
    waves = tuple(tuple(Fraction(x) for x in p) for p in waves)
    polarization = tuple(Fraction(x) for x in polarization)
    if sum(x*y for x,y in zip(waves[-1], polarization)):
        raise ValueError('The vector polarization must be transverse.')
    scale = lcm(*(x.denominator for p in waves for x in p))
    vector_scale = lcm(*(x.denominator for x in polarization))
    waves = tuple(tuple(int(scale*x) for x in p) for p in waves)
    vector = tuple(int(vector_scale*x) for x in polarization)
    return {counts: 24*shape_integral(counts, waves, vector)/(scale**19*vector_scale)
            for counts in shapes()}


if __name__ == '__main__':
    collinear = ((1,0,0),)*3+((0,1,0),(0,0,1))
    conserved = ((2,1,0),(-1,0,2),(1,-2,-1),(0,1,-1),(-2,0,0))
    assert all(sum(column) == 0 for column in zip(*conserved))
    for waves, vector, expected in (
        (collinear, (1,1,0), Fraction(1,10080)),
        (conserved, (0,1,0), -Fraction(59263,189000)),
    ):
        values = quartic(waves, vector)
        assert len(values) == 7 and sum(values.values()) == expected
        print('Ordered quartic', expected, flush=True)
        for counts, value in values.items():
            print(counts, value)
