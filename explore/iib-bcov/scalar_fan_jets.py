"""Normal derivatives of marked scalar letters in the fivefold and twofold fans."""
from fractions import Fraction
from itertools import permutations, product
from math import factorial, prod

from point_graph import solve_graph
from scalar_fan import chart, closed, fan, independent_integral, trees


def beta_moment(a, b, left, right):
    return Fraction(factorial(a + left - 1) * factorial(b + right - 1) * factorial(a + b - 1),
                    factorial(a - 1) * factorial(b - 1) * factorial(a + b + left + right - 1))


def ordered_moment(derivatives, dimension):
    prefix = derivatives[0]
    result = Fraction(1)
    for i, count in enumerate(derivatives[1:], 1):
        result *= beta_moment(dimension*i, dimension+1, prefix, count)
        prefix += count
    return result


def simplex_moment(derivatives, dimension):
    mass = dimension * len(derivatives)
    factors = prod(factorial(dimension + k - 1) // factorial(dimension - 1) for k in derivatives)
    return Fraction(factorial(mass - 1) * factors, factorial(mass + sum(derivatives) - 1))


def check_subdivided_currents(n):
    edges, cycles = fan(n)
    cuts = trees(cycles)
    lengths, _ = chart([Fraction(i+2, 2*i+5) for i in range(n)],
                       [Fraction(i+3, 3*i+7) for i in range(n)])
    shares = [Fraction(1)]
    for i in range(n):
        b = Fraction(i+3, 3*i+7)
        shares = [b*value for value in shares] + [1-b]
    reduced = [(-1, -1), (1, 0)] + [(0, 0)]*(n-1) + [(0, 1)]
    reduced_currents, _ = solve_graph(edges, reduced, lengths, 0)
    choices = product(cuts, repeat=3) if n <= 2 else (
        tuple(cuts[(i + j*(len(cuts)//3)) % len(cuts)] for j in range(3))
        for i in range(len(cuts)))
    checked = 0
    for selection in choices:
        sources = [[0]*5 for _ in range(n+2)]
        sources[1][3] = sources[n+1][4] = 1
        graph = []
        resistances = []
        projected = []
        for edge_id, (start, end) in enumerate(edges):
            sites = [start]
            for color, cut in enumerate(selection):
                if edge_id in cut:
                    sites.append(len(sources))
                    sources.append([int(j == color) for j in range(5)])
            sites.append(end)
            if edge_id <= n:
                projected.append((len(graph) + len(sites) - 2, edge_id))
            denominator = sum(range(1, len(sites)))
            for weight, (start, end) in enumerate(zip(sites, sites[1:]), 1):
                graph.append((start, end))
                resistances.append(lengths[edge_id] * weight / denominator)
        sources[0] = [-sum(column) for column in zip(*sources[1:])]
        currents, _ = solve_graph(graph, sources, resistances, 0)
        for e, spoke in projected:
            assert graph[e][1] == 0 and currents[e][3:] == reduced_currents[spoke]
            assert currents[e][4] == shares[spoke]
        checked += 1
    return checked


def check():
    print('Direct five-dimensional rooted current checks:',
          [check_subdivided_currents(n) for n in (1, 2, 3)])
    for n in range(1, 6):
        # In the two retained directions the bulk l density, including the
        # endpoint voltage, is the boundary l density pulled back by l -> l^2.
        ell = [Fraction(i+2, 2*i+5) for i in range(n)]
        b = [Fraction(i+3, 3*i+7) for i in range(n)]
        bulk_density = prod(ell[i]**(2*(n-1-i)+1) * b[i]**(2*i+1) * (1-b[i])**2
                            for i in range(n))
        boundary_density = prod((ell[i]**2)**(n-1-i) * b[i]**i * (1-b[i]) for i in range(n))
        share_product = prod(b[i]**(i+1) * (1-b[i]) for i in range(n))
        assert bulk_density == boundary_density * share_product * prod(ell)
        fixtures = [(0,)*(n+1), (3,) + (0,)*n, (2, 1) + (0,)*(n-1),
                    tuple((i*i+n*i+2) % 5 for i in range(n+1))]
        for degrees in fixtures:
            orders = set(permutations(degrees))
            for dimension in (1, 2):
                actual = sum(ordered_moment(order, dimension) for order in orders) / len(orders)
                assert actual == simplex_moment(degrees, dimension)
            raised = tuple(k+1 for k in degrees)
            bridge = Fraction(factorial(2*n+1), factorial(n))
            assert simplex_moment(degrees, 2) == bridge * simplex_moment(raised, 1)
            # Three ordinary repeated scalar amplitudes have Taylor divisor
            # (n!)^3. The e4 polarization measures only the endpoint voltage.
            bulk = (-1)**n * closed(n) * simplex_moment(degrees, 2) / ((3*n+1)*factorial(n)**3)
            paper = simplex_moment(raised, 1) / (factorial(n)**3 * (n+1))
            assert bulk == Fraction(-1, 2)**n * paper
        assert independent_integral(n, cotree_colors=0) == closed(n) / ((3*n+1)*factorial(n)**3)
        twice = (2,) + (0,)*n
        boundary = (n+1)*simplex_moment(twice, 1)
        bulk = (n+1)*simplex_moment(twice, 2)
        assert boundary - bulk == Fraction(n, (n+2)*(2*n+3))
        print(n, 'sum over the twice-differentiated letter, bulk and boundary:', bulk, boundary)


if __name__ == '__main__':
    check()
