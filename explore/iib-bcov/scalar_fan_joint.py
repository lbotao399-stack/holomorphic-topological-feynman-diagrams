"""Exact mixed endpoint-current moments of the boundary and fivefold scalar fans."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod

from flint import fmpq_mat
from point_graph import solve_graph
from scalar_fan import chart, fan
from scalar_fan_jets import beta_moment


def moment(jets, h, *, symmetric=False):
    """Normalized moment; each pair gives the first- and last-endpoint jet orders."""
    jets = tuple(map(tuple, jets))
    if h not in (1, 2) or not jets or any(min(pair) < 0 for pair in jets):
        raise ValueError('Supply nonnegative jet orders and h=1 (boundary) or 2 (bulk).')
    counts = Counter(jets)
    types = tuple(counts)
    order = tuple(types.index(pair) for pair in jets)
    n = len(jets) - 1
    total_a = sum(a for a, b in jets)
    total_b = sum(b for a, b in jets)
    beta = lru_cache(None)(beta_moment)

    @lru_cache(None)
    def integrate(remaining, assigned):
        left = sum(remaining)
        j = len(jets) - left
        previous_a = total_a - sum(k*a for k, (a, b) in zip(remaining, types))
        previous_b = total_b - sum(k*b for k, (a, b) in zip(remaining, types))
        result = F()
        for t in (range(len(types)) if symmetric else (order[j],)):
            if not remaining[t]:
                continue
            a, b = types[t]
            weight = F(remaining[t], left) if symmetric else F(1)
            if j:
                weight *= beta(h*j, h+1, previous_a+previous_b-assigned, a+b)
            if j == n:
                result += weight
                continue
            rest = list(remaining)
            rest[t] -= 1
            available = previous_a+a-assigned
            for take in range(available+1):
                # These labeled first-endpoint currents start their return to
                # ground here. All unassigned first currents pass this cut.
                extend = comb(available, take) * beta(h*(n-j), 1, total_a-assigned-take, take)
                result += weight * extend * integrate(tuple(rest), assigned+take)
        return result

    return integrate(tuple(counts[t] for t in types), 0)


def ordered_density(a, b):
    return prod(sum(min(x*v, u*y) for u, v in zip(a, b)) for x, y in zip(a, b))


def ordered_uniform_moment(jets):
    """Independent chamber integral on ordered uniform points in the unit interval."""
    size = len(jets)
    polynomial = [F(1)]
    normalization = factorial(size-1)**2 * factorial(size)
    for p, q in reversed(jets):
        normalization *= factorial(1+p+q)
        primitive = [F()]*(len(polynomial)+p+q+1)
        for j, coefficient in enumerate(polynomial):
            for k in range(q+1):
                degree = j+p+k+1
                primitive[degree] += coefficient*(-1)**k*comb(q,k)/degree
        polynomial = primitive
    normalization = F(normalization, factorial(size-1+sum(p for p,q in jets)) *
                      factorial(size-1+sum(q for p,q in jets)))
    return normalization*sum(polynomial)


def inverse_check(n):
    ell = [F(i+1, 2*i+3) for i in range(n)]
    split = [F(i+3, 2*i+5) for i in range(n)]
    times, _ = chart(ell, split)
    edges, _ = fan(n)
    sources = [(-1, -1), (1, 0)] + [(0, 0)]*(n-1) + [(0, 1)]
    currents, _ = solve_graph(edges, sources, times, 0)
    a, b = zip(*currents[:n+1])
    H = [x*sum(b[:i]) + y*sum(a[i:]) for i, (x, y) in enumerate(zip(a, b))]
    minors = [a[i]*b[i+1]-a[i+1]*b[i] for i in range(n)]
    assert min(minors) > 0
    inverse = [1/x for x in H] + [minor/(H[i]*H[i+1]) for i, minor in enumerate(minors)]
    assert all(r*times[0] == t*inverse[0] for r, t in zip(inverse, times))
    assert prod(H) == ordered_density(a, b) == prod(b)*prod(ell[i]**(n-i) for i in range(n))


def integrate_network(n, h):
    """Integrate the grounded electrical currents by rational tensor quadrature."""
    nodes = (F(1, 4), F(1, 2), F(3, 4))
    powers = fmpq_mat([[str(x**k) for x in nodes] for k in range(3)])
    def rule(a, b):
        rhs = fmpq_mat([[str(beta_moment(a, b, k, 0))] for k in range(3)])
        return tuple(zip(nodes, (F(str(row[0])) for row in powers.solve(rhs).tolist())))
    rules = [rule(h*(n-i), 1) for i in range(n)] + [rule(h*(i+1), h+1) for i in range(n)]
    edges, _ = fan(n)
    sources = [(-1, -1), (1, 0)] + [(0, 0)]*(n-1) + [(0, 1)]
    result = F()
    for grid in product(*rules):
        x, weights = zip(*grid)
        times, _ = chart(x[:n], x[n:])
        currents, _ = solve_graph(edges, sources, times, 0)
        result += prod(weights)*sum(a*b for a, b in currents[:n+1])
    return result


def check():
    for n in range(1, 9):
        inverse_check(n)
        for m in range(7):
            jets = ((1, m),) + ((0, 0),)*n
            bulk = (n+1)*moment(jets, 2, symmetric=True)
            paper = (n+1)*moment(jets, 1, symmetric=True)
            assert bulk == prod(F(4*k*(k+1), (2*k+m)*(2*k+m+2)) for k in range(1, n+1))
            assert paper == F(1, comb(n+m, m))
            if m == 1 and n <= 3:
                assert (bulk, paper) == (integrate_network(n, 2), integrate_network(n, 1))
                print('same-spoke mixed derivative', n, bulk, paper)
    for p, q in product(range(6), repeat=2):
        assert moment(((p,q), (0,0)), 2, symmetric=True) == F(12, (p+2)*(q+2)*(p+q+3))
    for n in range(1, 5):
        jets = tuple((i % 2, (i+1) % 3) for i in range(n+1))
        assert moment(jets, 1) == ordered_uniform_moment(jets)
        # Expand the ordered density into boundary displacement derivatives.
        comparison = F()
        for chosen in product(range(n+1), repeat=n+1):
            shifted = [list(pair) for pair in jets]
            for i, k in enumerate(chosen):
                shifted[max(i, k)][0] += 1
                shifted[min(i, k)][1] += 1
            comparison += moment(shifted, 1)
        assert moment(jets, 2) == F(2**n*factorial(2*n+1), factorial(n))*comparison
    print('Ordered differential factors, projective inverse, simplex integrals and mixed moments passed.')


if __name__ == '__main__':
    check()
