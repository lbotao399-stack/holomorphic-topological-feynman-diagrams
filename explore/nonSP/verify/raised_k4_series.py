"""Exact periods of iterated cotree subdivisions of the three-loop K4 case 0.

`graph_k4_series.py 0:0 0:1` regenerates the two input polynomials from Kirchhoff
currents. This replay derives each inserted factor from the fixed spanning tree,
checks the canonical lifted graphs directly, and evaluates both star charts.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, prod
from pathlib import Path
from random import Random
import sys

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(BASE.parent / 'ht-families' / 'verify'))
from dimension_raise import split_times, subdivide
from htk import brute
from htk.linalg import det
from k4exact import K4Case
from reduced_k4_series import beta_finite_parts, convergence, one_dimension, read_homogeneous


COTREE = (3, 6, 8)


def tree_weight(terms, layers, opposite):
    """Expand (g_A r_i prod_{j != i}(1-r_j))**layers in the verified numerator."""
    out = defaultdict(F)
    others = [i for i in range(3) if i != opposite]
    for (a, b, c, *ratios), coef in terms.items():
        for j in range(layers+1):
            for k in range(layers+1):
                powers = ratios.copy()
                powers[opposite] += layers
                powers[others[0]] += j
                powers[others[1]] += k
                out[(a+layers, b, c, *powers)] += (
                    coef*(-1)**(j+k)*comb(layers, j)*comb(layers, k))
    return {k:v for k,v in out.items() if v}


def check_sparsity(graph):
    slacks = []
    for size in range(2, graph.nv+1):
        for vertices in combinations(range(graph.nv), size):
            count = sum(u in vertices and v in vertices for u, v in graph.edges)
            if count:
                slacks.append(4*size-5-3*count)
    assert (len(slacks), min(slacks), slacks.count(0)) == (1874, 0, 14)


def edge_check(case, layers, opposite):
    rng = Random(719+case.Dc)
    for trial in range(3):
        point = case.random_point(rng, den=101)
        _, times = case.g_ref(point)
        original = case.G
        indices = COTREE
        sources = [case.lam.get(v, [0]*3) for v in range(original.root)]
        g_a, g_b, *r = point[-5:]
        q = 1-sum(g*s for g,s in zip([g_a,g_b,1-g_a-g_b],r))
        weight = g_a*r[opposite]*prod(1-r[j] for j in range(3) if j != opposite)/q
        U = brute.symanzik_U(original, times)
        assert prod(times[i] for i in indices)/U == weight
        density = brute.density(original, times, sources).rational()
        assert density
        for j in range(layers):
            graph = subdivide(original, indices)
            if j == trial == 0:
                check_sparsity(graph)
            minor = det([[row[i] for i in indices] for row in original.cycle_matrix()])
            assert minor == (-1)**(j+1)
            ratios = [F(rng.randrange(1, 23), 23) for _ in indices]
            times = split_times(times, indices, ratios)
            assert brute.symanzik_U(graph, times) == U
            sources = [row+[0] for row in sources]+[[0]*(3+j)+[1] for _ in indices]
            new_density = brute.density(graph, times, sources).rational()
            assert new_density*U == -density
            indices = tuple(range(original.ne, graph.ne))
            original, density = graph, new_density


def run(layers):
    assert layers > 0
    results = []
    for center_index, opposite in enumerate((1, 2)):
        case = K4Case(0, Dc_index=center_index)
        edge_check(case, layers, opposite)
        degree, parent = read_homogeneous(BASE/'out'/f'case0_D{center_index}_P.m')
        numerator = tree_weight(parent, layers, opposite)
        degree += layers
        power = 10+layers
        margin = convergence(numerator, degree, power)
        raw = beta_finite_parts(one_dimension(degree, numerator, power))
        answer = {k:v*(-1)**layers/factorial(layers)**3 for k,v in raw.items()}
        results.append(answer)
        print('H', 3+layers, 'star center', case.Dc, 'numerator terms', len(numerator),
              'absolute margin', margin, 'coefficients of 1,zeta(k)', answer, flush=True)
    assert results[0] == results[1]
    if layers == 1:
        assert results[0] == {0:-F(3119, 83980800), 3:F(1, 32400)}
    return results[0]


if __name__ == '__main__':
    for arg in sys.argv[1:] or ['1']:
        run(int(arg))
