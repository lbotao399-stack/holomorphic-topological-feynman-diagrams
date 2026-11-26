"""Genus-three fivefold trace/Feynman evaluation of the terminal planar-probe fusions."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json
import sys

IIB = Path(__file__).resolve().parents[2] / 'iib-bcov'
sys.path.insert(0, str(IIB))
from bulk_three_hexagon import potentials as hexagon
from elementary_scalar import coefficients_from_graphs
from genus_two_graph_ward import divergences
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers

BASIS = tuple(tuple(int(i == j) for j in range(5)) for i in range(5))


def layouts(first_count, second_count):
    """Nine labeled primaries enter the first hard pair, the second, or remain free."""
    for first in product(range(4), repeat=3):
        if sum(first) != first_count:
            continue
        for second in product(range(4), repeat=3):
            if sum(second) != second_count:
                continue
            free = tuple(3-i-j for i, j in zip(first, second))
            if min(free) < 0:
                continue
            weight = prod(comb(3, i)*comb(3-i, j) for i, j in zip(first, second))
            yield first, second, free, weight


def waves(first, second):
    """Unit-determinant shears put all transverse pair momenta on the positive planar legs."""
    free = tuple(3-i-j for i,j in zip(first, second))
    legs = [(1, 0, first[0]-4, first[1]-4, first[2]-3), (-1, 0, 0, 0, 0),
            (0, 1, *second), (0, -1, 0, 0, 0)]
    legs.extend(p for p, n in zip(BASIS[2:], free) for _ in range(n))
    legs.extend((BASIS[2], BASIS[3]))   # Absorbed root and transverse primary vector.
    assert len(legs) == 9 and not any(map(sum, zip(*legs)))
    return legs


def primary_weights(levels):
    split = comb(4, levels[0])*comb(2, levels[2])
    return {(a,b): split*w for a,b,_,w in layouts(4, 2)}


def one_level_weights(levels):
    """Four primary marks on e3 or e4 force the sole remaining psi into a rational tail."""
    terms = Counter()
    count = 0
    for distinguished, level in enumerate(levels):
        if not level:
            continue
        used = tuple(k-int(i == distinguished) for i,k in enumerate(levels))
        totals = sum(used[:2]), sum(used[2:])
        split = comb(totals[0], used[0])*comb(totals[1], used[2])
        for a, b, free, weight in layouts(*totals):
            for axis in (0, 1):
                if free[axis] != 3:
                    continue
                final = [list(a), list(b)]
                final[distinguished//2][axis] += 1
                terms[tuple(map(tuple, final))] += 4*split*weight
                count += split*weight
    assert sum(terms.values()) == 4*count
    return count, terms


def period(weights, values):
    return tuple(sum(w*values[k][j] for k,w in weights.items()) for j in (0,1))



def tree_fixture():
    # The third hard scalar still carries A1. The four equal charges include the root.
    return [(1,0,-4,-4,-1), (-1,0,0,2,0), BASIS[1], (0,-1,0,1,0),
            BASIS[2], BASIS[2], BASIS[2], BASIS[4], BASIS[2], BASIS[3]]


def fused_trees(p):
    groups = {}
    for i,j in combinations(range(10),2):
        joined = tuple(x+y for x,y in zip(p[i],p[j]))
        groups[i,j] = tuple(sorted([joined]+[q for h,q in enumerate(p) if h not in (i,j)]))
    return groups


def feynman_ward(p, trees, coefficients):
    from point_graph import solve_graph
    from square_three import G
    from flint import fmpq_mat
    # Restore the ten bivalent source vertices on both cubic cores, at unequal times.
    cores = (tuple(combinations(range(4),2)), ((0,1),(0,1),(1,2),(2,3),(2,3),(3,0)))
    for core in cores:
        for shift in range(6):
            edges = []
            for e,(a,b) in enumerate(core):
                chain = [a,*[4+i for i in range(10) if (i+shift)%6==e],b]
                edges.extend(zip(chain,chain[1:]))
            times = [3+(j*j+3*j)%17 for j in range(16)]
            currents,_ = solve_graph(edges, [(0,)*5]*4+p, times, 0)
            cycles = G(14,edges).C
            columns = fmpq_mat([[str(currents[j][k]*c[j]) for j in range(16)]
                                for c in cycles for k in range(5)])
            assert columns.rank()<15
    assert not any(x for row in hexagon(p) for x in row)
    reference = coefficients[trees[2,4]]
    for ch in (0,1):
        potentials = [[Fraction()]*5 for q in p]
        for (i,j),key in trees.items():
            v = coefficients[key][ch]
            collinear = int(p[i]==BASIS[2])+int(p[j]==BASIS[2])
            assert v == (0,1,4)[collinear]*reference[ch]
            for h in (i,j):
                for k in range(5):
                    potentials[h][k] += v*(p[i][k]+p[j][k])
        assert divergences(p,potentials) == tuple((22 if q==BASIS[2] else 4)*reference[ch] for q in p)
    return tuple(4*x for x in reference)

def check():
    groups = tuple(primary_weights((2,2,1,1)))
    points = [waves(*g) for g in groups]
    # One unsheared, fully fused descendant-boundary term, with e3 fused into mark 3.
    actual = [(1,0,-4,-4,-1), (-1,0,0,2,0), (0,1,1,0,0), (0,-1,0,1,0),
              BASIS[2], BASIS[2], BASIS[4], BASIS[2], BASIS[3]]
    example = ((0,2,2), (1,1,0))
    fixture = tree_fixture()
    trees = fused_trees(fixture)
    cuts = list(dict.fromkeys(trees.values()))
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++', '-std=c++17', '-O3', str(IIB/'operator'/'nine_sextics.cpp'), '-o', str(exe)], check=True)
        rows = orbit_integers(exe, [*points, actual, *cuts])
    assert rows.pop(len(points)) == rows[groups.index(example)]
    # Independently solve from fifteen integrated Feynman/Ward rows. All these wave
    # configurations have either three on a line or two distinct collinear pairs.
    fitted = coefficients_from_graphs(rows).tolist()
    data = json.loads((IIB/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words'])) == DOUBLE_SOFT
    channels = tuple(zip(*data['rational_and_zeta']))
    values = {}
    for g, row, independent in zip([*groups, *cuts], rows, fitted, strict=True):
        reduced = tuple(sum(x*Fraction(y) for x,y in zip(row,col)) for col in channels)
        assert reduced == tuple(Fraction(str(x)) for x in independent)
        values[g] = tuple(-x for x in reduced)   # g=3, root(v)=1, V3=-C.
    for levels in ((2,2,1,1), (2,2,2,0)):
        terminal = primary_weights(levels)
        count, descendant = one_level_weights(levels)
        print(levels, 'six fusions:',sum(terminal.values()),period(terminal, values), flush=True)
        print(levels, 'five fusions with an overfull axis:',count,period(descendant, values), flush=True)
    single = feynman_ward(fixture, trees, values)
    assert single == tuple(4*x for x in values[example]) == (Fraction(11180467393,23538138624000), -Fraction(220697,567567000))
    print('One ordinary A1 matrix element:',*single, flush=True)


if __name__ == '__main__':
    check()
