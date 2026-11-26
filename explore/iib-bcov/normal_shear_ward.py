"""Two- and three-loop Ward diagrams for the polynomial normal shear."""
from fractions import Fraction
from itertools import product

import networkx as nx

from bulk_eight import SHAPES, potentials
from bulk_three import evaluate as transverse_period
from planar_bulk import evaluate as planar_period
from primary_normal_jets import U, W, X, Q, Y, add, scaled
from square_three import skeleton_chains, tight_graphs


def projected_shape(shape):
    """The first vector is e^z4 d/dz1, the second is x² d/dy."""
    tangent_vector, shear, root = 0, 1, 7
    pieces = [Fraction(), Fraction()]
    for shear_sign, cubic_sign in product((-1, 1), repeat=2):
        fields = [Q, scaled(shear_sign, X), U, W, add(X, scaled(cubic_sign, Y)), Y, Y]
        total = tuple(map(sum, zip(*fields)))
        wave = (*fields, scaled(-1, total))
        phi = potentials(wave, shapes=(shape,))
        # With the tangent and shear odd parameters on the left in this order,
        # the cochain contraction has both terms and the same balancing slot.
        tangent_part = -total[4]*(phi[tangent_vector][0]-phi[root][0])
        shear_part = total[0]*(phi[shear][4]-phi[root][4])
        # The heat cofactor squared is doubly soft at each bivalent site.
        # Both y donors use two y jets, leaving one for the odd cubic.
        # Its even x jet therefore uses two x powers, as does the even shear.
        # Restore the two x² Taylor factorials and divide by 2! for the y donors.
        weight = Fraction(cubic_sign, 2)
        pieces[0] += weight*tangent_part
        pieces[1] += weight*shear_part
    return tuple(pieces)


def genus_three_topologies():
    # Each of the nine four-hard trivalent diagrams has rational and zeta(3)
    # columns. Both electrical numerators retain the original marked graph.
    by_core = {"double square": [], "tetrahedron": []}
    for graph in tight_graphs(3,3):
        network = nx.Graph(graph.edges)
        if max(degree for _, degree in network.degree) != 3:
            continue
        branches, paths = skeleton_chains(graph.edges, len(network))
        assert len(branches) == 4 and len(paths) == 6
        endpoints = {frozenset((nodes[0], nodes[-1])) for nodes, _ in paths}
        assert len(endpoints) in (4, 6)
        topology = "tetrahedron" if len(endpoints) == 6 else "double square"
        by_core[topology].append((transverse_period(graph.edges), planar_period(graph.edges)))

    assert [len(rows) for rows in by_core.values()] == [3, 6]
    all_rows = by_core["double square"]+by_core["tetrahedron"]
    assert tuple(sum(k[i] for k, _ in all_rows) for i in range(2)) == (
        -Fraction(4059837269,70614415872000), Fraction(1285913,27243216000))
    assert tuple(sum(j[i] for _, j in all_rows) for i in range(2)) == (
        -Fraction(2078212919,635529742848000), Fraction(131371,49037788800))
    result = {}
    for topology, rows in by_core.items():
        result[topology] = tuple(sum(2*k[i]/9+4*j[i] for k,j in rows) for i in range(2))
    assert result == {
        "double square": (-Fraction(347,6858432000), Fraction(0)),
        "tetrahedron": (-Fraction(370831,14370048000), Fraction(22219,1047816000)),
    }
    total = tuple(map(sum, zip(*result.values())))
    assert total == (-Fraction(7802719,301771008000), Fraction(22219,1047816000))
    print("Three-loop four-hard cores, rational and zeta(3) columns:", result, "sum:", total)


if __name__ == "__main__":
    totals = []
    expected = [
        (Fraction(-1,7290), Fraction(43,544320)),
        (Fraction(-17,51030), Fraction(29,136080)),
        (Fraction(-1,5670), Fraction(1,8640)),
        (Fraction(-1,5670), Fraction(1,7560)),
    ]
    for shape, target in zip(SHAPES, expected):
        parts = projected_shape(shape)
        assert parts == target, (shape, parts, target)
        print(shape, "tangent and shear voltage differences:", *parts, flush=True)
        totals.append(sum(parts))
    assert sum(totals) == -Fraction(11,38880)
    print("Shape totals:", *totals, "sum:", sum(totals), flush=True)
    genus_three_topologies()
