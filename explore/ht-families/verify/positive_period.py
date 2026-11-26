"""Two terminal choices and the positive four-loop comparison period."""
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path
import sys

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from htk.simplex import dirichlet


edges = [(0,1), (0,2), (0,4), (0,5), (1,2), (1,3), (1,4), (3,4), (4,5)]
times = tuple(Fraction(j+2,j+3) for j in range(len(edges)))
weights = (2,1,1,1,1,1,1,1,1)
power = Fraction(5,2)


def series(*children):
    return "series", children


def parallel(*children):
    return "parallel", children


def evaluate(node):
    """Return loop count, edge weight, resistance, U and integrated Dirichlet factors."""
    if isinstance(node, int):
        return 0, weights[node], times[node], 1, Fraction(1), 0
    operation, children = node
    values = [evaluate(child) for child in children]
    loops = sum(value[0] for value in values)
    weight = sum(value[1] for value in values)
    resistance = sum(value[2] for value in values)
    kirchhoff = prod(value[3] for value in values)
    parameters = [value[1] - power*value[0] for value in values]
    if operation == "parallel":
        loops += len(values)-1
        resistance = 1/sum(1/value[2] for value in values)
        kirchhoff *= prod(value[2] for value in values)/resistance
        parameters = [power-value for value in parameters]
    assert min(parameters) > 0
    coeff, pi_power = dirichlet(tuple(value-1 for value in parameters))
    return (loops, weight, resistance, kirchhoff,
            coeff*prod(value[4] for value in values),
            pi_power+sum(value[5] for value in values))


def main():
    bases = []
    for deleted in combinations(range(len(edges)), 4):
        tree = nx.Graph()
        tree.add_nodes_from(range(6))
        tree.add_edges_from(edge for j,edge in enumerate(edges) if j not in deleted)
        if nx.is_tree(tree):
            bases.append(deleted)
    assert len(bases) == 54
    kirchhoff = sum(prod(times[e] for e in basis) for basis in bases)
    masks = [sum(1<<e for e in basis) for basis in bases]
    slacks = []
    for subset in range(1, (1<<len(edges))-1):
        rank = max((basis & subset).bit_count() for basis in masks)
        weight = sum(value for e,value in enumerate(weights) if subset & (1<<e))
        slacks.append(power*rank - weight)
    assert len(slacks) == 510 and min(slacks) == Fraction(1,2)
    side04 = parallel(2,series(3,8))
    side14 = parallel(6,series(5,7))
    plans = [parallel(0,series(1,4),series(side04,side14)),
             parallel(1,series(4,parallel(0,series(side04,side14))))]
    for terminals,plan in zip(((0,1),(0,2)),plans):
        loops, weight, _, actual_u, coeff, exponent = evaluate(plan)
        assert (loops,weight,actual_u) == (4,10,kirchhoff)
        assert (coeff,exponent) == (Fraction(1,6),4)
        print(f"terminals {terminals}: U matches all 54 bases; period = {coeff} pi^{exponent}")
    print(f"{len(slacks)} proper edge subsets: minimum convergence exponent {min(slacks)}")


if __name__ == "__main__":
    main()
