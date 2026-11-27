"""Solve all source-covariance and integrated Kirchhoff equations exactly.

The unknowns are every coefficient of the homogeneous bare kernel and
every coefficient of all nine first-current polynomials.  No source is
set to zero in the result, and no current or coefficient is discarded.
"""
from fractions import Fraction
from hashlib import sha256
import argparse
import json
from pathlib import Path
import random
import sys

from flint import fmpq_mpoly_ctx
import networkx as nx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "code"), str(ROOT / "explore/circuit-general/cores")]
from htk.graph import Graph
import cores


def homogeneous(variables, degree):
    if variables == 1:
        yield (degree,)
    else:
        for first in range(degree + 1):
            for tail in homogeneous(variables - 1, degree - first):
                yield (first,) + tail


def rational(x):
    return Fraction(int(x.p), int(x.q)) if hasattr(x, "p") else Fraction(x)


class Elimination:
    def __init__(self, columns):
        self.columns = columns
        self.pivots = {}
        self.rows = 0
        self.hash = sha256()

    def add(self, row):
        self.rows += 1
        row = {int(k): rational(v) for k, v in row.items() if v}
        self.hash.update(json.dumps([[k, str(v)] for k, v in sorted(row.items())],
                                    separators=(",", ":")).encode() + b"\n")
        while row:
            column = min(row)
            if column not in self.pivots:
                pivot = row[column]
                self.pivots[column] = {k: v/pivot for k, v in row.items()}
                return True
            pivot = row[column]
            for k, v in self.pivots[column].items():
                value = row.get(k, Fraction(0)) - pivot*v
                if value:
                    row[k] = value
                else:
                    row.pop(k, None)
        return False

    def basis(self):
        result = []
        for free in sorted(set(range(self.columns)) - set(self.pivots)):
            vector = {free: Fraction(1)}
            for pivot, row in sorted(self.pivots.items(), reverse=True):
                value = -sum((a*vector.get(k, 0) for k, a in row.items() if k != pivot), Fraction(0))
                if value:
                    vector[pivot] = value
            result.append(vector)
        return result


def add(row, column, coefficient):
    row[column] = row.get(column, 0) + coefficient
    if not row[column]:
        del row[column]


def automorphisms(edges):
    graph = nx.Graph(edges)
    indices = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    result = []
    for mapping in nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
        permutation = tuple(mapping[i] for i in range(6))
        edge = tuple(indices[tuple(sorted((mapping[a], mapping[b])))] for a, b in edges)
        orientations = tuple(1 if (mapping[a], mapping[b]) == edges[edge[e]] else -1
                             for e, (a, b) in enumerate(edges))
        character = (-1) ** sum(a > b for i, a in enumerate(edge) for b in edge[i+1:])
        result.append((permutation, edge, orientations, character))
    return sorted(result)


def generate_constraints(name):
    edges = cores.FAMILIES[name]
    monomials = {degree: tuple(homogeneous(5, degree)) for degree in (4, 5)}
    index = {degree: {a: i for i, a in enumerate(monomials[degree])} for degree in (4, 5)}
    sizes = {degree: len(monomials[degree]) for degree in (4, 5)}
    columns = sizes[4] + 9*sizes[5]
    ctx = fmpq_mpoly_ctx.get(tuple("lambda%d" % v for v in range(5)))
    source = list(ctx.gens())
    full = source + [-sum(source)]
    system = Elimination(columns)
    all_actions = automorphisms(edges)
    for permutation, edge_permutation, orientations, character in all_actions:
        if permutation == tuple(range(6)):
            continue
        inverse = [permutation.index(i) for i in range(6)]
        arguments = [full[inverse[i]] for i in range(5)]
        transformations = {}
        for degree in (4, 5):
            transformations[degree] = []
            for exponent in monomials[degree]:
                polynomial = ctx.constant(1)
                for x, power in zip(arguments, exponent):
                    polynomial *= x**power
                transformations[degree].append({index[degree][tuple(map(int, a))]: rational(c)
                                                for a, c in polynomial.to_dict().items()})
        rows = [{} for _ in range(sizes[4])]
        for j, polynomial in enumerate(transformations[4]):
            for i, c in polynomial.items():
                add(rows[i], j, c)
            add(rows[j], j, -character)
        for row in rows:
            system.add(row)
        for old, new in enumerate(edge_permutation):
            rows = [{} for _ in range(sizes[5])]
            for j, polynomial in enumerate(transformations[5]):
                for i, c in polynomial.items():
                    add(rows[i], sizes[4]+new*sizes[5]+j, c)
                add(rows[j], sizes[4]+old*sizes[5]+j, -character*orientations[old])
            for row in rows:
                system.add(row)
    incidence = Graph(6, edges).incidence()
    for vertex in range(5):
        rows = [{} for _ in range(sizes[5])]
        for i in range(sizes[5]):
            for edge in range(9):
                add(rows[i], sizes[4]+edge*sizes[5]+i, incidence[edge][vertex])
        for i, exponent in enumerate(monomials[4]):
            raised = list(exponent)
            raised[vertex] += 1
            add(rows[index[5][tuple(raised)]], i, -1)
        for row in rows:
            system.add(row)
    return system, monomials, all_actions


def evaluate(vector, monomials, source):
    values = {}
    for degree, powers in monomials.items():
        values[degree] = []
        for exponent in powers:
            value = Fraction(1)
            for x, p in zip(source, exponent):
                value *= x**p
            values[degree].append(value)
    size4, size5 = len(monomials[4]), len(monomials[5])
    result = [sum((v*values[4][j] for j, v in vector.items() if j < size4), Fraction(0))]
    for edge in range(9):
        lower = size4 + edge*size5
        result.append(sum((v*values[5][j-lower] for j, v in vector.items()
                           if lower <= j < lower+size5), Fraction(0)))
    return result


def build(name):
    system, monomials, actions = generate_constraints(name)
    basis = system.basis()
    assert len(basis) == {"prism": 49, "K33": 5}[name]
    # Exact evaluation rows are selected until every free coefficient is
    # determined.  Every output polynomial is retained at every source.
    sampling = Elimination(len(basis))
    selected = []
    rng = random.Random(310731)
    for trial in range(100):
        point = tuple(rng.randint(-2, 3) for _ in range(5))
        values = [evaluate(b, monomials, point) for b in basis]
        previous = len(sampling.pivots)
        for output in range(10):
            sampling.add({j: value[output] for j, value in enumerate(values) if value[output]})
        if len(sampling.pivots) > previous:
            selected.append({"sources": point, "rank_after": len(sampling.pivots),
                             "evaluations": [[str(v) for v in value] for value in values]})
        if len(sampling.pivots) == len(basis):
            break
    if len(sampling.pivots) != len(basis):
        raise ArithmeticError("source evaluation matrix is not unisolvent")
    data = {"schema": "complete-source-current-covariance-v1", "core": name,
            "edges": cores.FAMILIES[name], "source_vertices": list(range(5)),
            "root_source": "minus the sum of all five original sources",
            "unknown_coefficients": system.columns, "constraint_rows": system.rows,
            "constraint_rank": len(system.pivots), "solution_dimension": len(basis),
            "constraint_sha256": system.hash.hexdigest(),
            "monomials": {str(k): v for k, v in monomials.items()},
            "basis": [[[int(k), str(v)] for k, v in sorted(b.items())] for b in basis],
            "automorphisms": actions, "samples": selected,
            "sample_rank": len(sampling.pivots)}
    path = HERE / (name + "_source_covariance.json")
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(name, "complete covariance dimension", len(basis), "rank", len(system.pivots),
          "of", system.columns, "samples", len(selected), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cores", nargs="*", default=["prism", "K33"])
    for name in parser.parse_args().cores:
        build(name)
