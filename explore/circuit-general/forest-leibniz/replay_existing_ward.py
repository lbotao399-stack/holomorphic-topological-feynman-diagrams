"""Replay a complete three-term graph Ward identity from existing full-source data.

All seven source vectors are retained. The calculation uses exact rational
coefficient arithmetic, so equality of the displayed expressions is proved
without assuming arithmetic independence of their named zeta constants.
The data and the independent lower-face integration predate this replay.
"""

from fractions import Fraction
from itertools import combinations
import json
from math import prod
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "explore/nonSP/verify/full_ward"
EDGES = ((0, 1), (0, 5), (0, 6), (1, 2), (1, 5),
         (1, 6), (2, 3), (3, 4), (4, 5), (4, 6))
SOURCES = ((1, 2), (-2, 1), (3, -1), (-1, -2), (2, 3), (-4, 1), (1, -4))
EXPECTED = (Fraction(-66661, 2304), Fraction(0), Fraction(0))
EXPECTED_SIGNED_TERMS = {
    (1, 2): (Fraction(105245, 1152), Fraction(286719, 1120), Fraction(-5373, 14)),
    (2, 3): (Fraction(61350797, 80640), Fraction(2377057, 1120), Fraction(-22570, 7)),
    (3, 4): (Fraction(-11841847, 13440), Fraction(-83243, 35), Fraction(50513, 14)),
}


def permutation_sign(sequence):
    return (-1) ** sum(left > right for i, left in enumerate(sequence)
                      for right in sequence[i + 1:])


def evaluate_record(data, index, sources):
    assert tuple(map(sum, zip(*sources))) == (0, 0)
    brackets = tuple(sources[i][0] * sources[j][1] - sources[i][1] * sources[j][0]
                     for i, j in data["pairs"])
    value = [Fraction(0), Fraction(0), Fraction(0)]
    for column, coefficients in data["coefficients"][index].items():
        monomial = prod(brackets[i] for i in data["seeds"][int(column)])
        for channel, coefficient in enumerate(coefficients):
            value[channel] += monomial * Fraction(coefficient)
    return tuple(value)


def replay():
    data = json.loads((DATA / "polynomials.json").read_text())
    targets = [nx.Graph(edges) for edges in data["edges"]]
    rows = []
    for position, (left, right) in enumerate(EDGES):
        collapsed = {vertex: left if vertex == right else vertex for vertex in range(7)}
        remaining = EDGES[:position] + EDGES[position + 1:]
        mapped = [tuple(sorted((collapsed[i], collapsed[j]))) for i, j in remaining]
        if len(set(mapped)) != len(mapped):
            continue  # parallel-current wedge is zero
        quotient = nx.Graph(mapped)
        if any(quotient.subgraph(vertices).number_of_edges() > 2 * size - 3
               for size in range(3, 6) for vertices in combinations(quotient, size)):
            continue  # violates the proper-subgraph Laman bound
        for index, target in enumerate(targets):
            matcher = nx.isomorphism.GraphMatcher(quotient, target)
            if not matcher.is_isomorphic():
                continue
            image = {vertex: matcher.mapping[collapsed[vertex]] for vertex in range(7)}
            target_edges = tuple(map(tuple, data["edges"][index]))
            ordered = [target_edges.index(tuple(sorted((image[i], image[j]))))
                       for i, j in remaining]
            sign = (-1) ** position * permutation_sign(ordered)
            sources = [[0, 0] for _ in range(6)]
            for vertex, source in enumerate(SOURCES):
                for component, coefficient in enumerate(source):
                    sources[image[vertex]][component] += coefficient
            value = evaluate_record(data, index, sources)
            signed = tuple(sign * coefficient for coefficient in value)
            assert signed == EXPECTED_SIGNED_TERMS[left, right]
            rows.append({"contracted_edge": [left, right], "catalog_index": index,
                         "edge_order_sign": sign, "merged_sources": sources,
                         "value": list(map(str, value)),
                         "signed_value": list(map(str, signed))})
            break
        else:
            raise AssertionError("An admissible quotient was absent from the existing catalog")
    assert len(rows) == 3
    total = tuple(sum(Fraction(row["signed_value"][channel]) for row in rows)
                  for channel in range(3))
    assert total == EXPECTED

    lower = json.loads((DATA / "lower_faces.json").read_text())
    assert tuple(map(tuple, lower["sources"])) == SOURCES
    face = next(record for record in lower["graphs"]
                if tuple(map(tuple, record["edges"])) == EDGES)
    assert face["faces"] == [[[0, 1, 4, 5, 6], ["-66661/2304", "0"]]]
    assert tuple(map(Fraction, face["rhs"])) + (Fraction(0),) == total
    return {"edges": EDGES, "sources": SOURCES, "terms": rows,
            "signed_sum": list(map(str, total)), "lower_face": face["faces"][0],
            "existing_data": ["full_ward/polynomials.json", "full_ward/lower_faces.json"],
            "scope": "graph Ward identity; complete theory-specific Wick assembly is separate"}


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
