"""Exact signed-automorphism zeros for complete four-source bulk3 polynomials.

Independent point controls use the original grounded vertex Laplacian and all
48 original current entries. No linear-minor or scalar-period evaluator is
imported by this checker.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import networkx as nx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "explore/circuit-general/fiveloop/enum_d3.json"
sys.path.insert(0, str(ROOT / "code"))
from htk import brute
from htk.graph import Graph
from htk.linalg import det


def parity(values):
    return -1 if sum(x > y for i, x in enumerate(values) for y in values[i + 1:]) % 2 else 1


def frame(row):
    values = [[Fraction(0)] * 3 for _ in range(row["nv"])]
    for i, vertex in enumerate(row["S0"][:3]):
        values[vertex][i] = Fraction(1)
    values[row["S0"][-1]] = [Fraction(-1)] * 3
    return values


def control(row, mapping, point):
    graph = Graph(row["nv"], row["edges"])
    sources = frame(row)
    transformed_sources = [None] * graph.nv
    for vertex in range(graph.nv):
        transformed_sources[mapping[vertex]] = sources[vertex]
    lookup = {frozenset(edge): i for i, edge in enumerate(graph.edges)}
    edge_image, edge_sign = [], []
    transformed_point = [None] * graph.ne
    for e, (left, right) in enumerate(graph.edges):
        image = lookup[frozenset((mapping[left], mapping[right]))]
        edge_image.append(image)
        edge_sign.append(1 if (mapping[left], mapping[right]) == graph.edges[image] else -1)
        transformed_point[image] = point[e]
    source_matrix = [transformed_sources[vertex] for vertex in row["S0"][:3]]
    source_determinant = det(source_matrix)
    original_density = brute.density(graph, point, sources).rational()
    permuted_density = brute.density(graph, transformed_point, sources).rational()
    changed_density = brute.density(graph, transformed_point, transformed_sources).rational()
    vertex_orientation = parity([mapping[v] for v in range(graph.nv)])
    edge_orientation = (-1) ** edge_sign.count(-1)
    if permuted_density != -original_density:
        raise ArithmeticError("The full original density did not reverse under the stated resistance substitution.")
    if changed_density != source_determinant ** 5 * permuted_density:
        raise ArithmeticError("The complete source-determinant law failed.")
    if changed_density != vertex_orientation * edge_orientation * original_density:
        raise ArithmeticError("The original vertex/edge orientation law failed.")
    original_currents = brute.currents(graph, point, sources)
    changed_currents = brute.currents(graph, transformed_point, transformed_sources)
    for edge, image in enumerate(edge_image):
        if changed_currents[image] != [edge_sign[edge] * value for value in original_currents[edge]]:
            raise ArithmeticError("An original current or source entry was lost by the automorphism.")
    return {"original_resistances": list(map(str, point)),
            "permuted_resistances": list(map(str, transformed_point)),
            "original_density": str(original_density), "permuted_density": str(permuted_density),
            "fully_transformed_density": str(changed_density),
            "source_matrix_determinant": str(source_determinant),
            "all_original_currents": [[str(value) for value in spinor] for spinor in original_currents],
            "all_transformed_currents": [[str(value) for value in spinor] for spinor in changed_currents]}


def build():
    rows = json.loads(DATA.read_text())["bulk3"]
    if len(rows) != 410:
        raise ArithmeticError("The complete bulk3 census changed.")
    witnesses, distribution, classified = [], Counter(), []
    for row in rows:
        graph = nx.Graph(row["edges"])
        if sorted(v for v in graph if graph.degree(v) == 2) != row["S0"] or len(row["S0"]) != 4:
            raise ArithmeticError("The complete four-source support changed.")
        witness = None
        for mapping in nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
            vertices = [mapping[v] for v in range(row["nv"])]
            source_permutation = [row["S0"].index(mapping[v]) for v in row["S0"]]
            reversed_edges = [e for e, (a, b) in enumerate(row["edges"]) if mapping[a] > mapping[b]]
            multiplier = parity(vertices) * (-1) ** len(reversed_edges) * parity(source_permutation) ** 5
            if multiplier == -1:
                points = [[Fraction(e + offset + 1, e + 1) for e in range(16)] for offset in (1, 3)]
                witness = {"index": row["index"], "core": row["core"],
                           "ordered_oriented_edges": row["edges"], "source_vertices": row["S0"],
                           "vertex_permutation": vertices, "source_permutation": source_permutation,
                           "vertex_sign": parity(vertices), "source_determinant_sign": parity(source_permutation),
                           "reversed_original_edges": reversed_edges, "bare_multiplier": -1,
                           "complete_bare_coefficient": {},
                           "independent_original_network_controls": [control(row, mapping, point) for point in points]}
                break
        distribution[row["core"], witness is not None] += 1
        classified.append({"index": row["index"], "core": row["core"],
                           "signed_automorphism_zero": witness is not None})
        if witness is not None:
            witnesses.append(witness)
    paths = (Path(__file__), DATA, ROOT / "code/htk/brute.py", ROOT / "code/htk/graph.py",
             ROOT / "code/htk/linalg.py", ROOT / "code/htk/CONVENTIONS.md")
    return {"schema": "five-loop-bulk3-signed-zeros-v1",
            "claim": "Complete bare source polynomials vanish; no claim that current moments vanish.",
            "all_classified_graphs": classified,
            "distribution": [{"core": core, "signed_zero": zero, "graphs": count}
                             for (core, zero), count in sorted(distribution.items())],
            "complete_bare_zeros": witnesses,
            "summary": {"classified_graphs": len(rows), "complete_bare_zeros": len(witnesses),
                        "zeros_outside_q8_available_families": sum(w["core"].startswith("V8_") and w["core"] != "V8_0"
                                                                  for w in witnesses),
                        "independent_density_transformation_points": 2 * len(witnesses),
                        "full_original_current_transformation_checks": 2 * 48 * len(witnesses)},
            "dependencies_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in paths}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    path = HERE / "symmetry_certificate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif json.loads(path.read_text()) != result:
        raise ArithmeticError("The complete signed-zero certificate changed.")
    print(json.dumps(result["summary"], indent=2))
