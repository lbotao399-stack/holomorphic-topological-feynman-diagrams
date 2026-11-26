"""Exact five- and six-loop chain-fibre certificates, with full-current replay.

Run from any working directory:
    python3 explore/circuit-general/series-fibres/run_examples.py

The output is a value of the density AFTER all bivalent-chain integrations
at explicitly retained core resistances.  The remaining projective core
integral is not performed.  Every cubature ordinate is independently replayed
using the original vertex-Laplacian and full current determinant.
"""
import argparse
import json
import math
from pathlib import Path
import random
import sys
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from htk import brute, catalog
from htk.dirichlet import product_simplex_rule
from htk.fibres import integrate_fibre
from htk.graph import Graph
from htk.series import SeriesReduction


def subdivide(nv, edges, lengths):
    result = []
    for (u, v), length in zip(edges, lengths):
        path = [u] + list(range(nv, nv + length - 1)) + [v]
        nv += length - 1
        result.extend(zip(path, path[1:]))
    return Graph(nv, result)


def sources(G, H, seed):
    rng = random.Random(seed)
    while True:
        result = [[rng.randint(-3, 3) for _ in range(H)] for _ in range(G.nv - 1)]
        if any(not any(row) for row in result):
            continue
        result.append([-sum(row[i] for row in result) for i in range(H)])
        if any(result[-1]):
            return result


def fraction_json(value):
    value = F(value)
    return {"numerator": value.numerator, "denominator": value.denominator}


def examples():
    cube = [(0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    yield ("five_loop_cube_H3", 3,
           subdivide(8, cube, [2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1]), 10503)
    complete = [(u, v) for u in range(5) for v in range(u + 1, 5)]
    yield ("six_loop_K5_H2", 2,
           subdivide(5, complete, [2 if i in (0, 7, 9) else 1 for i in range(10)]), 10602)


def certify(name, H, G, seed):
    if not catalog.tight_ok(G.to_networkx(), H):
        raise AssertionError("higher-loop benchmark violates graph admissibility")
    lam = sources(G, H, seed)
    reduction = SeriesReduction(G, lam)
    resistances = tuple(F(1 + i % 3) for i in range(len(reduction.chains)))
    result = integrate_fibre(reduction, resistances)
    if not result.value:
        raise AssertionError("the generic higher-loop benchmark must be nonzero")
    rule = product_simplex_rule(result.simplex_dimensions, result.degree)
    jacobian = math.prod(R ** (chain.length - 1)
                         for chain, R in zip(reduction.chains, resistances))
    independent = F(0)
    checks = 0
    for point, weight in rule.nodes:
        times = reduction.times(resistances, point)
        reference = brute.density(G, times, lam).rational() * jacobian
        contracted = reduction.fibre_density(resistances, point)
        if reference != contracted:
            raise AssertionError("full original determinant differs at cubature node")
        independent += weight * reference
        checks += 1
    if independent != result.value:
        raise AssertionError("independently replayed fibre integral differs")
    midpoint = tuple(tuple(F(1, chain.length) for _ in chain.edges) for chain in reduction.chains)
    times = reduction.times(resistances, midpoint)
    U = brute.symanzik_U(G, times)
    return {
        "name": name,
        "status": "exact rational chain-fibre integral; remaining core period unevaluated",
        "loops": G.loops,
        "holomorphic_dimension": H,
        "vertices": G.nv,
        "original_oriented_edges": G.edges,
        "full_balanced_sources": lam,
        "source_seed": seed,
        "chain_edge_ids": [chain.edges for chain in reduction.chains],
        "chain_vertices": [chain.vertices for chain in reduction.chains],
        "chain_orientation_signs": [chain.signs for chain in reduction.chains],
        "retained_core_resistances": [fraction_json(R) for R in resistances],
        "projective_gauge_chain": reduction.gauge,
        "symanzik_U": fraction_json(U),
        "original_determinant_size": G.ne - 1,
        "contracted_determinant_size": result.determinant_size,
        "fixed_source_columns": len(reduction.fixed_columns),
        "fixed_source_rank": reduction.fixed_rank,
        "simplex_dimensions": result.simplex_dimensions,
        "polynomial_total_degree_bound": result.degree,
        "cubature_nodes": result.nodes,
        "independent_full_determinant_checks": checks,
        "integrated_fibre_value": fraction_json(result.value),
        "remaining_projective_core_dimension": len(reduction.chains) - 1,
        "remaining_integral": "Integrate the fibre-density function of R over positive projective core resistances; the recorded value is only its value at the displayed core point.",
        "all_original_sources_retained": True,
        "graph_admissibility_exactly_checked": True,
        "all_exact_checks_passed": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "examples.json")
    parser.add_argument("--log", type=Path, default=HERE / "examples.log")
    args = parser.parse_args()
    records, lines = [], []
    for example in examples():
        record = certify(*example)
        records.append(record)
        value = record["integrated_fibre_value"]
        line = ("%s: determinant %d -> %d; %d exact cubature nodes; %d independent full-determinant checks; fibre value %d/%d; remaining core dimension %d"
                % (record["name"], record["original_determinant_size"],
                   record["contracted_determinant_size"], record["cubature_nodes"],
                   record["independent_full_determinant_checks"], value["numerator"],
                   value["denominator"], record["remaining_projective_core_dimension"]))
        print(line, flush=True)
        lines.append(line)
    output = {"scope": "T=0; exact integration of all series-chain fractions with every original source; remaining projective core integrals not evaluated", "examples": records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text("\n".join(lines) + "\nAll exact checks passed. These are fibre densities, not complete higher-loop periods.\n")


if __name__ == "__main__":
    main()
