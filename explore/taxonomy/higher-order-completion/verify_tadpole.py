"""Replay all 55 original-edge tadpole moments and independent exact controls."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
import tempfile
import time

import tadpole_moments as calculation
from scalar_cache import fixed_parameter_cache

HERE = Path(__file__).resolve().parent
EXPECTED_EDGES = ((0, 1), (0, 4), (1, 2), (1, 5), (2, 3),
                  (2, 5), (3, 4), (3, 5), (4, 5))
EXPECTED_SOURCES = ((1,), (0,), (0,), (0,), (0,), (-1,))


def scalar(record):
    assert set(record) <= {"sqrtpi^4"}
    return F(record.get("sqrtpi^4", 0))


def independent_subset_cycles(edges):
    out = []
    for mask in range(1, (1 << len(edges)) - 1):
        parent = list(range(6))

        def find(v):
            while parent[v] != v:
                v = parent[v]
            return v

        loops = 0
        selected = []
        for i, (u, v) in enumerate(edges):
            if mask & (1 << i):
                selected.append(i)
                a, b = find(u), find(v)
                if a == b:
                    loops += 1
                else:
                    parent[b] = a
        out.append((selected, loops))
    return out


def positive_semidefinite_rank(matrix):
    matrix = [list(map(F, row)) for row in matrix]
    rank = 0
    while matrix:
        assert all(matrix[i][j] == matrix[j][i]
                   for i in range(len(matrix)) for j in range(len(matrix)))
        pivot = next((i for i, row in enumerate(matrix) if row[i]), None)
        if pivot is None:
            assert all(not x for row in matrix for x in row)
            break
        assert matrix[pivot][pivot] > 0
        rest = [i for i in range(len(matrix)) if i != pivot]
        matrix = [[matrix[i][j] - matrix[i][pivot] * matrix[pivot][j] / matrix[pivot][pivot]
                   for j in rest] for i in rest]
        rank += 1
    return rank


def main():
    started = time.monotonic()
    path = HERE / "tadpole_all_edges.json"
    saved = json.loads(path.read_text())
    assert saved["complete"] is True
    assert saved["edges"] == list(map(list, EXPECTED_EDGES))
    assert saved["sources"] == list(map(list, EXPECTED_SOURCES))
    assert saved["scope"] == "all original-edge words of degree at most two"
    assert saved["H"] == saved["T"] == 1 and saved["pole"] == "9/2"
    assert saved["implementation_sha256"] == calculation.dependencies()
    values = {tuple(row["edges"]): scalar(row["value"]) for row in saved["moments"]}
    expected = {()} | {(e,) for e in range(9)} | set(combinations_with_replacement(range(9), 2))
    assert set(values) == expected and len(saved["moments"]) == len(expected) == 55
    for row in saved["moments"]:
        assert F(row["pole"]) == F(9, 2) + len(row["edges"])
    bare = values[()]
    first = [values[(e,)] for e in range(9)]
    second = [[values[tuple(sorted((e, f)))] for f in range(9)] for e in range(9)]
    incidence = [[int(u == v) - int(w == v) for u, w in EXPECTED_EDGES] for v in range(6)]
    for v in range(6):
        assert sum(incidence[v][e] * first[e] for e in range(9)) == EXPECTED_SOURCES[v][0] * bare
        for f in range(9):
            assert sum(incidence[v][e] * second[e][f] for e in range(9)) == EXPECTED_SOURCES[v][0] * first[f]
    reflection = {0: 0, 1: 4, 2: 3, 3: 2, 4: 1, 5: 5}
    reflected = []
    for u, v in EXPECTED_EDGES:
        a, b = reflection[u], reflection[v]
        reflected.append((EXPECTED_EDGES.index(tuple(sorted((a, b)))), 1 if a < b else -1))
    for word, value in values.items():
        image = tuple(sorted(reflected[e][0] for e in word))
        assert value == prod(reflected[e][1] for e in word) * values[image]
    covariance = [[second[e][f] / bare - first[e] * first[f] / bare ** 2
                   for f in range(9)] for e in range(9)]
    assert positive_semidefinite_rank(covariance) == 4
    assert positive_semidefinite_rank(second) == 5
    # The earlier seven nilpotency equations are controls, not input values.
    roots = (3, 5, 7, 8)
    a, b, c, d = roots
    at = lambda i, j: second[i][j]
    assert first[a] + first[d] == F(4, 4725)
    assert first[b] + first[c] == F(1, 1350)
    assert at(a, a) + at(d, d) == F(47, 145530)
    assert at(a, b) + at(c, d) == F(41, 242550)
    assert at(b, b) + at(c, c) == F(37, 145530)
    assert at(b, d) + 2 * at(a, d) + at(a, c) == F(43, 121275)
    assert at(b, c) - at(a, d) == -F(1, 53900)
    print("all 55 moments: source conservation, reflection, positivity, seven old equations", flush=True)

    ctx, polynomial, kirchhoff, currents, partitions, controls = calculation.derive_density()
    subsets = independent_subset_cycles(EXPECTED_EDGES)
    powers = {degree: set() for degree in range(3)}
    for word in sorted(expected):
        numerator = polynomial * prod(currents[e] for e in word)
        powers[len(word)].update(tuple(int(x) + 1 for x in exponent)
                                 for exponent in numerator.to_dict())
    minimum_twice_margin = None
    for degree, entries in powers.items():
        for alpha in entries:
            assert sum(alpha) == 4 * (F(9, 2) + degree)
            for selected, loops in subsets:
                margin = 2 * sum(alpha[e] for e in selected) - (9 + 2 * degree) * loops
                assert margin > 0, (degree, alpha, selected)
                minimum_twice_margin = margin if minimum_twice_margin is None else min(minimum_twice_margin, margin)
    print("strict projective convergence for every original numerator monomial", flush=True)

    # Compare the guarded cache to fresh original calls, including q switches.
    scalar_controls = []
    requests = []
    for degree, entries in powers.items():
        ordered = sorted(entries)
        requests.extend((F(9, 2) + degree, ordered[i])
                        for i in (0, len(ordered) // 2))
    cold = []
    for q, alpha in requests:
        calculation.cores._systems.clear()
        cold.append(calculation.cores.period(EXPECTED_EDGES, alpha, q))
    with fixed_parameter_cache(maxsize=128) as cache_stats:
        for index in (0, 1, 2, 3, 4, 5, 4, 4, 0, 5):
            q, alpha = requests[index]
            calculation.cores._systems.clear()
            actual = calculation.cores.period(EXPECTED_EDGES, alpha, q)
            assert actual == cold[index]
        for (q, alpha), value in zip(requests, cold):
            scalar_controls.append({"q": str(q), "alpha": alpha,
                                    "complete_period_value": {calculation.cores.enc(k): str(v) for k, v in value.items()}})
    assert cache_stats["last_reduction_cache"]["hits"] > 0
    print("six complete scalar periods agree with the original cold route", flush=True)

    calculation.cores._systems.clear()
    with tempfile.TemporaryDirectory(prefix="tadpole-exact-replay-") as directory:
        output = Path(directory) / "all.json"
        calculation.compute(output, root_only=False)
        fresh = json.loads(output.read_text())
    for field in ("edges", "sources", "H", "T", "pole", "bare_numerator", "kirchhoff",
                  "current_numerators", "pfaffian_partitions", "original_graph_controls",
                  "implementation_sha256", "complete", "moments", "scope", "scalar_periods_evaluated"):
        assert fresh[field] == saved[field], field
    root_saved = json.loads((HERE / "tadpole_moments.json").read_text())
    assert root_saved["complete"] is True
    assert root_saved["implementation_sha256"] == calculation.dependencies()
    root_expected = {()} | {(e,) for e in roots} | set(combinations_with_replacement(roots, 2))
    assert {tuple(row["edges"]) for row in root_saved["moments"]} == root_expected
    assert all(scalar(row["value"]) == values[tuple(row["edges"])] for row in root_saved["moments"])
    report = {"schema": "tadpole50-independent-controls-v1", "result": "passed",
              "certificate_sha256": sha256(path.read_bytes()).hexdigest(),
              "root_certificate_sha256": sha256((HERE / "tadpole_moments.json").read_bytes()).hexdigest(),
              "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "complete_original_edge_moments_reintegrated": 55,
              "distinct_scalar_periods": saved["scalar_periods_evaluated"],
              "strictly_convergent_original_monomials_by_degree": {str(k): len(v) for k, v in powers.items()},
              "smallest_convergence_margin": str(F(minimum_twice_margin, 2)),
              "all_vertex_first_current_identities": 6,
              "all_vertex_second_current_identities": 54,
              "all_reflection_identities": 55,
              "positive_second_moment_rank": 5, "positive_covariance_rank": 4,
              "earlier_nilpotency_constraints": 7,
              "cold_scalar_controls": scalar_controls, "cache_control": cache_stats,
              "seconds": round(time.monotonic() - started, 2)}
    calculation.atomic_json(HERE / "tadpole_verification.json", report)
    print("PASS: full direct reintegration and independent controls", flush=True)


if __name__ == "__main__":
    main()
