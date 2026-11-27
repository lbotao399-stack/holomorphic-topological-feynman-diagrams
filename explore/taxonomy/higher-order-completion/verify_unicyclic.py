"""Independent rational rank, kernel, color and provenance checks.

The default performs bounded independent controls and reintegrates one
complete composition.  --deep additionally regenerates all 74,458 exact
kernel integrals and compares every complete rational composition value.
"""
import argparse
import ast
from hashlib import sha256
import itertools
import json
from pathlib import Path
import tempfile
import time

import numpy as np
from flint import fmpq, fmpq_mat

import unicyclic_weights as production
import exact_holomorphic_kernel as kernel
import exact_colour
import fmode
import kern_h
import audit_spectators

HERE = Path(__file__).resolve().parent
GRAPHS = ["G2.1", "G2.2", "G2.3u", "G2.4u", "G2.5u"]


def main(deep=False):
    started = time.monotonic()
    path = HERE / "unicyclic_weights.json"
    data = json.loads(path.read_text())
    assert data["complete"] is True
    assert data["implementation_sha256"] == production.dependencies()
    assert data["operator_types"] == list(production.OPERATOR_TYPES)
    assert data["functional_types"] == list(production.FUNCTIONAL_TYPES)
    assert data["operator_jets"] == [[0, 1], [1, 2], [2]]
    pairs = [("Q0", g) for g in GRAPHS] + [("G1.1", "G1.1")] + [(g, "Q0") for g in GRAPHS]
    assert set(data["compositions"]) == {"|".join(pair) for pair in pairs}
    values = {key: list(map(fmpq, row["complete_rational_values"])) for key, row in data["compositions"].items()}
    assert all(len(row) == 8 for row in values.values())
    rows = [[values["Q0|" + g][i] + values[g + "|Q0"][i] for g in GRAPHS] for i in range(8)]
    right = [-value / 4 for value in values["G1.1|G1.1"]]
    matrix = fmpq_mat(rows)
    augmented = fmpq_mat([row + [value] for row, value in zip(rows, right)])
    assert matrix.rank() == augmented.rank() == 5
    selected = next(indices for indices in itertools.combinations(range(8), 5)
                    if fmpq_mat([rows[i] for i in indices]).det())
    minor = fmpq_mat([rows[i] for i in selected])
    solution = list(minor.solve(fmpq_mat([[right[i]] for i in selected])).entries())
    assert solution == list(map(fmpq, ("1/2", "1/2", "1", "1/2", "1/2")))
    assert all(sum(a * b for a, b in zip(row, solution)) == value for row, value in zip(rows, right))
    trees = fmpq_mat([row[:2] for row in rows])
    trees_augmented = fmpq_mat([row[:2] + [value] for row, value in zip(rows, right)])
    assert trees.rank() == 2 and trees_augmented.rank() == 3
    print("rank 5; unique full weights; all eight residuals zero; trees alone inconsistent", flush=True)

    bundles = [production.inputs(flavour) for flavour in (1, 2, 3, 4)]
    assert [bundle[-1] for bundle in bundles] == data["input_data"]
    for name, graph in bundles[0][4].items():
        assert data["graphs"][name] == {
            "edges": list(map(list, graph.edges)), "bulk_vertices": graph.nb,
            "root_spokes": graph.spokes, "output_slots": list(map(list, graph.outv)),
            "automorphism_order": graph.aut, "kernel_orientation": graph.gk.sign_n}
    if "provenance_only_revision" in data:
        module = ast.parse(Path(production.__file__).read_text())
        module.body = [node for node in module.body if not (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                                                            and node.name == "dependencies")]
        assert sha256(ast.dump(module, include_attributes=False).encode()).hexdigest() == data["provenance_only_revision"]["computational_module_ast_sha256"]

    # Independent integration method: factorial moments of the full expanded
    # integrand versus native exact antiderivatives at every simplex face.
    original_integral = kernel.integrate
    integral_controls = []

    def compare_integral(polynomial, blocks, dimension):
        result = original_integral(polynomial, blocks, dimension)
        if result:
            independent = kernel.integrate_monomials(polynomial, blocks, dimension)
            assert result == independent
            integral_controls.append({"complete_integrand_terms": len(polynomial),
                                      "simplex_blocks": blocks, "rational_integral": str(result)})
        return result

    kernel.integrate = compare_integral
    controls = []
    for outer, inner in pairs:
        colors, information, normalization, histories = production.color_coefficients(bundles[0], outer, inner)
        assert histories == data["compositions"][outer + "|" + inner]["all_Wick_histories"]
        accepted = 0
        for key in sorted(colors, key=repr):
            for (_, momenta), coefficient in sorted(colors[key].items()):
                value = kernel.key_kernel(information[key], momenta)
                if not value:
                    continue
                modular_inputs = [np.array(row, dtype=np.int64)[:, None] % kern_h.P for row in momenta]
                modular = int(fmode.key_kernel(information[key], modular_inputs, 1, 3)[0])
                reduced = int(value.p) % kern_h.P * pow(int(value.q), -1, kern_h.P) % kern_h.P
                assert reduced == modular
                controls.append({"composition": outer + "|" + inner, "history_key": repr(key),
                                 "all_five_momentum_spinors": momenta, "rational_value": str(value),
                                 "independent_modular_kernel": modular, "prime": kern_h.P})
                accepted += 1
                break
            if accepted == 3:
                break
        assert accepted == 3, (outer, inner, accepted)
        print(outer, inner, "three nonzero complete rational kernel controls", flush=True)
    kernel.integrate = original_integral
    assert len(controls) == len(integral_controls) == 33

    # A complete independent finite-functional replay, including all color
    # histories, all four flavour counts and all eight outputs.
    with tempfile.TemporaryDirectory(prefix="unicyclic-exact-replay-") as directory:
        output = Path(directory) / "fresh.json"
        production.compute(output, only=None if deep else {"Q0|G2.3u"})
        fresh = json.loads(output.read_text())
    for key, row in fresh["compositions"].items():
        assert row == data["compositions"][key], key
    assert fresh["input_data"] == data["input_data"]
    assert fresh["graphs"] == data["graphs"]
    audit_spectators.main()
    report = {"schema": "six-dimensional-exact-weights-independent-checks-v1", "result": "passed",
              "certificate_sha256": sha256(path.read_bytes()).hexdigest(),
              "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "all_eight_linear_equations": [[str(x) for x in row] for row in rows],
              "all_eight_right_hand_sides": list(map(str, right)),
              "matrix_rank": 5, "augmented_rank": 5,
              "nonzero_minor_rows": selected, "nonzero_minor_determinant": str(minor.det()),
              "graph_order": GRAPHS, "unique_weights": list(map(str, solution)),
              "all_residuals": ["0"] * 8,
              "trees_only_rank": 2, "trees_only_augmented_rank": 3,
              "nonzero_complete_kernel_controls": controls,
              "independent_factorial_moment_controls": integral_controls,
              "fully_reintegrated_compositions": list(fresh["compositions"]),
              "deep_full_reintegration": deep,
              "scope": "unique weights from the complete stated finite-functional equations; no universal all-order nilpotency claim",
              "seconds": round(time.monotonic() - started, 2)}
    production.atomic_json(HERE / "unicyclic_verification.json", report)
    print("PASS: exact corrected weights, full rank and independent controls", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deep", action="store_true")
    main(parser.parse_args().deep)
