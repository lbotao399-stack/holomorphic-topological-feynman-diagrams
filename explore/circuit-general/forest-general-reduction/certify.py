"""Reproduce complete sourced-ear kernels and full first-current spinors.

Run: python3 explore/circuit-general/forest-general-reduction/certify.py
Use --write to regenerate the certificate. Original four/five-loop graphs
are also integrated independently for the indicated controls. No finite
period table is generated or modified.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from numbers import Integral
from pathlib import Path
import sys

from fixtures import sourced_ear_forest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from htk.forest_ears import ForestEarMoments
from htk.forest_contractions import ForestContractions
from htk.forest_periods import CoreCurrentPeriods


CASES = ((False, 1), (False, 2), (True, 1), (True, 2),
         (True, 3), (True, 5), (True, 8))


def coefficients(value):
    if isinstance(value, (Integral, F)):
        return {"1": F(value)} if value else {}
    return {key: F(scalar) for key, scalar in value.items() if scalar}


def word_polynomial(graph, word):
    exponent = [0] * (2 * graph.ne)
    for (edge, component), power in word.items():
        exponent[2 * edge + component] = power
    return {tuple(exponent): F(1)}


def source_basis_spinor(first, second, sources):
    """Retain both components as coefficients of original sources lambda0,lambda1."""
    left, right = sources[:2]
    determinant = left[0] * right[1] - left[1] * right[0]
    if not determinant:
        raise ArithmeticError("the documented full-spinor source basis became dependent")
    result = {}
    for label in sorted(set(first) | set(second)):
        v0, v1 = first.get(label, F(0)), second.get(label, F(0))
        result[label] = {"lambda_vertex_0": (v0 * right[1] - v1 * right[0]) / determinant,
                         "lambda_vertex_1": (left[0] * v1 - left[1] * v0) / determinant}
    return result


def case(physical, ears):
    graph, sources, marked, interpretation = sourced_ear_forest(ears, physical=physical)
    words = ({}, {(graph.ne - 1, 0): 1}, {(graph.ne - 1, 1): 1})
    checks = []
    direct = CoreCurrentPeriods(graph, sources, marked, check=1) if ears <= 2 else None
    try:
        with ForestEarMoments(graph, sources, marked, check=1) as engine:
            steps = engine.steps
            for word in words:
                actual = coefficients(engine.moment(word, allow_series=False))
                record = {"powers": tuple((edge, component, power) for (edge, component), power in sorted(word.items())),
                          "coefficients": actual, "ear_stats": dict(engine.last_stats)}
                # Both generic first-current components receive complete
                # original-graph controls. Physical specializations retain
                # an independent original-graph bare control as well.
                if direct is not None and (not physical or not word):
                    expected = coefficients(direct.integrate(word_polynomial(graph, word), allow_series=False))
                    if expected != actual:
                        raise ArithmeticError("sourced-ear and original-graph integrals disagree")
                    record["independent_original_graph"] = {"coefficients": expected,
                                                             "stats": dict(direct.last_stats)}
                checks.append(record)
                print("%s %d loops %s: %s" % ("physical" if physical else "generic", graph.loops,
                                              record["powers"], actual), flush=True)
    finally:
        if direct is not None:
            direct.close()
    passive = ForestContractions(graph, sources, marked, check=0)
    return {"physical_source_pattern": physical, "ear_count": ears, "loops": graph.loops,
            "edges": tuple(graph.edges), "sources": sources, "marked": marked,
            "field_count_interpretation": interpretation,
            "previous_zero_source_contraction_steps": len(passive.steps),
            "ear_steps": steps, "checks": tuple(checks),
            "full_first_current_spinor_in_original_source_basis":
                source_basis_spinor(checks[1]["coefficients"], checks[2]["coefficients"], sources)}


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def certificate():
    paths = (Path(__file__), HERE / "fixtures.py", ROOT / "code/htk/forest_ears.py",
             ROOT / "code/htk/forest_periods.py", ROOT / "code/htk/_forest_period_worker.py")
    return serial({"scope": "Complete individual forest kernels and first-current spinors; an assembled Wick/color product is not asserted.",
                   "max_terms": 1000000, "max_source_samples": 10000,
                   "records": tuple(case(physical, ears) for physical, ears in CASES),
                   "source_hashes": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    actual = certificate()
    path = HERE / "ear_certificate.json"
    if arguments.write:
        path.write_text(json.dumps(actual, sort_keys=True, indent=2) + "\n")
        print("Wrote the complete sourced-ear certificate.", flush=True)
    else:
        if actual != json.loads(path.read_text()):
            raise ArithmeticError("the complete sourced-ear certificate changed")
        print("Every complete sourced-ear certificate field replays.", flush=True)


if __name__ == "__main__":
    main()
