"""Complete original-current Ward checks through degree two, with higher words.

The main control contains all 153 monomials of degree at most two in the
sixteen original current coordinates. Nothing is discarded using source
equations. Every signed quotient period is compared with an independent
rational spider-boundary integration. The default replays every exact
certificate field; --write regenerates the certificate.
"""
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
from time import perf_counter

from ward_three import EDGES, GENERIC_SOURCES, PHYSICAL_SOURCES, WardThree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def all_words(last):
    locations = [(edge, component) for edge in range(len(EDGES)) for component in range(2)]
    for degree in range(last + 1):
        for chosen in combinations_with_replacement(locations, degree):
            powers = {}
            for location in chosen:
                powers[location] = powers.get(location, 0) + 1
            yield powers


def serialize(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    return value


def calculate_case(name, sources, words):
    records = []
    period_words = 0
    with WardThree(sources, check=1, max_terms=1000000) as engine:
        graphs = [{key: quotient[key] for key in
                   ("position", "contracted_edge", "sign", "edges", "sources")}
                  for quotient in engine.quotients]
        for word in words:
            result = engine.verify(word)
            coefficients = [row["coefficients"] for row in result["terms"]]
            labels = sorted({label for row in coefficients for label in row if label != "1"})
            period_words += bool(labels)
            cancellations = {label: sum(row["sign"] * row["coefficients"].get(label, 0)
                                        for row in result["terms"]) for label in labels}
            assert not any(cancellations.values())
            records.append({
                "powers": result["powers"],
                "quotient_values_in_edge_order": coefficients,
                "signed_sum": result["signed_sum"],
                "independent_rational_boundary": result["rational_boundary"],
                "exact_period_coefficient_cancellations": cancellations,
            })
            if len(records) % 32 == 0:
                print(name, "complete words", len(records), flush=True)
    print(name, "verified", len(records), "words;", period_words,
          "have nonzero period coefficients before the signed sum", flush=True)
    return dict(name=name, sources=sources, original_edges=EDGES, quotients=graphs,
                word_count=len(records), words_with_individual_periods=period_words,
                complete_words=records)


def calculate():
    complete = list(all_words(2))
    assert len(complete) == 153
    generic = calculate_case("generic_complete_degree_two", GENERIC_SOURCES, complete)
    higher_words = [
        {(0, 0): 2, (2, 1): 1},
        {(1, 1): 1, (3, 0): 1, (6, 1): 1},
        {(2, 0): 1, (4, 1): 1, (7, 0): 1},
        {(7, 1): 3},
    ]
    higher = calculate_case("generic_mixed_degree_three", GENERIC_SOURCES, higher_words)
    physical_words = [{}] + [{(edge, component): 1}
                            for edge in range(len(EDGES)) for component in range(2)]
    physical = calculate_case("higher_superpotential_source_pattern", PHYSICAL_SOURCES, physical_words)
    scaled_sources = [[Fraction(2, 3) * entry for entry in row] for row in GENERIC_SOURCES]
    scaled_words = [{}, {(0, 0): 1}, {(2, 0): 1}, {(0, 0): 1, (2, 1): 1}]
    scaled = calculate_case("nonintegral_conserved_sources", scaled_sources, scaled_words)
    # Homogeneity is checked against the complete generic coefficients,
    # retaining the same original current word and every quotient.
    reference = {tuple(tuple(item) for item in row["powers"]): row
                 for row in generic["complete_words"]}
    for row in scaled["complete_words"]:
        degree = sum(power for _, _, power in row["powers"])
        factor = Fraction(2, 3) ** (6 + degree)
        old = reference[tuple(tuple(item) for item in row["powers"])]
        assert row["independent_rational_boundary"] == factor * old["independent_rational_boundary"]
        for current, previous in zip(row["quotient_values_in_edge_order"],
                                     old["quotient_values_in_edge_order"]):
            assert current == {label: factor * value for label, value in previous.items()}
    files = [HERE / "ward_three.py", HERE / "certify_words.py",
             ROOT / "code/htk/forest_periods.py", ROOT / "code/htk/_forest_period_worker.py"]
    return serialize(dict(
        meaning="Complete current-word graph Ward identities, before theory-specific Wick assembly",
        period_handling="Explicit rational coefficients cancel; no arithmetic independence is assumed",
        scope="All original-current monomials through degree two at the specified generic sources; additional exact controls",
        implementation_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in files},
        cases=[generic, higher, physical, scaled],
        total_verified_words=sum(case["word_count"] for case in (generic, higher, physical, scaled)),
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    start = perf_counter()
    result = calculate()
    destination = HERE / "current_words_certificate.json"
    if arguments.write:
        destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Wrote the complete current-word Ward certificate.", flush=True)
    else:
        assert result == json.loads(destination.read_text(encoding="utf-8")), "certificate differs"
        print("Every exact current-word Ward certificate field replays.", flush=True)
    print("Elapsed %.3f seconds." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
