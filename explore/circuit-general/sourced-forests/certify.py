"""Recompute complete sourced-forest values and the explicit rational series.

Run with --write to create the certificate; the default checks every exact
field against the stored file. Every result retains the full original graph,
balanced sources, original current word and exact coefficient dictionary.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk import ProductForestMoments, ExactPeriod, ExactForestSeries
from fixtures import (sourced_fork, paired_higher_five, irreducible_quartic_five,
                      nonplanar_quartic_six, higher_paired_comb)


def coefficients(value):
    if isinstance(value, Fraction):
        return {"1": str(value)} if value else {}
    if not isinstance(value, ExactPeriod):
        raise ArithmeticError("a fixture requiring a complete finite result returned a series")
    return {label: str(scalar) for label, scalar in sorted(value.items())}


def record(name, graph, sources, marked, words):
    with ProductForestMoments(graph, sources, marked, check=1, allow_series=False,
                              max_terms=1000000) as engine:
        moments = []
        for word in words:
            value = engine.moment(word)
            moments.append({"current_powers": [[edge, component, power]
                                               for (edge, component), power in sorted(word.items())],
                            "complete_value": coefficients(value), "statistics": engine.last_stats})
    print(name, [item["complete_value"] for item in moments], flush=True)
    return {"name": name, "loops": graph.loops, "edges": graph.edges,
            "full_balanced_sources": sources, "marked": marked,
            "convention": "htk K_Gamma with the listed original edge order and orientations",
            "moments": moments}


def zeta_bounds(power, last):
    partial = sum((Fraction(1, k**power) for k in range(1, last + 1)), Fraction(0))
    return (partial + Fraction(1, (power - 1) * (last + 1)**(power - 1)),
            partial + Fraction(1, (power - 1) * last**(power - 1)))


def calculate():
    examples = []
    graph, sources, marked = sourced_fork()
    examples.append(record("sourced_fork", graph, sources, marked,
                           [{}, {(0, 0): 1}, {(0, 0): 1, (6, 1): 1}]))
    graph, sources, marked = paired_higher_five()
    examples.append(record("paired_higher_five", graph, sources, marked,
                           [{}, {(6, 0): 1}, {(6, 0): 1, (10, 1): 2}, {(6, 0): 2, (0, 1): 1}]))
    graph, sources, marked = irreducible_quartic_five()
    hard = record("irreducible_quartic_five", graph, sources, marked, [{}])
    examples.append(hard)
    value = {label: Fraction(coefficient) for label, coefficient in hard["moments"][0]["complete_value"].items()}
    lower = upper = value["1"]
    for label in ("z3", "z5"):
        lo, hi = zeta_bounds(int(label[1:]), 32)
        coefficient = value[label]
        lower += coefficient * (lo if coefficient > 0 else hi)
        upper += coefficient * (hi if coefficient > 0 else lo)
    assert -3 < lower < upper < -2
    hard["strict_nonvanishing_certificate"] = {
        "method": "positive-series integral-test bounds for zeta(3) and zeta(5), through 32 terms",
        "rational_lower": str(lower), "rational_upper": str(upper),
        "consequence": "-3 < complete kernel < -2"}
    graph, sources, marked = nonplanar_quartic_six()
    examples.append(record("nonplanar_quartic_six", graph, sources, marked, [{}]))
    for count in (3, 5, 8):
        graph, sources, marked, word = higher_paired_comb(count)
        examples.append(record("higher_paired_comb_%d" % count, graph, sources, marked, [word]))
    graph, sources, marked = sourced_fork()
    with ProductForestMoments(graph, sources, marked, check=1, max_terms=1000000) as engine:
        series = engine.series({(0, 0): 1})
        assert isinstance(series, ExactForestSeries)
        terms = list(series.terms(2))
        exact_series = {
            "meaning": "complete infinite Gegenbauer series, with the explicit recurrence in DERIVATION.md",
            "core_edges": series.edges,
            "components": [{"denominator_power": q,
                            "numerator": [[list(exponent), str(coefficient)]
                                          for exponent, coefficient in sorted(numerator.items())]}
                           for q, numerator in series.components],
            "constant": str(series.constant),
            "certified_terms": [{"order": row.order, "term": str(row.term),
                                 "partial_sum": str(row.partial_sum),
                                 "absolute_tail_bound": str(row.error_bound)} for row in terms],
            "finite_period_control": {"1": "-2772371/43200", "z3": "2597/50"},
        }
    files = (ROOT / "code/htk/forest_contractions.py", ROOT / "code/htk/forest_periods.py",
             ROOT / "code/htk/_forest_period_worker.py", ROOT / "code/htk/product_forests.py",
             HERE / "fixtures.py", HERE / "certify.py", HERE / "comb_reference.py",
             HERE / "independent_duality.py", HERE / "independent_six.py")
    return {"method": "source-preserving local contractions, exact current brackets, and residual periods",
            "scope": "H=2,T=0 product forests at exact supplied rational sources",
            "implementation_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in files},
            "complete_finite_examples": examples, "exact_series_example": exact_series}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    start = perf_counter()
    result = json.loads(json.dumps(calculate()))
    path = HERE / "certificate.json"
    if arguments.write:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Wrote complete exact sourced-forest certificate.", flush=True)
    else:
        assert result == json.loads(path.read_text(encoding="utf-8")), "certificate differs"
        print("Every exact certificate field replays.", flush=True)
    print("Elapsed %.3f seconds." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
