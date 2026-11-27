"""Replay an actual unavailable finite period as an exact rational series.

Run with --write to create missing_period_certificate.json. The default
reconstructs the source-complete numerators, hashes every coefficient,
executes the actual order-zero series recurrence, verifies finite-only
rejection, and compares every deterministic field with the stored record.
Timings are printed only; no monomial dump or timing enters the certificate.
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

from htk import Graph, ProductForestMoments, ExactForestSeries
from htk.forest_periods import MissingForestPeriod
from htk.product_forests import UnsupportedProductForestReduction


EDGES = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
SOURCES = ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (0, 0), (-1, -1))
MARKED = (5, 6)
OPTIONS = {"check": 2, "max_terms": 1000000}


def numerator_digest(terms):
    """Digest every coefficient in a specified ordered edge-variable basis."""
    rows = sorted((tuple(map(int, exponent)), str(Fraction(str(coefficient))))
                  for exponent, coefficient in terms if coefficient)
    encoded = json.dumps([[list(exponent), coefficient] for exponent, coefficient in rows],
                         separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return {"nonzero_terms": len(rows), "sha256": hashlib.sha256(encoded).hexdigest()}


def calculate():
    graph = Graph(7, EDGES, "missing-K33-plus-edge-five-loop-product-forest")
    with ProductForestMoments(graph, SOURCES, MARKED, **OPTIONS) as engine:
        series = engine.moment()
        if not isinstance(series, ExactForestSeries):
            raise ArithmeticError("the missing-period control did not select the exact series")
        statistics = dict(engine.last_stats)
        if statistics["contraction"]["contraction_steps"] != 0:
            raise ArithmeticError("the irreducible missing-period fixture acquired a local contraction")
        if series.reason != "K33+e at q=7 not tabulated":
            raise ArithmeticError("the expected finite period is no longer the missing entry")
        # This diagnostic reads the cached exact numerator already constructed
        # by the public evaluation; no second determinant is computed.
        bare = numerator_digest(engine._core._bare.terms())
        components = [{"denominator_power": q, **numerator_digest(numerator.items())}
                      for q, numerator in series.components]
        row = next(series.terms(0))
        if (row.order != 0 or row.partial_sum != row.term or row.error_bound < 0
                or any(not isinstance(value, Fraction) for value in
                       (row.term, row.partial_sum, row.error_bound, row.lower, row.upper))):
            raise ArithmeticError("the series recurrence did not produce the exact certified row")
        if not row.lower <= 0 <= row.upper:
            raise ArithmeticError("the recorded order-zero nonvanishing status changed")
        exact_series = {
            "meaning": "complete infinite Gegenbauer series with its exact rational recurrence and tail bound",
            "core_edges": series.edges,
            "components": components,
            "constant": str(series.constant),
            "unavailable_finite_period_reason": series.reason,
            "certified_terms": [{"order": row.order, "term": str(row.term),
                                  "partial_sum": str(row.partial_sum),
                                  "absolute_tail_bound": str(row.error_bound),
                                  "rational_lower": str(row.lower), "rational_upper": str(row.upper)}],
            "interval_contains_zero": True,
            "integral_nonvanishing_proved": False,
        }

    with ProductForestMoments(graph, SOURCES, MARKED, allow_series=False, **OPTIONS) as engine:
        try:
            engine.moment()
        except UnsupportedProductForestReduction as error:
            if not isinstance(error.__cause__, MissingForestPeriod):
                raise ArithmeticError("finite-only rejection was not caused by a missing period") from error
            finite_only = {"rejected": True, "exception": type(error).__name__,
                           "message": str(error), "cause": str(error.__cause__)}
        else:
            raise ArithmeticError("finite-only evaluation incorrectly returned a result")

    implementation = (
        ROOT / "code/htk/forest_contractions.py",
        ROOT / "code/htk/forest_periods.py",
        ROOT / "code/htk/_forest_period_worker.py",
        ROOT / "code/htk/product_forests.py",
        ROOT / "explore/nonSP/verify/k4_periods/current_minor.py",
        ROOT / "explore/nonSP/verify/hepp_series/period.py",
        HERE / "certify_missing_period.py",
    )
    return {
        "schema_version": 1,
        "scope": "H=2,T=0, exact supplied balanced sources, bare current word",
        "convention": "htk K_Gamma with the listed original edge order and orientations",
        "original_graph": {"vertices": graph.nv, "loops": graph.loops, "edges": graph.edges},
        "full_balanced_sources": SOURCES,
        "marked": MARKED,
        "current_powers": [],
        "evaluation_options": OPTIONS,
        "original_bare_numerator": {"denominator_power": graph.loops + 2, **bare},
        "numerator_digest_encoding": (
            "UTF-8 compact JSON of [[exponent_tuple_as_list, reduced_Fraction_string], ...], "
            "sorted lexicographically by the full exponent tuple; every nonzero coefficient included. "
            "Bare variables follow original edges; series variables follow core_edges after exact "
            "chain integration and graph-orbit aggregation."),
        "exact_series": exact_series,
        "finite_only_control": finite_only,
        "statistics": statistics,
        "implementation_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in implementation},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    start = perf_counter()
    result = json.loads(json.dumps(calculate()))
    path = HERE / "missing_period_certificate.json"
    if arguments.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("Wrote exact missing-period series certificate.", flush=True)
    else:
        expected = json.loads(path.read_text(encoding="utf-8"))
        if result != expected:
            raise ArithmeticError("an exact missing-period certificate field differs")
        print("Every exact missing-period certificate field replays.", flush=True)
    print("Verified 10,577 original numerator terms, 7,368 core terms, the exact order-zero "
          "series term/tail, and finite-only rejection.", flush=True)
    print("Elapsed %.3f seconds." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
