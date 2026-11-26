"""Complete the full-source six-loop K5 example with certified core CRT.

Run from any directory:
  python3 explore/circuit-general/series-fibres/run_core.py --threads 4

The earlier examples.json stores fixed-resistance fibre values. This script
performs the remaining projective integral and writes a separate certificate.
Only exact integer/Fraction arithmetic enters the mathematical result.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk import brute, catalog
from htk.fibres import integrate_fibre
from htk.graph import Graph
from htk.series import SeriesReduction
from core_numerator import certified_core_numerator, integrate_core


def certify_nonzero(values):
    """Prove a displayed 1,zeta(3),zeta(5) value is nonzero by rational bounds.

    No independence assumption about the zeta values is made. For s>1,
    sum_{k=1}^N k**(-s) + 1/((s-1)*(N+1)**(s-1)) < zeta(s)
    < sum_{k=1}^N k**(-s) + 1/((s-1)*N**(s-1)).
    """
    if set(values) - {"1", "z3", "z5"}:
        raise ValueError("nonzero certificate supports only 1,zeta(3),zeta(5)")
    for cutoff in (16, 32, 64, 128, 256, 512, 1024):
        lower = upper = Fraction(values.get("1", 0))
        for exponent in (3, 5):
            partial = sum((Fraction(1, k ** exponent) for k in range(1, cutoff + 1)), Fraction(0))
            left = partial + Fraction(1, (exponent - 1) * (cutoff + 1) ** (exponent - 1))
            right = partial + Fraction(1, (exponent - 1) * cutoff ** (exponent - 1))
            coefficient = Fraction(values.get("z%d" % exponent, 0))
            lower += coefficient * (left if coefficient >= 0 else right)
            upper += coefficient * (right if coefficient >= 0 else left)
        if lower > 0 or upper < 0:
            scale = 10 ** 9
            rounded_lower = Fraction((lower.numerator * scale) // lower.denominator, scale)
            rounded_upper = Fraction(-((-upper.numerator * scale) // upper.denominator), scale)
            if rounded_lower > 0 or rounded_upper < 0:
                lower, upper = rounded_lower, rounded_upper
            return dict(lower=str(lower), upper=str(upper), excludes_zero=True, cutoff=cutoff,
                        method="exact rational zeta partial sums and integral-test tail bounds",
                        zeta_linear_independence_assumed=False)
    raise ArithmeticError("rational zeta bounds did not certify a nonzero value")


def _exact_polynomial(polynomial, conductances, loops):
    exponents = loops - 1 + np.asarray(
        polynomial.certificate["chain_lengths"], dtype=np.int64) - polynomial.powers
    selected = np.fromiter((bool(integer) for integer in polynomial.integers), dtype=bool,
                           count=len(polynomial.integers))
    exponents = exponents[selected]
    integers = [value for value in polynomial.integers if value]
    weights = [1] * len(exponents)
    for axis, value in enumerate(conductances):
        powers = [int(value) ** k for k in range(loops)]
        weights = [old * powers[int(exponent)] for old, exponent in zip(weights, exponents[:, axis])]
    return Fraction(sum(integer * weight for integer, weight in zip(integers, weights)),
                    polynomial.denominator)


def replay_fibres(reduction, polynomial):
    """Rational cubature replays, independent of modular Gauss interpolation."""
    records = []
    for conductances in ([1] * len(reduction.chains),
                         [1 + i % 3 for i in range(len(reduction.chains))]):
        resistances = [Fraction(1, value) for value in conductances]
        result = integrate_fibre(reduction, resistances)
        centre = [tuple(Fraction(1, chain.length) for _ in chain.edges)
                  for chain in reduction.chains]
        U = brute.symanzik_U(reduction.G, reduction.times(resistances, centre))
        jacobian = math.prod(R ** (chain.length - 1)
                             for R, chain in zip(resistances, reduction.chains))
        expected = (result.value / jacobian * U ** (reduction.n + reduction.H)
                    * math.prod(value ** (reduction.n - 1) for value in conductances))
        actual = _exact_polynomial(polynomial, conductances, reduction.n)
        if actual != expected:
            raise ArithmeticError("certified core polynomial disagrees with rational fibre cubature")
        records.append(dict(conductances=conductances, rational_cubature_nodes=result.nodes,
                            integrated_fibre_value=str(result.value), reciprocal_numerator=str(actual),
                            exact_equality=True))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--checkpoint-dir", type=Path,
                        default=Path(tempfile.gettempdir()) / "htk-certified-core")
    parser.add_argument("--output", type=Path, default=HERE / "six_loop_k5.json")
    arguments = parser.parse_args()
    fixture = next(record for record in json.loads((HERE / "examples.json").read_text())["examples"]
                   if record["name"] == "six_loop_K5_H2")
    graph = Graph(fixture["vertices"], fixture["original_oriented_edges"])
    reduction = SeriesReduction(graph, fixture["full_balanced_sources"])
    if not catalog.tight_ok(graph.to_networkx(), reduction.H):
        raise ArithmeticError("the full-source benchmark is not an admissible graph")
    start = time.monotonic()
    polynomial = certified_core_numerator(reduction, threads=arguments.threads,
                                          checkpoint_dir=arguments.checkpoint_dir)
    polynomial.certificate["chain_lengths"] = [chain.length for chain in reduction.chains]
    fibre_replays = replay_fibres(reduction, polynomial)
    print("Two independent rational fibre replays passed.", flush=True)
    values, period_certificate = integrate_core(polynomial)
    if "file" in period_certificate:
        period_certificate["file"] = str(Path(period_certificate["file"]).relative_to(ROOT))
    record = dict(name="six_loop_K5_H2", status="complete bare T=0 projective kernel",
                  convention="code/htk/CONVENTIONS.md: K_Gamma(lambda; z=0)",
                  loops=reduction.n, holomorphic_dimension=reduction.H,
                  vertices=graph.nv, original_oriented_edges=graph.edges,
                  full_balanced_sources=[[int(value) if Fraction(value).denominator == 1 else str(value)
                                          for value in row] for row in reduction.full_sources],
                  all_original_sources_retained=True, graph_admissibility_exactly_checked=True,
                  core_edges=polynomial.edges, q=polynomial.q,
                  chain_edge_ids=[chain.edges for chain in reduction.chains],
                  original_determinant_size=graph.ne - 1,
                  contracted_determinant_size=reduction.determinant_size,
                  finite_field_cubature_nodes=27,
                  exact_value={label: str(value) for label, value in sorted(values.items())},
                  nonzero_certificate=certify_nonzero(values),
                  polynomial_certificate=polynomial.certificate,
                  independent_rational_fibre_replays=fibre_replays,
                  period_certificate=period_certificate,
                  remaining_projective_integrals=0,
                  elapsed_seconds=round(time.monotonic() - start, 3))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(record, indent=2) + "\n")
    print("Complete six-loop K5 kernel: " + json.dumps(record["exact_value"], sort_keys=True), flush=True)
    print("Certificate: " + str(arguments.output), flush=True)


if __name__ == "__main__":
    main()
