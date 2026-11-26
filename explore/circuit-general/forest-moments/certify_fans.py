"""Replay the complete eighth-order fan jet at 128 loops.

Every original path and spoke displacement is nonzero. Coefficients are
exact Fractions, including the source phase from the path edges. The file
stores every coefficient, not numerical samples or a coefficient hash.

Run: python3 explore/circuit-general/forest-moments/certify_fans.py
Use --write to regenerate the certificate after an intentional change.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "code"))

from htk import fan_taylor


def rational(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def calculate():
    loops, order = 128, 8
    displacements = tuple((Fraction(edge % 5 + 1, 7), Fraction(edge % 7 + 1, 11))
                          for edge in range(2 * loops + 1))
    result = fan_taylor(loops, (2, -1), (3, 4), displacements, order)
    assert result.peak_states <= result.state_bound == 165
    assert all(result.coefficients)
    return {
        "convention": "htk K_Gamma; coefficient of u^degree in K(lambda; u*z)",
        "graph": "htk.graph.fan(128), original edge order and orientations",
        "loops": loops,
        "order": order,
        "first_source": [2, -1],
        "last_source": [3, 4],
        "hub_source": [-5, -3],
        "other_vertex_sources": "zero",
        "displacement_formula": ["(edge_index % 5 + 1)/7", "(edge_index % 7 + 1)/11"],
        "coefficients": [rational(value) for value in result.coefficients],
        "coefficients_divided_by_bare_kernel": [rational(value / result.coefficients[0])
                                               for value in result.coefficients],
        "peak_states": result.peak_states,
        "proved_state_bound": result.state_bound,
        "exact_transfer_additions": result.transfers,
        "global_ratio_coordinates": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    start = perf_counter()
    result = calculate()
    path = HERE / "fan_certificate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        action = "wrote"
    else:
        assert result == json.loads(path.read_text(encoding="utf-8")), "fan certificate differs"
        action = "verified"
    print("%s 128-loop fan: all 9 exact coefficients through degree 8; "
          "%d/%d states; %d transfer additions; %.3f s" %
          (action, result["peak_states"], result["proved_state_bound"],
           result["exact_transfer_additions"], perf_counter() - start))


if __name__ == "__main__":
    main()
