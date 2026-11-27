"""The complete open-color three-loop quartic product.

Inputs are beta1^i1*(partial_{+a} beta2^i2)*gamma5 and b_A*beta3^i3.
Every input color index and the output c^B*gamma4^i4*gamma5 is retained.
The four-flavor superpotential tensor and all four representations are
arbitrary. gamma5 is a free gauge singlet.

The Wick census and full geometric first-current moments are reused from
the fingerprinted, independently replayable nonrenormalizable-product
certificate. No geometric integral or special color contraction is inferred
from a selected entry. Here each matching is reduced to an ordered generator
word, before choosing a representation or superpotential tensor.

In the physical differential-operator formula, MM words mean products of
the manuscript generators t_A,t_B, while fR means i*C^C_AB*t_C. The common
factor is +g^2/[128*(8*pi^2)^3]. Thus the representation relation replaces
fR on one W leg by the AB word minus the BA word. All twelve unreduced
coefficients and all ten reduced coefficients are retained exactly.

The complete zeta terms are coefficients times
([t_A,t_B]-i*C^C_AB*t_C) on W legs 1 and 2. No gauge-invariance identity
for W, chosen color contraction, or independence of periods is used.
The result is for this complete first-current channel, not all products.

Default CLI: regenerate and compare the complete deterministic certificate.
--write: write that certificate after an intentional implementation change.
"""
from collections import defaultdict
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREVIOUS = HERE.parent / "nonrenormalizable-product"
sys.path.insert(0, str(PREVIOUS))
import nonabelian_wick as wick

base = wick.base
GEOMETRY = PREVIOUS / "nonabelian_wick_certificate.json"


def matching_action(record):
    """Trace all color lines from W to their original external endpoints.

    A labels the ghost paired to the original b_A; B labels the surviving
    ghost. The action matrix has the external matter index as its row.
    Hence following a W-to-input path gives the reverse of matrix-product
    order. fR uses antisymmetry in the two ghost indices, without Jacobi.
    The returned scalar includes the gauge vertex's -1/2, if present.
    """
    pairs = record["pairs"]
    partners = {left: right for left, right in pairs}
    partners.update({right: left for left, right in pairs})
    if record["cubic_types"] == ("M", "M"):
        external = partners[("b", 0, 4, 0)]
        output = next(field for field in record["surviving_fields"] if field[0] == "c")
        letters = {external[2]: "A", output[2]: "B"}
        if set(letters) != {1, 2}:
            raise ArithmeticError("both distinct matter-vertex ghosts must be retained")
        words, seen = [], []
        for flavor in range(1, 5):
            field, path = ("gamma", flavor, 0, 0), []
            while (field in partners and partners[field][0] == "beta"
                   and partners[field][2] in (1, 2)):
                vertex = partners[field][2]
                if vertex in seen:
                    raise ArithmeticError("a color line contains an unaccounted cycle")
                seen.append(vertex)
                path.append(letters[vertex])
                field = ("gamma", flavor, vertex, 0)
            if path:
                words.append((flavor, tuple(reversed(path))))
        if sorted(seen) != [1, 2]:
            raise ArithmeticError("a matter generator was omitted from the full tensor")
        return ("MM", tuple(words)), Fraction(1)

    if set(record["cubic_types"]) != {"M", "G"}:
        raise ArithmeticError("only the exhaustively admissible MM, MG, GM cases are supported")
    gauge = 1 if record["cubic_types"][0] == "G" else 2
    matter = 3 - gauge
    flavor = record["matter_flavors"][matter - 1]
    external = partners[("b", 0, 4, 0)]
    if external[0] != "c" or external[2] != gauge:
        raise ArithmeticError("the external b source was not routed through the gauge vertex")
    if partners[("b", 0, gauge, 0)] != ("c", 0, matter, 0):
        raise ArithmeticError("the gauge output was not routed to the matter generator")
    if external[3] not in (0, 1):
        raise ArithmeticError("the original ordered gauge-ghost slot was lost")
    return ("fR", flavor), Fraction(-1, 2) * (1 if external[3] == 0 else -1)


def canonical_matching(record):
    """Relabel only the two integration slots; keep every source and current."""
    candidates = []
    for mapping in ((0, 1, 2, 3, 4), (0, 2, 1, 3, 4)):
        oriented = tuple((mapping[a], mapping[b]) for a, b in record["ordered_edges"])
        edges = tuple(sorted(tuple(sorted(edge)) for edge in oriented))
        permutation = tuple(edges.index(tuple(sorted(edge))) for edge in oriented)
        locations = tuple(mapping[next(vertex for species, flavor, vertex, slot
                                       in record["surviving_fields"] if species == name)]
                          for name in ("gamma", "c", "gamma5"))
        original_edge = record["current_edge"]
        orientation = 1 if oriented[original_edge][0] < oriented[original_edge][1] else -1
        sign = record["wick_sign"] * base._parity(permutation) * orientation
        candidates.append(((edges, locations), permutation[original_edge], sign))
    return min(candidates)


def grouped_matchings():
    records, rejected = wick.uncolored_matchings()
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    for record in records:
        key, edge, sign = canonical_matching(record)
        word, coefficient = matching_action(record)
        result[key][edge][word] += sign * coefficient
    counts = {kinds: sum(record["cubic_types"] == kinds for record in records)
              for kinds in (("M", "M"), ("M", "G"), ("G", "M"), ("G", "G"))}
    if counts != {("M", "M"): 20, ("M", "G"): 4, ("G", "M"): 4, ("G", "G"): 0}:
        raise ArithmeticError("the complete original Wick census changed")
    return dict(sorted(result.items())), records, rejected


def read_geometry():
    certificate = json.loads(GEOMETRY.read_text(encoding="utf-8"))
    for name, expected in certificate["implementation_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ArithmeticError("the geometric certificate has a stale implementation fingerprint: " + name)
    nodes = {tuple(sample["source_spinors_gamma4_c_gamma5"][2])
             for sample in certificate["source_polynomial_certificate"]["samples"]}
    if nodes != {(a, b) for a in range(5) for b in range(5 - a)}:
        raise ArithmeticError("the complete degree-four source grid is required")
    return certificate


def evaluate_sample(sample, groups=None):
    """Assemble the entire formal tensor, including every period label."""
    if groups is None:
        groups = grouped_matchings()[0]
    moments = {(tuple(map(tuple, row["ordered_edges"])),
                tuple(row["survivor_vertices_gamma4_c_gamma5"])):
               row["full_first_current_moments"] for row in sample["graphs"]}
    result = defaultdict(lambda: [defaultdict(Fraction), defaultdict(Fraction)])
    for key, edges in groups.items():
        for edge, words in edges.items():
            for word, coefficient in words.items():
                result[word]
                for coordinate in range(2):
                    for label, value in moments[key][str(edge)][coordinate].items():
                        result[word][coordinate][label] += coefficient * Fraction(value)
    return {word: tuple({label: value for label, value in coordinate.items() if value}
                        for coordinate in spinor) for word, spinor in sorted(result.items())}


def covariant_tensor(geometry, groups=None):
    samples = {tuple(sample["source_spinors_gamma4_c_gamma5"][2]): evaluate_sample(sample, groups)
               for sample in geometry["source_polynomial_certificate"]["samples"]}
    words = sorted(next(iter(samples.values())))
    if any(sorted(values) != words for values in samples.values()):
        raise ArithmeticError("a source specialization removed a formal tensor word")
    labels = {label for values in samples.values() for spinor in values.values()
              for coordinate in spinor for label in coordinate} | {"z3"}
    if not labels <= {"1", "z3"}:
        raise ArithmeticError("the available three-loop geometric expression contains a new period label")
    result = {}
    for word in words:
        channels = {}
        for label in sorted(labels):
            polynomials = tuple(base.interpolate({node: values[word][coordinate].get(label, 0)
                                                   for node, values in samples.items()})
                                for coordinate in range(2))
            restored = base.covariant_restore(polynomials)
            if restored:
                channels[label] = restored
        result[word] = channels
    return result


def _add_terms(destination, terms, scalar=1):
    for vector, powers, coefficient in terms:
        destination[vector, powers] += scalar * coefficient


def apply_representation_relation(tensor):
    """Use i*C^C_AB*t_C = t_A*t_B-t_B*t_A, retaining the entire tensor."""
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    for word, channels in tensor.items():
        if word[0] == "MM":
            replacements = ((word, 1),)
        elif word[0] == "fR":
            flavor = word[1]
            replacements = ((("MM", ((flavor, ("A", "B")),)), 1),
                            (("MM", ((flavor, ("B", "A")),)), -1))
        else:
            raise ValueError("unknown formal tensor word")
        for replacement, scalar in replacements:
            for label, terms in channels.items():
                _add_terms(result[replacement][label], terms, scalar)
    return {word: {label: tuple((vector, powers, coefficient)
                               for (vector, powers), coefficient in sorted(terms.items()) if coefficient)
                   for label, terms in sorted(channels.items()) if any(terms.values())}
            for word, channels in sorted(result.items())}


def _serial_tensor(tensor):
    return [{"generator_word": word, "complete_covariant_coefficient": channels}
            for word, channels in sorted(tensor.items())]


def _serial_values(values):
    return [{"generator_word": word, "complete_spinor": spinor}
            for word, spinor in sorted(values.items())]


def certify():
    geometry = read_geometry()
    groups, records, rejected = grouped_matchings()
    tensor = covariant_tensor(geometry, groups)
    reduced = apply_representation_relation(tensor)
    if len(tensor) != 12 or len(reduced) != 10:
        raise ArithmeticError("the complete universal tensor count changed")
    if any(set(channels) - {"1"} for channels in reduced.values()):
        raise ArithmeticError("the representation identity did not cancel the complete nonrational expression")
    defects = []
    allowed_nonrational_words = set()
    for flavor in (1, 2):
        ab = ("MM", ((flavor, ("A", "B")),))
        ba = ("MM", ((flavor, ("B", "A")),))
        fr = ("fR", flavor)
        allowed_nonrational_words.update((ab, ba, fr))
        first = {(vector, powers): coefficient for vector, powers, coefficient in tensor[ab]["z3"]}
        for other in (ba, fr):
            if {(vector, powers): coefficient for vector, powers, coefficient in tensor[other]["z3"]} != {
                    key: -value for key, value in first.items()}:
                raise ArithmeticError("the exact representation-defect factorization failed")
        defects.append({"W_leg": flavor, "full_spinor_coefficient": tensor[ab]["z3"],
                        "exact_tensor_factor": "[t_A,t_B]-i*C^C_AB*t_C"})
    if any("z3" in channels and word not in allowed_nonrational_words for word, channels in tensor.items()):
        raise ArithmeticError("an additional nonrational tensor word was omitted")
    replays = []
    for sample in geometry["independent_full_source_replays"]:
        values = evaluate_sample(sample, groups)
        vectors = sample["source_spinors_gamma4_c_gamma5"]
        for word, channels in tensor.items():
            restored = [defaultdict(Fraction), defaultdict(Fraction)]
            for label, terms in channels.items():
                spinor = base.evaluate_covariant(terms, vectors)
                for coordinate in range(2):
                    if spinor[coordinate]:
                        restored[coordinate][label] = spinor[coordinate]
            if tuple(map(dict, restored)) != values[word]:
                raise ArithmeticError("a full-source, full-word replay disagrees with the geometric data")
        replays.append({"source_spinors_gamma4_c_gamma5": vectors,
                        "all_formal_tensor_coefficients": _serial_values(values)})
    implementations = (Path(__file__).resolve(), PREVIOUS / "quartic_wick.py",
                       PREVIOUS / "nonabelian_wick.py", GEOMETRY)
    return {
        "schema": 1,
        "scope": "complete first-current g^2*hbar^3 coefficient with arbitrary quartic W tensor and arbitrary representations; all input and output color indices open",
        "inputs": ["beta1^i1*(partial_{+a}beta2^i2)*gamma5", "b_A*beta3^i3"],
        "output": "c^B*gamma4^i4*gamma5",
        "free_color_indices_in_order": ["i1", "i2", "A", "i3", "B", "i4"],
        "spectator": "gamma5 is a free gauge singlet and all of its source is retained",
        "representation_words": {"MM": "ordered products of t_A,t_B on the stated W legs",
                                 "fR": "i*C^C_AB*t_C on the stated W leg",
                                 "matrix_orientation": "the external matter index is the row; the original W index is the column"},
        "normalization": {"W_placements_over_factorial": "1/2", "graph_sign": -1,
                          "head_derivative_factor": -2, "source_substitution": "lambda=-partial/2 on every surviving field",
                          "source_current_total_degree": 7,
                          "physical_multiplier_for_displayed_t_words": "+g^2/(128*(8*pi^2)^3)"},
        "implementation_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in implementations},
        "geometric_data_reused_exactly": str(GEOMETRY.relative_to(ROOT)),
        "matching_counts": {"MM": 20, "MG": 4, "GM": 4, "GG": 0},
        "complete_matching_ledger": [record | {"universal_action": matching_action(record),
                                               "canonical_source_current_sign": canonical_matching(record)}
                                     for record in records],
        "excluded_matchings_with_witnesses": rejected,
        "canonical_graphs": [{"ordered_edges": key[0], "source_vertices_gamma4_c_gamma5": key[1],
                              "all_routed_current_words": [{"edge": edge,
                                                           "tensor_weights": [{"generator_word": word, "coefficient": coefficient}
                                                                              for word, coefficient in sorted(words.items())]}
                                                          for edge, words in sorted(edges.items())]}
                             for key, edges in sorted(groups.items())],
        "source_polynomial_proof": {"total_degree": 7, "spectator_source_degree_bound": 4,
                                    "normalization_sources": [[1, 0], [0, 1]],
                                    "complete_grid": [(a, b) for a in range(5) for b in range(5 - a)],
                                    "bound": "rank at most three for the spectator-source block of the bare current determinant, plus one inserted current",
                                    "restoration": "det(M)^3*M, then polynomial continuation to dependent sources"},
        "covariant_bracket_order": geometry["covariant_bracket_order"],
        "covariant_vector_order": geometry["covariant_vector_order"],
        "complete_unreduced_tensor": _serial_tensor(tensor),
        "complete_representation_defects": defects,
        "complete_rational_tensor": _serial_tensor(reduced),
        "independent_full_source_replays": replays,
        "all_nonrational_terms_cancel_by_representation_relation": True,
        "W_gauge_invariance_used": False,
        "period_independence_assumed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    certificate = base._jsonable(certify())
    target = HERE / "quartic_tensor_certificate.json"
    if args.write:
        target.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print("Wrote the complete universal tensor certificate.")
    else:
        if not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != certificate:
            raise ArithmeticError("the complete regenerated universal tensor certificate differs")
        print("Verified all 28 matchings, 12 open tensor words, and the complete 10-word rational result.")


if __name__ == "__main__":
    main()
