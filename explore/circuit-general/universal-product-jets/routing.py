"""Exact original-input current routing for every commutator matching.

Each current is represented by its six integer coefficients on
lambda_gamma4, lambda_c, lambda_gamma5, y_WA, y_WB, y_AB. This identity
holds for the full lower spinor, not only a chosen entry. All four
contracted input legs and both ordered gauge-ghost matchings are retained.
"""
from fractions import Fraction
import argparse
import hashlib
import json

import jet_tensor as jets
from htk.linalg import inverse

INPUTS = (("beta", 1, 3, 0), ("beta", 2, 3, 0), ("beta", 3, 4, 0), ("b", 0, 4, 0))
INPUT_NAMES = ("beta1", "beta2", "beta3", "b")
BASIS = ("lambda_gamma4", "lambda_c", "lambda_gamma5", "y_WA", "y_WB", "y_AB")


def exact_currents(record):
    """Solve the original tree incidence equations over the entire source basis."""
    graph = jets.universal.wick.Graph(5, record["ordered_edges"])
    chords = tuple(graph.edges.index(edge) for edge in ((0, 3), (0, 4), (3, 4)))
    tree = tuple(edge for edge in range(7) if edge not in chords)
    graph._check_tree(tree)
    incidence = graph.incidence()
    rows = [[Fraction(0)] * 6 for vertex in range(4)]
    rows[0][0] = 1
    ghost_vertex = next(field[2] for field in record["surviving_fields"] if field[0] == "c")
    rows[ghost_vertex][1] = 1
    rows[3][2] = 1
    currents = [None] * 7
    for position, edge in enumerate(chords):
        current = [Fraction(0)] * 6
        current[3 + position] = 1
        currents[edge] = tuple(current)
        for vertex in range(4):
            rows[vertex][3 + position] -= incidence[edge][vertex]
    inv = inverse([[incidence[edge][vertex] for edge in tree] for vertex in range(4)])
    for position, edge in enumerate(tree):
        currents[edge] = tuple(sum(inv[position][vertex] * rows[vertex][column] for vertex in range(4))
                               for column in range(6))
    for vertex in range(4):
        expected = [Fraction(0)] * 6
        if vertex == 0:
            expected[0] = 1
        if vertex == ghost_vertex:
            expected[1] = 1
        if vertex == 3:
            expected[2] = 1
        if tuple(sum(incidence[edge][vertex] * currents[edge][column] for edge in range(7))
                 for column in range(6)) != tuple(expected):
            raise ArithmeticError("an original full-source Kirchhoff equation failed")
    return tuple(currents), tree, chords


def certify():
    records = jets.universal.grouped_matchings()[1]
    direct = (Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0))
    other_A = (Fraction(0), Fraction(0), Fraction(-1), Fraction(-1), Fraction(0), Fraction(1))
    beta3 = (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(0))
    b = (Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(-1), Fraction(-1))
    ledger = []
    for record in records:
        word, scalar = jets.universal.matching_action(record)
        if word[0] == "fR":
            slot = word[1]
        elif len(word[1]) == 1 and len(word[1][0][1]) == 2:
            slot = word[1][0][0]
        else:
            continue
        if slot not in (1, 2):
            raise ArithmeticError("an additional commutator slot appeared")
        currents, tree, chords = exact_currents(record)
        maps = []
        for name, field in zip(INPUT_NAMES, INPUTS):
            pair = next(pair for pair in record["pairs"] if field in pair)
            edge = tuple(sorted((pair[0][2], pair[1][2])))
            if edge[1] != field[2]:
                raise ArithmeticError("an original input derivative was not at its stated oriented head")
            position = record["ordered_edges"].index(edge)
            maps.append({"input": name, "original_field": field, "original_edge": edge,
                         "edge_position": position, "full_spinor_affine_coefficients": currents[position],
                         "head_derivative_factor": -2})
        expected = (other_A, direct, beta3, b) if slot == 1 else (direct, other_A, beta3, b)
        if tuple(row["full_spinor_affine_coefficients"] for row in maps) != expected:
            raise ArithmeticError("a complete original input current differs across the commutator family")
        ledger.append({"cubic_types": record["cubic_types"], "matter_flavors": record["matter_flavors"],
                       "ordered_edges": record["ordered_edges"], "pairs": record["pairs"],
                       "wick_sign": record["wick_sign"], "generator_word": word, "color_prefactor": scalar,
                       "commutator_slot": slot, "full_source_vertices": {"lambda_gamma4": 0,
                           "lambda_c": next(field[2] for field in record["surviving_fields"] if field[0] == "c"),
                           "lambda_gamma5": 3, "dependent_balanced_source": 4},
                       "original_tree_edges": tree, "common_chord_edges": chords,
                       "all_original_input_current_maps": maps})
    if len(ledger) != 16:
        raise ArithmeticError("the complete commutator matching census changed")
    paths = (jets.HERE / "routing.py", jets.universal.HERE / "quartic_tensor.py",
             jets.universal.PREVIOUS / "nonabelian_wick.py")
    return {"schema": 1, "scope": "all simultaneous finite input-current words on beta1,beta2,beta3,b; every original source and free spinor index retained",
            "affine_spinor_basis_in_order": BASIS,
            "complete_original_matching_count": len(ledger), "all_original_matching_maps": ledger,
            "spectator_derivatives": "derivatives on the separate free gamma5 are external source factors and are never reassigned to a propagator",
            "geometric_period_relation_not_assumed": True,
            "implementation_sha256": {str(path.relative_to(jets.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = jets.base._jsonable(certify())
    target = jets.HERE / "routing_certificate.json"
    if args.write:
        target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Wrote every simultaneous-input current map on all sixteen original commutator matchings.")
    elif not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != result:
        raise ArithmeticError("the complete regenerated current-routing certificate differs")
    else:
        print("Verified all four complete input currents and all source rows on all sixteen original matchings.")


if __name__ == "__main__":
    main()
