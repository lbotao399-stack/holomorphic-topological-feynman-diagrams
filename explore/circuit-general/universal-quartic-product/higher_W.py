"""Complete higher-degree topology classification for the same input family.

The superpotential is a polarized pure-matter polynomial. Its first three
polarization slots are 1,2,3; further slots are 4,6,7,... so that gamma5
remains the free neutral spectator in beta1*(partial beta2)*gamma5. These
slot labels are computational labels, not an assumption of distinct
physical matter species. Diagonal specialization, with the exact rational
polarization normalization and all differentiated-slot assignments, gives
the original polynomial. All surviving fields and their original vertices
are retained. Setting gamma5 to 1 is allowed.

For any degree >=3, all valid ordered Wick matchings are one of:
  * 16 members of the universal representation-commutator cancellation;
  * 8 source-free cubic single-star contractions to an SP graph;
  * 4*(degree-3) generic source-complete SP matchings.

The finite all-degree argument is not a numerical extrapolation. If both
matter vertices act on flavors other than the three input beta flavors,
the two beta fields at A both contract directly to W, giving a forbidden
parallel edge. Thus at most one extra W flavor can occur among the matter
vertices. Relabeling that one flavor gives exactly the degree-four ledger;
the no-extra-flavor case is the degree-three ledger. Extra uncontracted W
fields only add their sources at the original W vertex. The eight single
star centres are cubic vertices 1 or 2, whose source remains exactly zero.

The resulting rationality statement is for this input family and pure
W(gamma), at first derivative order. It is not an all-input or W(gamma,w)
claim. See DERIVATION.md for the separate all-loop termination argument.
"""
from collections import Counter
from itertools import product
import argparse
import hashlib
import json

import quartic_tensor as universal
from htk.graph import Graph
from htk.spred import plan_generic, NotSeriesParallel
from htk.star import find_single, SingleMove


def polarized_flavors(degree):
    if isinstance(degree, bool) or not isinstance(degree, int) or degree < 3:
        raise ValueError("the polarized superpotential degree must be an integer at least three")
    return (1, 2, 3) + ((4,) + tuple(range(6, degree + 2)) if degree >= 4 else ())


def matchings(degree):
    """Independent exhaustive enumeration with every surviving W field retained."""
    W_flavors = polarized_flavors(degree)
    valid, rejected = [], []
    for kinds in product(("M", "G"), repeat=2):
        for flavors in product(*(W_flavors if kind == "M" else (0,) for kind in kinds)):
            fields = [("gamma", flavor, 0, 0) for flavor in W_flavors]
            for vertex, (kind, flavor) in enumerate(zip(kinds, flavors), 1):
                if kind == "M":
                    fields += [("c", 0, vertex, 0), ("gamma", flavor, vertex, 0), ("beta", flavor, vertex, 0)]
                else:
                    fields += [("b", 0, vertex, 0), ("c", 0, vertex, 0), ("c", 0, vertex, 1)]
            fields += [("beta", 1, 3, 0), ("beta", 2, 3, 0), ("gamma5", 5, 3, 0),
                       ("b", 0, 4, 0), ("beta", 3, 4, 0)]
            antifields = tuple(field for field in fields if field[0] in ("b", "beta"))
            ordinary = tuple(field for field in fields if field[0] in ("c", "gamma"))

            def recurse(index, used, pairs):
                if index == len(antifields):
                    yield pairs
                    return
                left = antifields[index]
                for right in ordinary:
                    if right in used or left[2] == right[2]:
                        continue
                    if left[0] == "b" and right[0] != "c":
                        continue
                    if left[0] == "beta" and (right[0] != "gamma" or left[1] != right[1]):
                        continue
                    yield from recurse(index + 1, used | {right}, pairs + ((left, right),))

            for pairs in recurse(0, set(), ()):
                edges = tuple(sorted(tuple(sorted((left[2], right[2]))) for left, right in pairs) + [(3, 4)])
                record = {"cubic_types": kinds, "matter_flavors": flavors, "initial_fields": tuple(fields),
                          "pairs": pairs, "ordered_edges": edges}
                failure = universal.base._admissibility(edges)
                if failure is not None:
                    rejected.append(record | {"vanishing_witness": failure})
                    continue
                remaining, sign, trace = list(fields), 1, []
                for edge in reversed(edges[:-1]):
                    pair = next((a, b) for a, b in pairs if tuple(sorted((a[2], b[2]))) == edge)
                    odd = next(field for field in pair if field[0] in ("beta", "c"))
                    preceding = sum(field[0] in ("beta", "c") for field in remaining[:remaining.index(odd)])
                    pairing_sign = -1 if odd[0] == "beta" else 1
                    sign *= (-1) ** preceding * pairing_sign
                    trace.append({"edge": edge, "odd_derivative": odd, "preceding_odd_fields": preceding,
                                  "pairing_sign": pairing_sign})
                    for field in pair:
                        remaining.remove(field)
                current_edge = next(tuple(sorted((a[2], b[2]))) for a, b in pairs if ("beta", 2, 3, 0) in (a, b))
                if current_edge[1] != 3:
                    raise ArithmeticError("the original head derivative was not preserved")
                survivors = Counter(field[0] for field in remaining)
                if survivors != Counter({"gamma": degree - 3, "c": 1, "gamma5": 1}):
                    raise ArithmeticError("the complete output field inventory changed")
                valid.append(record | {"wick_sign": sign, "left_derivative_trace": trace,
                                       "surviving_fields": tuple(remaining), "current_edge": edges.index(current_edge),
                                       "derivative_factor": -2})
    return tuple(valid), tuple(rejected)


def classify(record):
    """Keep all sources symbolic; contract only demonstrably empty cubic vertices."""
    graph = Graph(5, record["ordered_edges"])
    try:
        plan = plan_generic(graph)
    except NotSeriesParallel:
        pass
    else:
        return {"kind": "SP", "generic_SP_plan": plan,
                "source_support": sorted({field[2] for field in record["surviving_fields"]})}
    occupied = {field[2] for field in record["surviving_fields"]}
    for choice in find_single(graph):
        if choice[0] in (1, 2) and choice[0] not in occupied:
            move = SingleMove(graph, *choice)
            return {"kind": "single_to_SP", "source_free_cubic_centre": choice[0], "single_move": choice,
                    "remaining_original_sources": tuple(record["surviving_fields"]),
                    "reduced_edges": move.red.edges, "generic_SP_plan": plan_generic(move.red)}
    word, factor = universal.matching_action(record)
    allowed = [("MM", ((flavor, order),)) for flavor in (1, 2) for order in (("A", "B"), ("B", "A"))]
    allowed += [("fR", flavor) for flavor in (1, 2)]
    if word not in allowed:
        raise ArithmeticError("an unresolved higher-degree tensor word appeared")
    W_survivors = tuple(field for field in record["surviving_fields"] if field[0] == "gamma")
    if any(field[2] != 0 for field in W_survivors):
        raise ArithmeticError("a surviving W source moved in a commutator matching")
    return {"kind": "commutator", "generator_word": word, "color_prefactor": factor,
            "all_surviving_W_fields_at_vertex_zero": W_survivors,
            "complete_W_source": "sum of the independent sources of all listed surviving W fields",
            "ghost_vertex": next(field[2] for field in record["surviving_fields"] if field[0] == "c"),
            "spectator_vertex": 3}


def _commutator_signature(record):
    return (record["cubic_types"], record["matter_flavors"], record["pairs"], record["ordered_edges"],
            record["wick_sign"], record["current_edge"],
            next(field[2] for field in record["surviving_fields"] if field[0] == "c"))


def certify():
    replays = []
    reference = None
    representative_ledgers = []
    for degree in range(3, 10):
        records, rejected = matchings(degree)
        classified = [(record, classify(record)) for record in records]
        counts = Counter(result["kind"] for record, result in classified)
        if counts != Counter({"SP": 4 * (degree - 3), "single_to_SP": 8, "commutator": 16}):
            raise ArithmeticError("the complete higher-degree census changed")
        signatures = tuple(sorted(_commutator_signature(record) for record, result in classified
                                  if result["kind"] == "commutator"))
        if reference is None:
            reference = signatures
        elif signatures != reference:
            raise ArithmeticError("a higher-degree source-complete Wick sign or current routing changed")
        for record in records:
            if sum(flavor > 3 for kind, flavor in zip(record["cubic_types"], record["matter_flavors"])
                   if kind == "M") > 1:
                raise ArithmeticError("two extra matter flavors evaded the parallel-edge obstruction")
        if degree in (3, 4):
            representative_ledgers.append({"degree": degree,
                                           "all_valid_matchings": [record | {"exact_reduction": result}
                                                                    for record, result in classified],
                                           "all_excluded_matchings_with_witnesses": rejected})
        replays.append({"degree": degree, "W_flavors": polarized_flavors(degree),
                        "complete_matching_count": len(records), "classification_counts": dict(counts),
                        "all_original_commutator_signs_and_currents_match": True})
    paths = (universal.HERE / "higher_W.py", universal.HERE / "quartic_tensor.py",
             universal.ROOT / "code/htk/spred.py", universal.ROOT / "code/htk/star.py")
    return {"schema": 1,
            "scope": "arbitrary polarized pure-matter W, the fixed three-beta/one-b first-current input family, all source spinors retained",
            "spectator": "gamma5 is a free neutral field absent from W; polarization slots beyond 1,2,3,4 are labeled 6,7,... and do not denote additional physical species",
            "undoing_polarization": "The symmetric polarized tensor is normalized by 1/(degree!) times the fully multilinear coefficient. Its diagonal gives the original homogeneous polynomial. Differentiating the diagonal sums all differentiated-slot assignments with their exact integer multiplicities. The cancellation is linear in the W tensor and holds before every such specialization.",
            "implementation_sha256": {str(path.relative_to(universal.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in paths},
            "finite_all_degree_argument": "At most one matter vertex can act on a W flavor outside the three input flavors: otherwise both A betas contract directly to W and produce a forbidden parallel W--A edge. The no-extra-flavor and one-extra-flavor cases are exactly the two complete representative ledgers. Every additional uncontracted W field retains its source at vertex0.",
            "representative_complete_ledgers": representative_ledgers,
            "independent_higher_degree_replays": replays,
            "all_orders_for_these_inputs": {"antifield_bound": "n+r<=4 for n interactions and r pure-W vertices",
                                             "cubic_only_even_field_bound": "2*n<=n+1, since neutral gamma5 cannot contract",
                                             "third_order": "exactly one W and two cubic vertices; complete representation-commutator cancellation",
                                             "orders_at_most_two": "all geometric graphs are series-parallel",
                                             "orders_at_least_four": "zero by the two field-count bounds",
                                             "coefficient_ring": "rational coefficients in the original interaction coefficients and W tensor, after the universal (8*pi^2)^(-n) normalization",
                                             "not_claimed": "an explicit evaluation of every lower-order coefficient; all input families; W(gamma,w); quantum-closed or nontrivial observable classes"}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    certificate = universal.base._jsonable(certify())
    target = universal.HERE / "higher_W_certificate.json"
    if args.write:
        target.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print("Wrote the complete higher-degree classification certificate.")
    else:
        if not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != certificate:
            raise ArithmeticError("the complete regenerated higher-degree classification differs")
        print("Verified the complete representative ledgers and all higher-degree source/current replays.")


if __name__ == "__main__":
    main()
