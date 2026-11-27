"""Complete SU(2) adjoint quartic product, with all color and current factors.

W=(gamma_1.gamma_3)*(gamma_2.gamma_4), with a free neutral gamma_5.
A=(beta_1.partial_{+a}beta_2)*gamma_5, B=b.beta_3.
The dots contract the complete adjoint metric. Both inputs are invariant
under constant gauge transformations; Q-closure is not assumed.

t_a[j,i]=-i*epsilon(a,j,i), so the matter interaction has real coefficient
+g*epsilon(a,i,j)*c^a*gamma_j*beta^i. The gauge interaction is
-g*epsilon(a,b,c)*b_a*c^b*c^c/2. Every epsilon and every metric is contracted,
and all nine surviving (c^a,gamma_4^b) color entries are checked before
factoring out delta_ab. The two gauge-vertex case is exhaustively excluded
by its two parallel W--A edges, rather than omitted.

The fixed-W-slot sum contains all ordered MM, MG, GM, and GG choices.
Its common physical multiplier is +g^2/(8*pi^2)^3, since
(-1)^3 * (3/3!) * (-2) = +1 and the real vertex coefficients are included.
After lambda=-partial/2 on the three surviving fields this becomes
-g^2/[128*(8*pi^2)^3]. The coefficient is at hbar^3.

The same complete fifteen-node source interpolation as quartic_wick.py
proves the result for arbitrary sources. No period relation or coefficient
independence is presumed. The default CLI replays the entire certificate;
--write regenerates it after an intentional change.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
from pathlib import Path

import quartic_wick as base
from htk.graph import Graph
from htk.forest_periods import CoreCurrentPeriods

HERE = Path(__file__).resolve().parent
PAIRING = ((1, 3), (2, 4))


def epsilon(a, b, c):
    return 0 if len({a, b, c}) < 3 else (-1) ** ((a > b) + (a > c) + (b > c))


def _fields(kinds, flavors):
    fields = [("gamma", flavor, 0, 0) for flavor in range(1, 5)]
    for vertex, (kind, flavor) in enumerate(zip(kinds, flavors), 1):
        if kind == "M":
            fields += [("c", 0, vertex, 0), ("gamma", flavor, vertex, 0), ("beta", flavor, vertex, 0)]
        else:
            fields += [("b", 0, vertex, 0), ("c", 0, vertex, 0), ("c", 0, vertex, 1)]
    return tuple(fields + [("beta", 1, 3, 0), ("beta", 2, 3, 0), ("gamma5", 5, 3, 0),
                           ("b", 0, 4, 0), ("beta", 3, 4, 0)])


def uncolored_matchings():
    """Exhaust every cubic-vertex type, flavor, and ordered Wick matching."""
    valid, rejected = [], []
    for kinds in product(("M", "G"), repeat=2):
        for flavors in product(*(range(1, 5) if kind == "M" else (0,) for kind in kinds)):
            fields = _fields(kinds, flavors)
            antifields = tuple(field for field in fields if field[0] in ("b", "beta"))
            ordinary = tuple(field for field in fields if field[0] in ("c", "gamma"))

            def match(index, used, pairs):
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
                    yield from match(index + 1, used | {right}, pairs + ((left, right),))

            for pairs in match(0, set(), ()):
                edges = tuple(sorted(tuple(sorted((a[2], b[2]))) for a, b in pairs) + [(3, 4)])
                failure = base._admissibility(edges)
                if failure is not None:
                    rejected.append({"cubic_types": kinds, "matter_flavors": flavors,
                                     "pairs": pairs, "ordered_edges": edges, "vanishing_witness": failure})
                    continue
                remaining, sign, trace = list(fields), 1, []
                for edge in reversed(edges[:-1]):
                    pair = next((a, b) for a, b in pairs if tuple(sorted((a[2], b[2]))) == edge)
                    odd = next(field for field in pair if field[0] in ("beta", "c"))
                    preceding = sum(field[0] in ("beta", "c") for field in remaining[:remaining.index(odd)])
                    pairing_sign = -1 if odd[0] == "beta" else 1
                    sign *= (-1) ** preceding * pairing_sign
                    trace.append({"edge": edge, "odd_derivative": odd,
                                  "preceding_odd_fields": preceding, "pairing_sign": pairing_sign})
                    for field in pair:
                        remaining.remove(field)
                if set((species, flavor) for species, flavor, vertex, slot in remaining) != {
                        ("c", 0), ("gamma", 4), ("gamma5", 5)}:
                    raise ArithmeticError("a surviving output field was omitted")
                edge = next(tuple(sorted((a[2], b[2]))) for a, b in pairs if ("beta", 2, 3, 0) in (a, b))
                if edge[1] != 3:
                    raise ArithmeticError("the differentiated input is not at the oriented head")
                valid.append({"cubic_types": kinds, "matter_flavors": flavors, "initial_fields": fields,
                              "pairs": pairs, "ordered_edges": edges, "wick_sign": sign,
                              "left_derivative_trace": trace, "surviving_fields": tuple(remaining),
                              "current_edge": edges.index(edge), "derivative_factor": -2})
    return tuple(valid), tuple(rejected)


def complete_color_tensor(record):
    """Contract every internal color and compute all nine output entries."""
    fields = tuple(field for field in record["initial_fields"] if field[0] != "gamma5")
    parent = {field: field for field in fields}

    def find(field):
        while parent[field] != field:
            parent[field] = parent[parent[field]]
            field = parent[field]
        return field

    def join(left, right):
        parent[find(right)] = find(left)

    metrics = [(('gamma', left, 0, 0), ('gamma', right, 0, 0)) for left, right in PAIRING]
    metrics += [(('beta', 1, 3, 0), ('beta', 2, 3, 0)), (('b', 0, 4, 0), ('beta', 3, 4, 0))]
    metrics += list(record["pairs"])
    for left, right in metrics:
        join(left, right)
    roots = sorted({find(field) for field in fields})
    ghost = next(field for field in record["surviving_fields"] if field[0] == "c")
    gamma = next(field for field in record["surviving_fields"] if field[:2] == ("gamma", 4))
    ghost_root, gamma_root = find(ghost), find(gamma)
    epsilons, prefactor = [], F(1)
    for vertex, (kind, flavor) in enumerate(zip(record["cubic_types"], record["matter_flavors"]), 1):
        if kind == "M":
            tensor = (('c', 0, vertex, 0), ('beta', flavor, vertex, 0), ('gamma', flavor, vertex, 0))
        else:
            tensor = (('b', 0, vertex, 0), ('c', 0, vertex, 0), ('c', 0, vertex, 1))
            prefactor *= F(-1, 2)
        epsilons.append(tuple(find(field) for field in tensor))
    matrix = []
    for out_c in range(3):
        row = []
        for out_gamma in range(3):
            if ghost_root == gamma_root and out_c != out_gamma:
                row.append(F(0))
                continue
            fixed = {ghost_root: out_c, gamma_root: out_gamma}
            free = tuple(root for root in roots if root not in fixed)
            value = 0
            for colors in product(range(3), repeat=len(free)):
                assignment = fixed | dict(zip(free, colors))
                term = 1
                for tensor in epsilons:
                    term *= epsilon(*(assignment[root] for root in tensor))
                value += term
            row.append(prefactor * value)
        matrix.append(tuple(row))
    if any(matrix[a][b] != (matrix[0][0] if a == b else 0) for a in range(3) for b in range(3)):
        raise ArithmeticError("the complete output tensor is not the stated invariant metric")
    return {"complete_output_color_matrix": tuple(matrix), "common_metric_coefficient": matrix[0][0],
            "all_metric_pairs": tuple(metrics), "identified_color_roots": tuple(roots),
            "epsilon_color_roots": tuple(epsilons), "gauge_vertex_prefactor": prefactor}


def _canonical(record):
    candidates = []
    for mapping in ((0, 1, 2, 3, 4), (0, 2, 1, 3, 4)):
        oriented = tuple((mapping[a], mapping[b]) for a, b in record["ordered_edges"])
        edges = tuple(sorted(tuple(sorted(edge)) for edge in oriented))
        permutation = tuple(edges.index(tuple(sorted(edge))) for edge in oriented)
        locations = tuple(mapping[next(vertex for species, flavor, vertex, slot in record["surviving_fields"]
                                       if species == name)] for name in ("gamma", "c", "gamma5"))
        edge = record["current_edge"]
        orientation = 1 if oriented[edge][0] < oriented[edge][1] else -1
        coefficient = record["wick_sign"] * record["color"]["common_metric_coefficient"]
        coefficient *= base._parity(permutation) * orientation
        candidates.append(((edges, locations), permutation[edge], coefficient))
    return min(candidates)


def grouped_matchings():
    uncolored, rejected = uncolored_matchings()
    records = tuple(record | {"color": complete_color_tensor(record)} for record in uncolored)
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(F)))
    for record in records:
        key, edge, coefficient = _canonical(record)
        kind = "MM" if record["cubic_types"] == ("M", "M") else "MG"
        groups[key][kind][edge] += coefficient
    return dict(sorted(groups.items())), records, rejected


def evaluate(vectors, *, check=1):
    groups, records, rejected = grouped_matchings()
    totals = {kind: [defaultdict(F), defaultdict(F)] for kind in ("MM", "MG", "complete")}
    evaluations = []
    for (edges, locations), kinds in groups.items():
        sources = base._sources(locations, vectors)
        engine = CoreCurrentPeriods(Graph(5, edges), sources, marked=(3, 4), check=check)
        moments = {}
        try:
            for edge in sorted({edge for terms in kinds.values() for edge in terms}):
                spinor = []
                for coordinate in range(2):
                    exponent = [0] * 14
                    exponent[2 * edge + coordinate] = 1
                    spinor.append(base._periods(engine.integrate({tuple(exponent): F(1)}, allow_series=False)))
                moments[edge] = tuple(spinor)
            for kind, terms in kinds.items():
                for edge, coefficient in terms.items():
                    for coordinate in range(2):
                        for label, value in moments[edge][coordinate].items():
                            totals[kind][coordinate][label] += coefficient * value
                            totals["complete"][coordinate][label] += coefficient * value
        finally:
            engine.close()
        evaluations.append({"ordered_edges": edges, "survivor_vertices_gamma4_c_gamma5": locations,
                            "full_balanced_sources": sources, "complete_color_weighted_currents": dict(kinds),
                            "full_first_current_moments": moments})
    clean = {kind: tuple({label: coefficient for label, coefficient in coordinate.items() if coefficient}
                        for coordinate in spinor) for kind, spinor in totals.items()}
    return {"source_spinors_gamma4_c_gamma5": tuple(map(tuple, vectors)), "graphs": evaluations,
            "exact_spinors_by_vertex_type": clean}


def certify(*, progress=False):
    groups, records, rejected = grouped_matchings()
    counts = {kind: sum(record["cubic_types"] == kind for record in records)
              for kind in product(("M", "G"), repeat=2)}
    if counts != {("M", "M"): 20, ("M", "G"): 4, ("G", "M"): 4, ("G", "G"): 0}:
        raise ArithmeticError("the exhaustive nonabelian Wick census changed")
    double_gauge = tuple(record for record in rejected if record["cubic_types"] == ("G", "G"))
    if not double_gauge or any(record["ordered_edges"].count((0, 3)) != 2 for record in double_gauge):
        raise ArithmeticError("the double-gauge exclusion does not retain both beta contractions")
    nodes = tuple((a, b) for a in range(5) for b in range(5 - a))
    samples = {}
    for index, node in enumerate(nodes):
        samples[node] = evaluate(((1, 0), (0, 1), node))
        if progress:
            print("complete nonabelian source grid", index + 1, "/", len(nodes), flush=True)
    covariants = {}
    for kind in ("MM", "MG", "complete"):
        covariants[kind] = {}
        labels = {label for sample in samples.values() for coordinate in sample["exact_spinors_by_vertex_type"][kind]
                  for label in coordinate} | {"z3"}
        for label in sorted(labels):
            polynomials = tuple(base.interpolate({node: sample["exact_spinors_by_vertex_type"][kind][coordinate].get(label, 0)
                                                   for node, sample in samples.items()}) for coordinate in range(2))
            terms = base.covariant_restore(polynomials)
            if terms:
                covariants[kind][label] = terms
    replays = []
    for vectors in (((1, 2), (3, -1), (-2, 4)), ((2, -1), (1, 3), (4, -2)),
                    ((1, 2), (2, 4), (3, -1)), ((1, 2), (0, 0), (-2, 4)),
                    ((1, 2), (3, -1), (0, 0))):
        direct = evaluate(vectors)
        for kind in ("MM", "MG", "complete"):
            restored = [defaultdict(F), defaultdict(F)]
            for label, terms in covariants[kind].items():
                spinor = base.evaluate_covariant(terms, vectors)
                for coordinate in range(2):
                    if spinor[coordinate]:
                        restored[coordinate][label] = spinor[coordinate]
            if tuple(map(dict, restored)) != direct["exact_spinors_by_vertex_type"][kind]:
                raise ArithmeticError("the independent full-spinor replay disagrees")
        replays.append(direct)
    no_undifferentiated_ghost = all(vector == 1 or powers[0] + powers[1] > 0
                                  for terms in covariants["complete"].values()
                                  for vector, powers, coefficient in terms)
    if not no_undifferentiated_ghost or "z3" in covariants["complete"]:
        raise ArithmeticError("the stated complete cancellation failed")
    implementations = (Path(__file__).resolve(), HERE / "quartic_wick.py", base.ROOT / "code/htk/forest_periods.py",
                       base.ROOT / "code/htk/_forest_period_worker.py", base.ROOT / "code/htk/conventions.py")
    return {"schema": 1, "theory": "SU(2), four adjoint chiral flavors, W=(gamma1.gamma3)*(gamma2.gamma4), free neutral gamma5",
            "inputs": ["(beta1.partial_{+a} beta2)*gamma5", "b.beta3"],
            "all_surviving_fields": ["gamma4", "c", "gamma5"],
            "output_color_tensor": "delta_ab c^a gamma4^b gamma5; all nine entries checked for every matching",
            "representation_convention": "t_a[j,i]=-i*epsilon(a,j,i)",
            "scope": "complete g^2*hbar^3 cochain coefficient for these color-singlet inputs; full local gauge invariance, Q-closure, and BPS nontriviality are not asserted",
            "implementation_sha256": {str(path.relative_to(base.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in implementations},
            "normalization": {"W_placements_over_factorial": "1/2", "I_over_K": 1,
                              "graph_sign": -1, "input_head_derivative_factor": -2,
                              "matter_real_tensor": "+epsilon(a,i,j)", "gauge_real_tensor": "-epsilon(a,b,c)/2",
                              "physical_coefficient_from_raw_sources": "+g^2/(8*pi^2)^3",
                              "source_substitution": "lambda_gamma4,c,gamma5=-partial_gamma4,c,gamma5/2",
                              "physical_coefficient_from_covariant_derivatives": "-g^2/(128*(8*pi^2)^3)"},
            "complete_fixed_W_matchings": records, "matching_counts": counts,
            "excluded_matchings_with_witnesses": rejected,
            "double_gauge_exclusion": "the two differentiated/undifferentiated beta inputs at A both contract only to W, giving parallel W--A edges",
            "source_polynomial_certificate": {"source_degree": 7, "spectator_source_degree_bound": 4,
                                               "normalization_sources": [[1, 0], [0, 1]], "complete_lower_grid": nodes,
                                               "samples": [samples[node] for node in nodes],
                                               "degree_bound_reason": "one source has rank at most n=3 in the current determinant; the one inserted current adds one",
                                               "restoration": "det(M)^3*M under invertible GL(2), followed by polynomial continuation"},
            "covariant_bracket_order": ["<lambda_gamma4,lambda_c>", "<lambda_gamma5,lambda_c>", "<lambda_gamma4,lambda_gamma5>"],
            "covariant_vector_order": ["lambda_gamma4", "lambda_c", "lambda_gamma5"],
            "complete_covariant_spinors_by_vertex_type": covariants,
            "independent_full_source_replays": replays,
            "zeta3_cancels_identically": "z3" not in covariants["complete"],
            "vanishes_identically_at_lambda_c_zero": no_undifferentiated_ghost}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    certificate = base._jsonable(certify(progress=True))
    print("complete source-polynomial zeta3 cancellation:", certificate["zeta3_cancels_identically"], flush=True)
    print("generic full spinors:", certificate["independent_full_source_replays"][0]["exact_spinors_by_vertex_type"], flush=True)
    target = HERE / "nonabelian_wick_certificate.json"
    if args.write:
        target.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print("Wrote", target.name, flush=True)
    else:
        if not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != certificate:
            raise ArithmeticError("the complete regenerated nonabelian certificate differs")
        print("Verified the complete stored nonabelian certificate.", flush=True)


if __name__ == "__main__":
    main()
