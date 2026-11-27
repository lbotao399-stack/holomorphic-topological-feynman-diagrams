"""Complete three-loop abelian quartic product with a differentiated input.

Interaction: -i*g*c*sum(q_i*gamma_i*beta_i) + gamma_1*...*gamma_4.
Inputs: A=beta_1*(partial_{+a} beta_2)*gamma_5, B=b*beta_3;
gamma_5 is free and neutral. At order g**2*hbar**3 every surviving field
is retained: c, gamma_4, and gamma_5.
The four charges are symbolic in the calculation; physical choices obey
sum(q_i)=sum(q_i**3)=0. The complete graph sum is evaluated, including
left-derivative signs and the current on the actually differentiated edge.

The first-current polynomial has degree at most four in the gamma_5 source.
Fifteen exact source assignments therefore determine it after an invertible
GL(2) normalization of the gamma_4 and c source spinors. Covariant restoration
retains both spinor coordinates and extends polynomially to dependent sources.
See the generated certificate for every Wick matching and every sample.

Normalization from sources/LLC_main.tex, eq-2-31, eq-W-jet-edge, eq-D-lambda:
  fixed W slot: 3/3! = 1/2; I=K at three loops;
  graph sign -1; (-i)**2=-1; differentiated head A contributes -2*y.
Thus multiply the raw htk current polynomial by -g**2/(8*pi**2)**3,
then set all three sources to -partial/2 on their own surviving fields.
Since the polynomial has degree seven, equivalently multiply its expression
in physical partial derivatives by g**2/(128*(8*pi**2)**3).

The default command compares the complete regenerated exact certificate.
Run with --write to replace it after an intentional change. No table is built.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk.graph import Graph
from htk.forest_periods import CoreCurrentPeriods


def _parity(order):
    return (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:])


def _admissibility(edges):
    if len(set(edges)) != len(edges):
        return {"reason": "parallel edges"}
    for size in range(2, 6):
        for vertices in combinations(range(5), size):
            count = sum(a in vertices and b in vertices for a, b in edges)
            if count > 2 * size - 3:
                return {"reason": "Laman subset", "vertices": vertices, "edges": count}
    return None


def _fields(first, second, W_vertex):
    matter_vertices = tuple(v for v in range(3) if v != W_vertex)
    result = []
    for vertex in range(3):
        if vertex == W_vertex:
            result += [("gamma", flavor, vertex) for flavor in range(1, 5)]
        else:
            flavor = (first, second)[matter_vertices.index(vertex)]
            result += [("c", 0, vertex), ("gamma", flavor, vertex), ("beta", flavor, vertex)]
    return result + [("beta", 1, 3), ("beta", 2, 3), ("gamma5", 5, 3), ("b", 0, 4), ("beta", 3, 4)]


def wick_matchings(W_vertex=0):
    """Exhaust all beta-gamma and b-c matchings before the Laman test.

    Every field is an identifiable (species, flavor, original vertex) triple.
    The sign trace applies the manuscript's left derivatives in reverse edge
    order. No gauge invariance identity or period cancellation enters here.
    """
    if W_vertex not in (0, 1, 2):
        raise ValueError("W_vertex must be an interaction slot 0, 1, or 2")
    valid, rejected = [], []
    matter_vertices = tuple(v for v in range(3) if v != W_vertex)
    for first in range(1, 5):
        for second in range(1, 5):
            fields = _fields(first, second, W_vertex)
            beta = tuple(field for field in fields if field[0] == "beta")
            gamma = tuple(field for field in fields if field[0] == "gamma")

            def match(index, used, pairs):
                if index == len(beta):
                    yield pairs
                    return
                left = beta[index]
                for right in gamma:
                    if right not in used and left[1] == right[1] and left[2] != right[2]:
                        yield from match(index + 1, used | {right}, pairs + ((left, right),))

            for matter_pairs in match(0, set(), ()):
                for ghost_vertex in matter_vertices:
                    pairs = matter_pairs + ((("b", 0, 4), ("c", 0, ghost_vertex)),)
                    edges = tuple(sorted(tuple(sorted((a[2], b[2]))) for a, b in pairs) + [(3, 4)])
                    failure = _admissibility(edges)
                    if failure is not None:
                        rejected.append({"matter_flavors": (first, second), "pairs": pairs,
                                         "edges": edges, "vanishing_witness": failure})
                        continue
                    remaining = list(fields)
                    sign, trace = 1, []
                    for edge in reversed(edges[:-1]):
                        pair = next((a, b) for a, b in pairs if tuple(sorted((a[2], b[2]))) == edge)
                        odd = next(field for field in pair if field[0] in ("beta", "c"))
                        preceding = sum(field[0] in ("beta", "c")
                                        for field in remaining[:remaining.index(odd)])
                        pairing_sign = -1 if odd[0] == "beta" else 1
                        sign *= (-1) ** preceding * pairing_sign
                        trace.append({"edge": edge, "odd_derivative": odd,
                                      "preceding_odd_fields": preceding, "pairing_sign": pairing_sign})
                        for field in pair:
                            remaining.remove(field)
                    if set((s, i) for s, i, v in remaining) != {("c", 0), ("gamma", 4), ("gamma5", 5)}:
                        raise ArithmeticError("a surviving field was lost")
                    differentiated = ("beta", 2, 3)
                    differentiated_edge = next(tuple(sorted((a[2], b[2]))) for a, b in pairs
                                               if differentiated in (a, b))
                    if differentiated_edge[1] != 3:
                        raise ArithmeticError("the input derivative is not at the oriented head")
                    valid.append({"W_vertex": W_vertex, "matter_flavors": (first, second),
                                  "matter_vertices": matter_vertices, "initial_fields": fields,
                                  "pairs": pairs, "ordered_edges": edges, "wick_sign": sign,
                                  "left_derivative_trace": trace, "surviving_fields": tuple(remaining),
                                  "current_edge": edges.index(differentiated_edge),
                                  "derivative_factor": -2})
    return tuple(valid), tuple(rejected)


def _canonical(record):
    """Only relabel the two identical even matter slots; retain every current."""
    W_vertex = record["W_vertex"]
    m0, m1 = record["matter_vertices"]
    candidates = []
    for order in ((m0, m1), (m1, m0)):
        mapping = {W_vertex: 0, order[0]: 1, order[1]: 2, 3: 3, 4: 4}
        oriented = tuple((mapping[a], mapping[b]) for a, b in record["ordered_edges"])
        unordered = tuple(tuple(sorted(edge)) for edge in oriented)
        edges = tuple(sorted(unordered))
        permutation = tuple(edges.index(edge) for edge in unordered)
        locations = tuple(mapping[next(v for species, flavor, v in record["surviving_fields"]
                                           if species == name)] for name in ("gamma", "c", "gamma5"))
        edge = record["current_edge"]
        orientation = 1 if oriented[edge][0] < oriented[edge][1] else -1
        # D=2 vertex relabelling/reversal signs are +1; the edge order is odd.
        coefficient = record["wick_sign"] * _parity(permutation) * orientation
        candidates.append(((edges, locations), permutation[edge], coefficient))
    return min(candidates)


def grouped_matchings(W_vertex=0):
    valid, rejected = wick_matchings(W_vertex)
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    for record in valid:
        key, edge, sign = _canonical(record)
        charges = tuple(sorted(record["matter_flavors"]))
        groups[key][charges][edge] += sign
    clean = {}
    for key, charge_terms in groups.items():
        clean[key] = {charges: {edge: value for edge, value in terms.items() if value}
                      for charges, terms in charge_terms.items() if any(terms.values())}
    return dict(sorted(clean.items())), valid, rejected


def _periods(value):
    if hasattr(value, "items"):
        return {str(label): Fraction(coefficient) for label, coefficient in value.items() if coefficient}
    return {"1": Fraction(value)} if value else {}


def _sources(locations, vectors):
    if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
        raise ValueError("the complete gamma_4, c, and gamma_5 source spinors are required")
    if any(isinstance(value, (float, complex)) for vector in vectors for value in vector):
        raise TypeError("exact rational source values are required")
    rows = [[Fraction(0), Fraction(0)] for _ in range(5)]
    for vertex, vector in zip(locations, vectors):
        for coordinate in range(2):
            rows[vertex][coordinate] += Fraction(vector[coordinate])
    for coordinate in range(2):
        rows[4][coordinate] = -sum(rows[v][coordinate] for v in range(4))
    return tuple(map(tuple, rows))


def evaluate(vectors, *, check=1):
    """Return every exact graph moment and the complete symbolic charge sum."""
    groups, valid, rejected = grouped_matchings()
    total = defaultdict(lambda: [defaultdict(Fraction), defaultdict(Fraction)])
    evaluations = []
    for (edges, locations), charge_terms in groups.items():
        sources = _sources(locations, vectors)
        engine = CoreCurrentPeriods(Graph(5, edges), sources, marked=(3, 4), check=check)
        moments = {}
        try:
            for edge in sorted({edge for terms in charge_terms.values() for edge in terms}):
                spinor = []
                for coordinate in range(2):
                    exponent = [0] * 14
                    exponent[2 * edge + coordinate] = 1
                    spinor.append(_periods(engine.integrate({tuple(exponent): Fraction(1)}, allow_series=False)))
                moments[edge] = tuple(spinor)
            for charges, terms in charge_terms.items():
                # Keep charge monomials whose value vanishes at this source.
                # The interpolation lower set includes such degenerate nodes.
                total[charges]
                for edge, sign in terms.items():
                    for coordinate in range(2):
                        for label, coefficient in moments[edge][coordinate].items():
                            total[charges][coordinate][label] += sign * coefficient
        finally:
            engine.close()
        evaluations.append({"ordered_edges": edges, "survivor_vertices_gamma4_c_gamma5": locations,
                            "full_balanced_sources": sources,
                            "symbolic_charge_current_terms": [
                                {"charge_pair": charges, "edge_coefficients": dict(terms)}
                                for charges, terms in sorted(charge_terms.items())],
                            "full_first_current_moments": moments})
    clean = {charges: tuple({label: coefficient for label, coefficient in coordinate.items() if coefficient}
                           for coordinate in spinor)
             for charges, spinor in total.items()}
    return {"source_spinors_gamma4_c_gamma5": tuple(map(tuple, vectors)), "graphs": evaluations,
            "symbolic_charge_coefficients": clean}


def specialize(evaluation, charges):
    if len(charges) != 4:
        raise ValueError("four matter charges are required")
    result = [defaultdict(Fraction), defaultdict(Fraction)]
    for (first, second), spinor in evaluation["symbolic_charge_coefficients"].items():
        for coordinate in range(2):
            for label, coefficient in spinor[coordinate].items():
                result[coordinate][label] += Fraction(charges[first - 1]) * Fraction(charges[second - 1]) * coefficient
    return tuple({label: coefficient for label, coefficient in entry.items() if coefficient} for entry in result)


def _binomial_polynomial(degree):
    coefficients = [Fraction(1)]
    for root in range(degree):
        update = [Fraction(0)] * (len(coefficients) + 1)
        for power, coefficient in enumerate(coefficients):
            update[power] -= root * coefficient
            update[power + 1] += coefficient
        coefficients = update
    return tuple(coefficient / factorial(degree) for coefficient in coefficients)


def interpolate(samples, degree=4):
    """Exact Newton differences on the entire two-variable lower set."""
    nodes = tuple((a, b) for a in range(degree + 1) for b in range(degree + 1 - a))
    if set(samples) != set(nodes):
        raise ValueError("the complete total-degree interpolation grid is required")
    output = defaultdict(Fraction)
    for a, b in nodes:
        difference = sum((-1) ** (a + b - i - j) * comb(a, i) * comb(b, j) * samples[i, j]
                         for i in range(a + 1) for j in range(b + 1))
        for i, left in enumerate(_binomial_polynomial(a)):
            for j, right in enumerate(_binomial_polynomial(b)):
                output[i, j] += difference * left * right
    return {exponent: value for exponent, value in output.items() if value}


def covariant_restore(polynomials):
    """Restore a complete spinor in gamma_4, c, gamma_5 and their brackets.

    Entries are (surviving-vector index, bracket exponents, coefficient), with
    brackets (<gamma_4,c>, <gamma_5,c>, <gamma_4,gamma_5>).
    Their exponents always total three. The degree-four part must be the
    gamma_5 source spinor times a homogeneous cubic.
    """
    first, second = (dict(polynomial) for polynomial in polynomials)
    terms = []
    if first.get((0, 4), 0) or second.get((4, 0), 0):
        raise ArithmeticError("the interpolated spinor violates GL(2) covariance")
    for a in range(4):
        b = 3 - a
        coefficient = first.pop((a + 1, b), 0)
        if coefficient != second.pop((a, b + 1), 0):
            raise ArithmeticError("the degree-four spinor does not point along the spectator source")
        if coefficient:
            terms.append((2, (0, a, b), coefficient))
    for vector, polynomial in enumerate((first, second)):
        for (a, b), coefficient in sorted(polynomial.items()):
            if a + b > 3:
                raise ArithmeticError("a noncovariant term remained")
            terms.append((vector, (3 - a - b, a, b), coefficient))
    return tuple(terms)


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def evaluate_covariant(terms, vectors):
    vectors = tuple(tuple(map(Fraction, vector)) for vector in vectors)
    brackets = (_bracket(vectors[0], vectors[1]), _bracket(vectors[2], vectors[1]),
                _bracket(vectors[0], vectors[2]))
    result = [Fraction(0), Fraction(0)]
    for vector, powers, coefficient in terms:
        scalar = Fraction(coefficient)
        for bracket, power in zip(brackets, powers):
            scalar *= bracket ** power
        for coordinate in range(2):
            result[coordinate] += scalar * vectors[vector][coordinate]
    return tuple(result)


def certify(*, progress=False):
    groups, matchings, rejected = grouped_matchings()
    # Independent slot enumeration makes the factor 3/3! explicit.
    for W_vertex in (1, 2):
        if grouped_matchings(W_vertex)[0] != groups:
            raise ArithmeticError("the three superpotential placements disagree")
    if len(matchings) != 20:
        raise ArithmeticError("the exhaustive Wick census changed")
    nodes = tuple((a, b) for a in range(5) for b in range(5 - a))
    samples = {}
    for index, node in enumerate(nodes):
        samples[node] = evaluate(((1, 0), (0, 1), node))
        if progress:
            print("complete source grid", index + 1, "/", len(nodes), flush=True)
    charge_pairs = sorted({pair for charge_terms in groups.values() for pair in charge_terms})
    labels = sorted({label for sample in samples.values()
                     for spinor in sample["symbolic_charge_coefficients"].values()
                     for coordinate in spinor for label in coordinate} | {"z3"})
    covariants = {}
    for charges in charge_pairs:
        covariants[charges] = {}
        for label in labels:
            polynomials = tuple(interpolate({node: sample["symbolic_charge_coefficients"][charges][coordinate].get(label, 0)
                                              for node, sample in samples.items()}) for coordinate in range(2))
            terms = covariant_restore(polynomials)
            if terms:
                covariants[charges][label] = terms
    replays = []
    for vectors, charges in (
            (((1, 2), (3, -1), (-2, 4)), (1, 1, -1, -1)),
            (((2, -1), (1, 3), (4, -2)), (1, 2, -1, -2)),
            (((1, 2), (2, 4), (3, -1)), (1, -1, 3, -3)),
            (((1, 2), (3, -1), (-2, 4)), (1, -1, 0, 0)),
            (((1, 2), (0, 0), (-2, 4)), (1, -1, 0, 0))):
        if sum(charges) or sum(charge ** 3 for charge in charges):
            raise ArithmeticError("a replay theory has a gauge anomaly or a noninvariant superpotential")
        direct = evaluate(vectors)
        for pair in charge_pairs:
            restored = [defaultdict(Fraction), defaultdict(Fraction)]
            for label, terms in covariants[pair].items():
                value = evaluate_covariant(terms, vectors)
                for coordinate in range(2):
                    if value[coordinate]:
                        restored[coordinate][label] = value[coordinate]
            if tuple(map(dict, restored)) != direct["symbolic_charge_coefficients"][pair]:
                raise ArithmeticError("the independent complete-source replay disagrees")
        replays.append({"charges": charges, "evaluation": direct, "complete_raw_spinor": specialize(direct, charges)})
    implementation = (Path(__file__).resolve(), ROOT / "code/htk/forest_periods.py",
                      ROOT / "code/htk/_forest_period_worker.py", ROOT / "code/htk/conventions.py")
    no_undifferentiated_ghost = all(vector == 1 or powers[0] + powers[1] > 0
                                  for channels in covariants.values() for terms in channels.values()
                                  for vector, powers, coefficient in terms)
    if not no_undifferentiated_ghost:
        raise ArithmeticError("an undifferentiated surviving ghost remained")
    return {"schema": 1, "theory": "U(1), four interacting chiral fields, W=gamma1*gamma2*gamma3*gamma4, neutral free gamma5",
            "inputs": ["beta1*(partial_{+a} beta2)*gamma5", "b*beta3"],
            "all_surviving_fields": ["gamma4", "c", "gamma5"],
            "scope": "complete g^2*hbar^3 cochain coefficient for these inputs; no general rationality theorem or BPS nontriviality claim",
            "implementation_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in implementation},
            "normalization": {"W_placements_over_factorial": "1/2", "I_over_K": 1,
                              "graph_sign": -1, "two_matter_vertex_i_factors": -1,
                              "input_head_derivative_factor": -2,
                              "physical_coefficient_from_raw_sources": "-g^2/(8*pi^2)^3",
                              "source_substitution": "lambda_gamma4,c,gamma5 = -partial_gamma4,c,gamma5/2",
                              "physical_coefficient_from_covariant_derivatives": "g^2/(128*(8*pi^2)^3)"},
            "W_slot_equality": [0, 1, 2], "complete_fixed_W_matchings": matchings,
            "excluded_matchings_with_witnesses": rejected,
            "source_polynomial_certificate": {"source_degree": 7, "spectator_source_degree_bound": 4,
                                               "normalization_sources": [[1, 0], [0, 1]],
                                               "degree_bound_reason": "a source at A contributes rank at most n=3 to the current determinant; the one inserted current adds at most one",
                                               "restoration": "GL(2) covariance: det(M)^3*M; polynomial continuation includes dependent sources",
                                               "complete_lower_grid": nodes,
                                               "samples": [samples[node] for node in nodes]},
            "covariant_bracket_order": ["<lambda_gamma4,lambda_c>", "<lambda_gamma5,lambda_c>", "<lambda_gamma4,lambda_gamma5>"],
            "covariant_vector_order": ["lambda_gamma4", "lambda_c", "lambda_gamma5"],
            "complete_covariant_charge_polynomial": covariants,
            "independent_full_source_replays": replays,
            "vanishes_identically_at_lambda_c_zero": no_undifferentiated_ghost,
            "zeta3_cancels_identically": all("z3" not in terms for terms in covariants.values())}


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {(",".join(map(str, key)) if isinstance(key, tuple) else str(key)): _jsonable(item)
                for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    certificate = _jsonable(certify(progress=True))
    print("complete symbolic-charge/source zeta3 cancellation:", certificate["zeta3_cancels_identically"], flush=True)
    for replay in certificate["independent_full_source_replays"]:
        print("charges", replay["charges"], "complete raw spinor", replay["complete_raw_spinor"], flush=True)
    target = HERE / "quartic_wick_certificate.json"
    if args.write:
        target.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print("Wrote", target.name, flush=True)
    else:
        if not target.is_file() or json.loads(target.read_text(encoding="utf-8")) != certificate:
            raise ArithmeticError("the complete regenerated certificate differs from the stored certificate")
        print("Verified the complete stored certificate.", flush=True)


if __name__ == "__main__":
    main()
