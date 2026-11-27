"""Recompute an integrated three-terminal T=1 response with every source/current.

The new response calculation is independent of the direct bulk numerator.
The scalar period tables and complete covariance constraints are declared
inputs.  The direct complete answers are consulted only after integration.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from math import prod
from pathlib import Path
import argparse
import gzip
import json
import sys
import time

from star_response import StarResponse, integer_determinant, rational, scalar
import source_orbits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMPLETION = ROOT / "explore/circuit-general/topological-completion"


def encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, separators=(",", ":")).encode()


def digest(value):
    return sha256(canonical(value)).hexdigest()


def polynomial_record(polynomial):
    rows = [[list(map(int, a)), str(rational(c))]
            for a, c in sorted(polynomial.to_dict().items())]
    return {"terms": len(rows), "sha256": digest(rows)}


def edge_actions(edges):
    edge_index = {tuple(sorted(edge)): i for i, edge in enumerate(edges)}
    result = set()
    for vertex in permutations(range(6)):
        transformed = [tuple(sorted((vertex[a], vertex[b]))) for a, b in edges]
        if all(edge in edge_index for edge in transformed):
            result.add(tuple(edge_index[edge] for edge in transformed))
    return tuple(sorted(result))


class Integration:
    def __init__(self, response):
        self.response = response
        self.actions = edge_actions(response.edges)
        self.inputs = set()
        self.table_data = {}
        self.boundaries = []
        for mask in range(1, (1 << 9)-1):
            selected = tuple(i for i in range(9) if mask & (1 << i))
            parent = list(range(6))
            def find(v):
                while v != parent[v]:
                    v = parent[v]
                return v
            rank = 0
            for e in selected:
                a, b = map(find, response.edges[e])
                if a != b:
                    parent[a] = b
                    rank += 1
            loops = len(selected)-rank
            if loops:
                self.boundaries.append((selected, loops))
        self.convergence_checks = set()

    @lru_cache(maxsize=None)
    def representative(self, alpha):
        return min(tuple(alpha[j] for j in action) for action in self.actions)

    def table(self, power):
        if power not in self.table_data:
            name = self.response.name+"_q"+str(power.numerator)+"o"+str(power.denominator)+".json"
            paths = (COMPLETION / "tables" / name,
                     ROOT / "explore/circuit-general/cores/tables" / name)
            path = next((p for p in paths if p.exists()), None)
            if path is None:
                raise FileNotFoundError("required exact scalar level: "+name)
            data = json.loads(path.read_text())
            if tuple(map(tuple, data["edges"])) != self.response.edges:
                raise ArithmeticError("scalar table edge order disagrees with the response")
            if Fraction(data["q"]) != power:
                raise ArithmeticError("scalar table pole disagrees with the complete numerator")
            self.inputs.add(path)
            self.table_data[power] = {tuple(alpha): {label: Fraction(c) for label, c in value.items()}
                                      for alpha, value in data["values"]}
        return self.table_data[power]

    def converges(self, alpha, power):
        key = alpha, power
        if key in self.convergence_checks:
            return
        if min(alpha) < 1 or sum(alpha) != 4*power:
            raise ArithmeticError("nonpositive or nonprojective scalar term")
        for selected, loops in self.boundaries:
            if power.denominator*sum(alpha[e] for e in selected) <= power.numerator*loops:
                raise ArithmeticError("an individual scalar boundary is divergent")
        self.convergence_checks.add(key)

    def integrate(self, numerator, power):
        coefficients = {}
        for exponent, coefficient in numerator.to_dict().items():
            alpha = tuple(int(x)+1 for x in exponent)
            representative = self.representative(alpha)
            self.converges(representative, power)
            coefficients[representative] = coefficients.get(representative, Fraction(0))+rational(coefficient)
        coefficients = {a: c for a, c in coefficients.items() if c}
        table = self.table(power) if coefficients else {}
        value = {}
        for alpha, coefficient in coefficients.items():
            self.converges(alpha, power)
            for label, c in table[alpha].items():
                value[label] = value.get(label, Fraction(0))+coefficient*c
        value = {label: c for label, c in value.items() if c}
        return {"pole": str(power), "complete_numerator": polynomial_record(numerator),
                "scalar_orbits": [[list(a), str(c)] for a, c in sorted(coefficients.items())],
                "value": encode(value)}

    def all_currents(self, density):
        targets = [self.integrate(density.numerator, density.power)]
        for current in density.currents:
            targets.append(self.integrate(self.response.bounded(density.numerator*current), density.power+1))
        return targets


def solve_complete_source_polynomial(covariance, samples):
    # Every unknown is retained.  The exact homogeneous source/covariance
    # certificate supplies a basis of the entire 1204-coefficient solution
    # space, rather than an ansatz chosen from observed numerical values.
    basis = [{int(k): Fraction(v) for k, v in row} for row in covariance["basis"]]
    pivots = {}
    for point, sample in zip(covariance["samples"], samples):
        source = tuple(map(Fraction, point["sources"]))
        weights4 = [prod(x**p for x, p in zip(source, exponent))
                    for exponent in covariance["monomials"]["4"]]
        weights5 = [prod(x**p for x, p in zip(source, exponent))
                    for exponent in covariance["monomials"]["5"]]
        evaluations = []
        for vector in basis:
            result = [Fraction(0)]*10
            for column, coefficient in vector.items():
                if column < len(weights4):
                    output, value = 0, weights4[column]
                else:
                    output = 1+(column-len(weights4))//len(weights5)
                    value = weights5[(column-len(weights4))%len(weights5)]
                result[output] += coefficient*value
            evaluations.append(result)
        if encode(evaluations) != point["evaluations"]:
            raise ArithmeticError("a claimed source evaluation row does not evaluate the complete basis")
        for output in range(10):
            row = {j: v[output] for j, v in enumerate(evaluations) if v[output]}
            value = {label: Fraction(c) for label, c in sample["targets"][output]["value"].items()}
            while row:
                column = max(row)
                if column not in pivots:
                    pivot = row[column]
                    pivots[column] = ({k: v/pivot for k, v in row.items()},
                                      {k: v/pivot for k, v in value.items()})
                    break
                factor = row[column]
                previous, answer = pivots[column]
                for k, v in previous.items():
                    row[k] = row.get(k, 0)-factor*v
                    if not row[k]:
                        del row[k]
                for k, v in answer.items():
                    value[k] = value.get(k, 0)-factor*v
                    if not value[k]:
                        del value[k]
            else:
                if value:
                    raise ArithmeticError("an excess exact source evaluation is inconsistent")
    if len(pivots) != len(basis):
        raise ArithmeticError("the full source/current polynomial is underdetermined")
    coefficients = [{} for _ in basis]
    for column, (row, answer) in sorted(pivots.items()):
        value = dict(answer)
        for k, v in row.items():
            if k != column:
                for label, c in coefficients[k].items():
                    value[label] = value.get(label, 0)-v*c
        coefficients[column] = {k: v for k, v in value.items() if v}
    collected = {}
    for vector, coefficient in zip(basis, coefficients):
        for column, c in vector.items():
            for label, v in coefficient.items():
                key = column, label
                collected[key] = collected.get(key, 0)+c*v
    size4, size5 = len(covariance["monomials"]["4"]), len(covariance["monomials"]["5"])
    outputs = [[] for _ in range(10)]
    for (column, label), coefficient in sorted(collected.items()):
        if not coefficient:
            continue
        if column < size4:
            output, exponent = 0, covariance["monomials"]["4"][column]
        else:
            output = 1+(column-size4)//size5
            exponent = covariance["monomials"]["5"][(column-size4)%size5]
        outputs[output].append([exponent, label, str(coefficient)])
    return outputs, len(pivots)


def evaluate_answer(outputs, source):
    result = []
    for polynomial in outputs:
        value = {}
        for exponent, label, coefficient in polynomial:
            value[label] = value.get(label, 0)+Fraction(coefficient)*prod(x**k for x, k in zip(source, exponent))
        result.append({label: c for label, c in value.items() if c})
    return encode(result)


def chart_control(response, density):
    # Response triangle resistances occupy the same three slots as the
    # removed star.  This is a positive bijection, with no missing fibre.
    a, b, c = response.star
    r = [response.x[e] for e in (a, b, c)]
    total = sum(r)
    numerators = [r[0]*r[1], r[0]*r[2], r[1]*r[2]]
    derivative = [[p.derivative(e)*total-p for e in (a, b, c)] for p in numerators]
    jacobian_numerator = response.determinant(derivative)
    if jacobian_numerator != -prod(r)*total**3:
        raise ArithmeticError("the signed three-terminal response Jacobian is wrong")
    def transform(polynomial, clearing):
        answer = response.zero
        powers = [total**k for k in range(clearing+1)]
        for exponent, coefficient in polynomial.to_dict().items():
            exponent = tuple(map(int, exponent))
            star_degree = exponent[a]+exponent[b]+exponent[c]
            if star_degree > clearing:
                raise ArithmeticError("the complete response chart needs an extra pole")
            term = response.context.constant(coefficient)*powers[clearing-star_degree]
            for e, power in enumerate(exponent):
                if e in response.star:
                    term *= numerators[response.star.index(e)]**power
                else:
                    term *= response.x[e]**power
            answer = response.bounded(answer+term)
        return answer
    denominator = transform(response.U, 2)
    triangle_U = response.divide(denominator, total)
    terminals = [next(v for v in response.edges[e] if v != response.centre)
                 for e in response.star]
    triangle_edges = [(e, response.edges[e]) for e in range(9) if e not in response.star]
    triangle_edges += [(e, pair) for e, pair in zip(response.star, combinations(terminals, 2))]
    tree_sum = response.zero
    for chosen in combinations(range(9), 4):
        parent = list(range(6))
        def find(v):
            while parent[v] != v:
                v = parent[v]
            return v
        for edge in chosen:
            u, v = map(find, triangle_edges[edge][1])
            if u == v:
                break
            parent[u] = v
        else:
            tree_sum += prod(response.x[triangle_edges[e][0]] for e in range(9) if e not in chosen)
    if triangle_U != tree_sum:
        raise ArithmeticError("the response triangle determinant normalization is wrong")
    clearing = int(2*density.power)-3
    numerator = response.bounded(prod(r)*transform(density.numerator, clearing))
    # Its total degree includes the absolute response Jacobian.  Signed
    # integration is determined by the original ordered projective form.
    expected = 6*density.power-9
    if any(sum(map(int, exponent)) != expected for exponent in numerator.to_dict()):
        raise ArithmeticError("the integrated response chart has incorrect degree")
    return {"star_edge_slots": list(response.star), "jacobian_sign": -1,
            "jacobian_absolute": "R01*R02*R12/(R01+R02+R12)^3",
            "cleared_denominator": polynomial_record(denominator),
            "triangle_U": polynomial_record(triangle_U),
            "complete_transformed_numerator": polynomial_record(numerator),
            "pole": str(density.power), "all_chart_divisions_exact": True}


def build(name="prism"):
    start = time.monotonic()
    covariance_path = COMPLETION / (name+"_source_covariance.json")
    covariance = json.loads(covariance_path.read_text())
    # Replay the full constraint rank and every basis vector.  This is a
    # declared algebraic input shared with the direct calculation, whereas
    # all density numerators and integrals below are newly reconstructed.
    sys.path.insert(0, str(COMPLETION))
    import source_covariance
    system, _, _ = source_covariance.generate_constraints(name)
    if len(system.pivots) != covariance["constraint_rank"] or system.columns != 1204:
        raise ArithmeticError("the complete covariance rank changed")
    basis = [[[int(k), str(v)] for k, v in sorted(row.items())] for row in system.basis()]
    if basis != covariance["basis"]:
        raise ArithmeticError("the full covariance solution basis changed")
    response = StarResponse(name)
    integration = Integration(response)
    samples = []
    for index, sample in enumerate(covariance["samples"]):
        density = response.density(sample["sources"])
        targets = integration.all_currents(density)
        samples.append({"source": sample["sources"], "density": polynomial_record(density.numerator),
                        "pole": str(density.power), "exact_U_divisions": density.exact_U_divisions,
                        "all_original_current_numerators": [polynomial_record(p) for p in density.currents],
                        "targets": targets})
        print(name, "response integration", index+1, "of", len(covariance["samples"]),
              "bare and all nine currents complete; seconds", round(time.monotonic()-start, 2), flush=True)
    outputs, rank = solve_complete_source_polynomial(covariance, samples)
    canonical_value = evaluate_answer(outputs, (1, 1, 0, -1, 0))[0]
    if name == "prism" and canonical_value != {"sqrtpi^4": "311/630", "sqrtpi^8": "-17/360"}:
        raise ArithmeticError("the displayed canonical prism value is wrong")
    generic_source = (1, 2, 3, 4, 5)
    generic = response.density(generic_source)
    generic_targets = integration.all_currents(generic)
    if [t["value"] for t in generic_targets] != evaluate_answer(outputs, generic_source):
        raise ArithmeticError("the extra generic source point disagrees with the full polynomial")
    controls = []
    for times, source in ((tuple(range(1, 10)), generic_source),
                          ((2, 5, 3, 7, 11, 13, 17, 19, 23), generic_source),
                          (tuple(range(1, 10)), (0, 2, 3, 4, 5))):
        control = response.contact_control(times, source)
        if not control["single_contact_coefficient"]:
            raise ArithmeticError("the retained contact control unexpectedly vanished")
        from htk import brute
        actual = brute.density(response.graph, times, [[x] for x in source], T=1, check_basic=True)
        if actual.coeff != control["complete_coefficient"] or actual.radicand != control["U"]:
            raise ArithmeticError("the full response differs from the original exact density")
        controls.append({"times": times, "sources": source, **control})
    chart = chart_control(response, generic)
    # These references are first consulted after every independent
    # response integral and the full source reconstruction are complete.
    from htk.topological_polynomial import BulkPolynomial
    reference = BulkPolynomial(response.graph, [[x] for x in generic_source], H=1, T=1)
    if generic.power != reference.power or generic.numerator.to_dict() != reference.numerator.to_dict():
        raise ArithmeticError("complete response numerator differs from direct Pfaffian")
    if any(a.to_dict() != b[0].to_dict() for a, b in zip(generic.currents, reference.current_numerators)):
        raise ArithmeticError("an original current differs from the direct computation")
    reference_path = COMPLETION / (name+"_complete_polynomials.json")
    reference_answer = json.loads(reference_path.read_text())
    if outputs != reference_answer["bare_and_all_nine_current_polynomials"]:
        raise ArithmeticError("the complete integrated response differs from the direct source polynomial")
    # Delete a different vertex, including incoming original star edges;
    # compare entire time polynomials, not selected rational-time values.
    alternate = StarResponse(name, centre=1)
    alternate_density = alternate.density(generic_source)
    if alternate_density.power != generic.power or alternate_density.numerator != generic.numerator:
        raise ArithmeticError("changing the three-terminal centre changes the full kernel")
    if any(a != b for a, b in zip(alternate_density.currents, generic.currents)):
        raise ArithmeticError("changing the centre changes an original current")
    data = {"schema": "integrated-three-terminal-markov-response-v1", "core": name,
            "H": 1, "T": 1, "edges": response.edges, "deleted_vertex": response.centre,
            "all_independent_source_vertices": list(range(5)),
            "root_source": "minus the sum of all five independent original sources",
            "source_unknown_coefficients": 1204, "constraint_rank": len(system.pivots),
            "complete_covariance_dimension": len(basis), "exact_evaluation_rank": rank,
            "all_source_bare_and_nine_current_polynomials": outputs,
            "generic_all_source_control": {"sources": generic_source, "targets": generic_targets},
            "response_chart": chart, "contact_controls": controls,
            "whole_polynomial_controls": {"direct_pfaffian": True, "all_nine_currents": True,
                                           "different_deleted_vertex": 1},
            "strict_scalar_boundary_checks": len(integration.convergence_checks),
            "response_tree_terms": len(response.U), "nonzero_contact_covariance_pairs": response.contact_pairs,
            "complete_exterior_splits": len(response.exterior), "maximum_polynomial_support": response.peak_terms,
            "all_contact_terms_retained": True, "termwise_scalar_integration_certified": True,
            "source_reduction_has_no_discarded_coefficient": True}
    inputs = integration.inputs | {covariance_path, reference_path,
                                  COMPLETION / "source_covariance.py",
                                  ROOT / "code/htk/topological_polynomial.py",
                                  ROOT / "code/htk/brute.py", ROOT / "code/htk/graph.py"}
    return data, samples, inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    orbit_certificate = source_orbits.verify(args.write)
    data, samples, inputs = build("prism")
    data["full_pointwise_source_convergence_certificate_sha256"] = digest(orbit_certificate)
    sample_bytes = canonical(samples)+b"\n"
    sample_file = HERE / "response_samples.json.gz"
    output_file = HERE / "complete_prism_response.json"
    certificate_file = HERE / "certificate.json"
    sources = inputs | {HERE / "star_response.py", HERE / "verify.py",
                        HERE / "source_orbits.py", HERE / "source_orbit_certificate.json",
                        HERE / "README.md", HERE / "DERIVATION.md"}
    certificate = {"schema": "integrated-three-terminal-response-certificate-v1",
                   "complete_result_sha256": digest(data),
                   "complete_sample_sha256": sha256(sample_bytes).hexdigest(),
                   "sample_count": len(samples), "all_integrated_outputs_per_sample": 10,
                   "every_source_and_original_current_retained": True,
                   "sources": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                               for path in sorted(sources)}}
    if args.write:
        sample_file.write_bytes(gzip.compress(sample_bytes, mtime=0))
        output_file.write_bytes(canonical(data)+b"\n")
        certificate_file.write_bytes(canonical(certificate)+b"\n")
        print("Wrote complete integrated response, all exact samples, and source bindings.", flush=True)
    else:
        if gzip.decompress(sample_file.read_bytes()) != sample_bytes:
            raise ArithmeticError("the complete response sample replay changed")
        if json.loads(output_file.read_text()) != json.loads(canonical(data)):
            raise ArithmeticError("the complete response/source polynomial replay changed")
        if json.loads(certificate_file.read_text()) != certificate:
            raise ArithmeticError("the integrated response certificate or a bound input changed")
        print("Replayed the complete integrated response and every exact source/current coefficient.", flush=True)
    print("Prism: 14 unisolvent sources; 140 integrated outputs; 1204 original coefficients; exact rank 49.")
    print("Pointwise source convergence: 162 automorphism rows; full exact degree-four/five ranks 70 and 126.")
    print("Nonzero contact controls include an unsourced deleted vertex; every double-contact term vanishes.")
    print("Complete result coefficient counts:", [len(p) for p in data["all_source_bare_and_nine_current_polynomials"]])
    print("canonical_certificate_sha256:", digest(certificate), flush=True)


if __name__ == "__main__":
    main()
