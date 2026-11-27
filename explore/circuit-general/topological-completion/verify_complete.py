"""Independent exact checks of the new half-integer tables and full kernels."""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import gzip
import json
from pathlib import Path
import sys
import time

from flint import fmpq, fmpq_mpoly_ctx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "code"), str(ROOT / "explore/circuit-general/cores")]
from htk.graph import Graph
from htk.topological_polynomial import BulkPolynomial, rational
from htk.bulk_kernels import evaluate_bulk
import cores
import gutil
from source_covariance import automorphisms


def read(path):
    return json.loads(path.read_text())


def scalar_table(name, pole):
    filename = name + "_q" + cores.qname(pole) + ".json"
    choices = (HERE / "tables" / filename, ROOT / "explore/circuit-general/cores/tables" / filename)
    path = next(p for p in choices if p.exists())
    return {tuple(a): {k: Fraction(v) for k, v in value.items()} for a, value in read(path)["values"]}


def modular_rank(rows, columns, prime=2147483647):
    pivots = {}
    witnesses = []
    for index, (coefficients, _) in sorted(enumerate(rows), key=lambda item: len(item[1][0])):
        row = {}
        for k, v in coefficients:
            v = Fraction(v)
            if v.denominator % prime == 0:
                raise ArithmeticError("rank certificate prime divides an input denominator")
            row[int(k)] = v.numerator*pow(v.denominator, -1, prime) % prime
        row = {k: v for k, v in row.items() if v}
        while row:
            column = min(row)
            if column not in pivots:
                inverse = pow(row[column], -1, prime)
                pivots[column] = {k: v*inverse % prime for k, v in row.items()}
                witnesses.append(index)
                break
            factor = row[column]
            for k, v in pivots[column].items():
                value = (row.get(k, 0)-factor*v) % prime
                if value:
                    row[k] = value
                else:
                    row.pop(k, None)
        if len(pivots) == columns:
            break
    assert len(pivots) == columns, (len(pivots), columns)
    return {"prime": prime, "rank": len(pivots), "independent_original_rows": witnesses}


def verify_table(path):
    data = read(path)
    if "rows" in data:
        rows = data["rows"]
    else:
        packed = (path.parent / data["relation_file"]).read_bytes()
        raw = gzip.decompress(packed)
        assert sha256(packed).hexdigest() == data["relation_compressed_sha256"]
        assert sha256(raw).hexdigest() == data["relation_uncompressed_sha256"]
        rows = json.loads(raw)
    values = [{k: Fraction(v) for k, v in val.items()} for _, val in data["values"]]
    for coefficients, rhs in rows:
        answer = {}
        for column, scalar in coefficients:
            for label, value in values[column].items():
                answer[label] = answer.get(label, 0) + Fraction(scalar)*value
        assert {k: v for k, v in answer.items() if v} == {k: Fraction(v) for k, v in rhs.items()}
    rank = modular_rank(rows, len(values))
    return {"table": path.name, "orbits": len(values), "exact_original_rows": len(rows),
            "rank_certificate": rank}


def dimension_shift(name, pole):
    lower, upper = scalar_table(name, pole), scalar_table(name, pole+1)
    edges = cores.FAMILIES[name]
    autos = gutil.edge_automorphisms(edges)
    cotrees = gutil.cotrees(edges)
    for alpha, expected in lower.items():
        actual = {}
        for cotree in cotrees:
            shifted = tuple(a+(e in cotree) for e, a in enumerate(alpha))
            representative = gutil.orbit_rep(shifted, autos)
            for label, scalar in upper[representative].items():
                actual[label] = actual.get(label, 0) + scalar
        assert {k: v for k, v in actual.items() if v} == expected
    return {"core": name, "lower_pole": str(pole), "independent_identities": len(lower),
            "summands_per_identity": len(cotrees)}


def verify_polynomial_identities(name):
    data = read(HERE / (name + "_complete_polynomials.json"))
    ctx = fmpq_mpoly_ctx.get(tuple("source%d" % i for i in range(5)))
    source = list(ctx.gens())
    full = source + [-sum(source)]
    values = []
    for output, terms in enumerate(data["bare_and_all_nine_current_polynomials"]):
        by_label = {}
        for exponent, label, scalar in terms:
            assert sum(exponent) == (4 if output == 0 else 5)
            assert tuple(exponent) not in by_label.setdefault(label, {})
            by_label[label][tuple(exponent)] = Fraction(scalar)
        values.append({label: ctx.from_dict({a: fmpq(c.numerator, c.denominator)
                                             for a, c in coefficients.items()})
                       for label, coefficients in by_label.items()})
    labels = sorted({label for output in values for label in output})
    graph = Graph(6, data["edges"])
    incidence = graph.incidence()
    for vertex in range(5):
        for label in labels:
            assert sum(incidence[edge][vertex]*values[1+edge].get(label, 0) for edge in range(9)) == source[vertex]*values[0].get(label, 0)
    actions = automorphisms(tuple(map(tuple, graph.edges)))
    for permutation, edge_map, orientation, character in actions:
        arguments = [full[permutation.index(i)] for i in range(5)]
        for label in labels:
            bare = values[0].get(label, ctx.constant(0))
            assert bare.compose(*arguments) == character*bare
            for old, new in enumerate(edge_map):
                transformed = values[1+new].get(label, ctx.constant(0)).compose(*arguments)
                assert transformed == character*orientation[old]*values[1+old].get(label, 0)
    return {"core": name, "homogeneous_outputs": 10, "all_graph_automorphisms": len(actions),
            "independent_source_conservation_identities": 5*len(labels), "constant_labels": labels}


def holdout(name):
    """New full-source values from the independent Gaussian-minor algorithm."""
    graph = Graph(6, cores.FAMILIES[name])
    sources = [[1], [2], [3], [4], [5]]
    density = BulkPolynomial(graph, sources, H=1, T=1, method="gaussian_minor")
    density.check_at((2, 3, 5, 7, 11, 13, 17, 19, 23))
    autos = gutil.edge_automorphisms(tuple(graph.edges))
    orbit_cache = {}
    outputs = []
    for word in [()] + [((edge, 0),) for edge in range(9)]:
        numerator, pole = density.inserted(word)
        coefficients = {}
        for exponent, scalar in numerator.to_dict().items():
            alpha = tuple(int(x)+1 for x in exponent)
            if alpha not in orbit_cache:
                orbit_cache[alpha] = gutil.orbit_rep(alpha, autos)
            representative = orbit_cache[alpha]
            coefficients[representative] = coefficients.get(representative, Fraction(0)) + rational(scalar)
        table = scalar_table(name, pole)
        actual = {}
        for alpha, coefficient in coefficients.items():
            if coefficient:
                for label, scalar in table[alpha].items():
                    actual[label] = actual.get(label, 0) + coefficient*scalar
        actual = {k: v for k, v in actual.items() if v}
        expected = evaluate_bulk(graph, sources, H=1, T=1, word=word)
        assert actual == dict(expected.coefficients)
        outputs.append({"word": word, "value": {k: str(v) for k, v in actual.items()},
                        "complete_scalar_orbits": sum(bool(v) for v in coefficients.values())})
    return {"core": name, "independent_density_method": "Gaussian linear minor",
            "sources": sources, "all_outputs": outputs}


def known_prism_controls():
    graph = Graph(6, [(0,1),(0,2),(0,3),(1,4),(1,5),(2,3),(2,5),(3,4),(4,5)])
    cases = [([1,2,-1,3,-2,-3], {"sqrtpi^4": Fraction(-808,63), "sqrtpi^8": Fraction(2041,1680)}),
             ([1,0,-1,0,0,0], {}),
             ([1,1,0,-1,0,-1], {"sqrtpi^4": Fraction(-33,70), "sqrtpi^8": Fraction(17,360)})]
    for source, expected in cases:
        actual = evaluate_bulk(graph, [[x] for x in source], H=1, T=1)
        assert dict(actual.coefficients) == expected
    return {"historical_prism_values": len(cases), "different_edge_order_and_vertex_map": True}


def main():
    start = time.monotonic()
    results = {"tables": [], "dimension_shifts": [], "polynomial_identities": [], "holdouts": []}
    for path in sorted((HERE / "tables").glob("*.json")):
        result = verify_table(path)
        results["tables"].append(result)
        print(path.name, result["orbits"], "full rank; exact rows", result["exact_original_rows"], flush=True)
    for name, pole in (("K33", Fraction(7,2)), ("K33", Fraction(9,2)), ("prism", Fraction(9,2))):
        result = dimension_shift(name, pole)
        results["dimension_shifts"].append(result)
        print("independent dimension shift", result, flush=True)
    for name in ("prism", "K33"):
        results["polynomial_identities"].append(verify_polynomial_identities(name))
        results["holdouts"].append(holdout(name))
        print(name, "all source identities and all ten independent Gaussian-minor holdouts pass", flush=True)
    results["historical_controls"] = known_prism_controls()
    results["seconds"] = round(time.monotonic()-start, 2)
    results["passed"] = True
    (HERE / "verification.json").write_text(json.dumps(results, separators=(",", ":")) + "\n")
    print("ALL EXACT CHECKS PASS", results["seconds"], "seconds", flush=True)


if __name__ == "__main__":
    main()
