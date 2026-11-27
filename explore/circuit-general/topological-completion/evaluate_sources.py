"""All-source four-loop bare kernels and all original first currents."""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "code"), str(ROOT / "explore/circuit-general/cores")]
from htk.graph import Graph
from htk.topological_polynomial import BulkPolynomial, rational
import cores
import gutil
from source_covariance import Elimination, evaluate


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    temporary.replace(path)


def polynomial_record(p):
    rows = [[list(map(int, e)), str(rational(c))] for e, c in sorted(p.to_dict().items())]
    return {"terms": len(rows), "sha256": sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()}


def generate(name):
    covariance = json.loads((HERE / (name + "_source_covariance.json")).read_text())
    graph = Graph(6, covariance["edges"])
    autos = gutil.edge_automorphisms(tuple(graph.edges))
    @lru_cache(maxsize=None)
    def orbit(alpha):
        return gutil.orbit_rep(alpha, autos)
    start = time.monotonic()
    for index, sample in enumerate(covariance["samples"]):
        output = HERE / "samples" / (name + "_%02d.json" % index)
        if output.exists():
            continue
        sources = [[int(x)] for x in sample["sources"]]
        density = BulkPolynomial(graph, sources, H=1, T=1, method="pfaffian")
        density.check_at(tuple(range(1, 10)))
        density.check_at((2, 5, 3, 7, 11, 13, 17, 19, 23))
        row = {"schema": "complete-topological-current-sample-v1", "core": name,
               "sample": index, "sources": sources, "density": polynomial_record(density.numerator),
               "U": polynomial_record(density.U), "pole": str(density.power),
               "all_original_current_numerators": [polynomial_record(p[0]) for p in density.current_numerators],
               "complete_polynomial_circuit_identities": True,
               "independent_grounded_laplacian_controls": 2, "targets": []}
        for word in [()] + [((edge, 0),) for edge in range(9)]:
            numerator, power = density.inserted(word)
            coefficients = {}
            for exponent, coefficient in numerator.to_dict().items():
                representative = orbit(tuple(int(x)+1 for x in exponent))
                coefficients[representative] = coefficients.get(representative, Fraction(0)) + rational(coefficient)
            row["targets"].append({"word": word, "pole": str(power),
                                   "complete_numerator": polynomial_record(numerator),
                                   "scalar_orbits": [[list(a), str(c)] for a, c in sorted(coefficients.items()) if c]})
        write(output, row)
        print(name, index, "complete", "sources", sample["sources"], "density terms", len(density.numerator),
              "all nine currents", "seconds", round(time.monotonic()-start, 2), flush=True)


@lru_cache(maxsize=None)
def table(name, pole):
    filename = name + "_q" + cores.qname(Fraction(pole)) + ".json"
    paths = (HERE / "tables" / filename, ROOT / "explore/circuit-general/cores/tables" / filename)
    path = next((p for p in paths if p.exists()), None)
    if path is None:
        raise FileNotFoundError("required complete period table: " + filename)
    data = json.loads(path.read_text())
    return {tuple(a): {label: Fraction(value) for label, value in val.items()} for a, val in data["values"]}


def solve(name):
    covariance_path = HERE / (name + "_source_covariance.json")
    covariance = json.loads(covariance_path.read_text())
    basis = [{int(k): Fraction(v) for k, v in row} for row in covariance["basis"]]
    rows = []
    answers = []
    sample_records = []
    for index, sample in enumerate(covariance["samples"]):
        path = HERE / "samples" / (name + "_%02d.json" % index)
        data = json.loads(path.read_text())
        all_values = []
        for target in data["targets"]:
            scalar = table(name, target["pole"]) if target["scalar_orbits"] else {}
            value = {}
            for alpha, coefficient in target["scalar_orbits"]:
                for label, entry in scalar[tuple(alpha)].items():
                    value[label] = value.get(label, Fraction(0)) + Fraction(coefficient)*entry
            value = {k: v for k, v in value.items() if v}
            all_values.append(value)
        for output in range(10):
            rows.append({j: Fraction(v[output]) for j, v in enumerate(sample["evaluations"]) if Fraction(v[output])})
            answers.append(all_values[output])
        sample_records.append({"sources": sample["sources"],
                               "bare_and_all_nine_currents": [{k: str(v) for k, v in val.items()} for val in all_values]})
    # Inhomogeneous rational elimination with a vector of exact named
    # constants as its right hand side.  Every surplus sample is checked.
    pivots = {}
    for row, answer in zip(rows, answers):
        row, answer = dict(row), dict(answer)
        while row:
            column = min(row)
            if column not in pivots:
                value = row[column]
                pivots[column] = ({k: v/value for k, v in row.items()}, {k: v/value for k, v in answer.items()})
                break
            factor = row[column]
            old, result = pivots[column]
            for k, v in old.items():
                row[k] = row.get(k, 0) - factor*v
                if not row[k]:
                    del row[k]
            for k, v in result.items():
                answer[k] = answer.get(k, 0) - factor*v
                if not answer[k]:
                    del answer[k]
        else:
            if answer:
                raise ArithmeticError("complete source covariance disagrees with an exact integral")
    if len(pivots) != len(basis):
        raise ArithmeticError("the complete source polynomial is underdetermined")
    coefficients = [{} for _ in basis]
    for column, (row, answer) in sorted(pivots.items(), reverse=True):
        value = dict(answer)
        for k, factor in row.items():
            if k == column:
                continue
            for label, scalar in coefficients[k].items():
                value[label] = value.get(label, 0) - factor*scalar
        coefficients[column] = {k: v for k, v in value.items() if v}
    full = {}
    for vector, value in zip(basis, coefficients):
        for column, factor in vector.items():
            for label, scalar in value.items():
                key = (column, label)
                full[key] = full.get(key, 0) + factor*scalar
    full = {k: v for k, v in full.items() if v}
    size4, size5 = len(covariance["monomials"]["4"]), len(covariance["monomials"]["5"])
    outputs = [[] for _ in range(10)]
    for (column, label), coefficient in sorted(full.items()):
        if column < size4:
            output, exponent = 0, covariance["monomials"]["4"][column]
        else:
            output = 1 + (column-size4)//size5
            exponent = covariance["monomials"]["5"][(column-size4)%size5]
        outputs[output].append([exponent, label, str(coefficient)])
    data = {"schema": "complete-four-loop-topological-polynomials-v1", "core": name,
            "H": 1, "T": 1, "edges": covariance["edges"], "source_vertices": list(range(5)),
            "root_source": "minus the sum of all five original sources",
            "bare_and_all_nine_current_polynomials": outputs,
            "source_covariance_dimension": len(basis), "exact_fit_rank": len(pivots),
            "all_evaluated_samples": sample_records,
            "covariance_certificate_sha256": sha256(covariance_path.read_bytes()).hexdigest()}
    write(HERE / (name + "_complete_polynomials.json"), data)
    print(name, "FULL source polynomials", [len(x) for x in outputs], "exact fit rank", len(pivots), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("generate", "solve"))
    parser.add_argument("core", choices=("prism", "K33"))
    args = parser.parse_args()
    (generate if args.action == "generate" else solve)(args.core)
