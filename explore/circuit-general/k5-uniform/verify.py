"""Replay the guarded K5 recurrence, all new q=9 values, and the missing input relation.

No production graph reducer, vertex-row generator, determinant routine, or
period evaluator is imported.  Existing exact tables are read as stated
lower-level inputs.  Every computed value and source digest is replayed.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import argparse
import json

import numpy as np

import k5_dimension as dimension


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
K5 = HERE.parent / "k5"
CERTIFICATE = HERE / "certificate.json"
NEW_VALUES = HERE / "q9_interior.json"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value):
    return sha256(canonical_json(value)).hexdigest()


def encode_value(value):
    return {label: str(coefficient) for label, coefficient in sorted(value.items()) if coefficient}


@lru_cache(None)
def canon5(weights):
    return dimension.canonical(weights, 5)


@lru_cache(None)
def canon6(weights):
    return dimension.canonical(weights, 6)


@lru_cache(None)
def read_momentum(q):
    contents = json.loads((K5 / f"tables/K5_q{q}.json").read_text())
    require(contents["q"] == q and tuple(map(tuple, contents["edges"])) == dimension.E5,
            "input table conventions changed")
    return {tuple(alpha): {label: F(value) for label, value in coefficients.items()}
            for alpha, kind, coefficients in contents["values"]}


@lru_cache(None)
def read_completed(q):
    return {tuple(beta): {label: F(value) for label, value in coefficients.items()}
            for beta, coefficients in json.loads((K5 / f"data/J_q{q}_tower.json").read_text())}


def multiply(first, second):
    result = defaultdict(F)
    for powers, coefficient in first.items():
        for other, value in second.items():
            result[tuple(a + b for a, b in zip(powers, other))] += coefficient * value
    return {powers: value for powers, value in result.items() if value}


def independent_gram_polynomial():
    """4x4 Gram construction, independent of the helper's 6x6 determinant."""
    zero = (0,) * 10

    def variable(edge, coefficient=F(1)):
        powers = list(zero)
        powers[dimension.I5[tuple(sorted(edge))]] = 1
        return {tuple(powers): coefficient}

    matrix = []
    for i in range(4):
        row = []
        for j in range(4):
            if i == j:
                row.append(variable((i, 4)))
            else:
                entry = {}
                for edge, coefficient in [((i, 4), F(1, 2)), ((j, 4), F(1, 2)),
                                          ((i, j), F(-1, 2))]:
                    entry.update(variable(edge, coefficient))
                row.append(entry)
        matrix.append(row)
    answer = defaultdict(F)
    for permutation in permutations(range(4)):
        inversions = sum(permutation[i] > permutation[j] for i in range(4)
                         for j in range(i + 1, 4))
        term = {zero: F((-1) ** inversions)}
        for i, j in enumerate(permutation):
            term = multiply(term, matrix[i][j])
        for powers, value in term.items():
            answer[powers] += value
    return {powers: value for powers, value in answer.items() if value}


def determinant_checks():
    independent = independent_gram_polynomial()
    require(independent == dict(dimension.gram_terms()), "the two complete determinants differ")
    for edge in range(15):
        require(any(powers[edge] == 0 for powers, value in dimension.cayley_menger_terms()),
                "edge obstruction has no determinant witness")
    for triangle in dimension.TRIANGLES6:
        require(any(sum(powers[i] for i in triangle) == 0
                    for powers, value in dimension.cayley_menger_terms()),
                "triangle obstruction has no determinant witness")
    return {"complete_cayley_menger_terms": len(dimension.cayley_menger_terms()),
            "independent_gram_terms": len(independent),
            "guard_necessity_edge_witnesses": 15,
            "guard_necessity_triangle_witnesses": 20,
            "polynomial_sha256": digest([[list(powers), str(value)]
                                         for powers, value in sorted(independent.items())])}


def bounded_compositions(length, total, cap):
    if length == 0:
        if total == 0:
            yield ()
        return
    for first in range(max(0, total - cap * (length - 1)), min(cap, total) + 1):
        for tail in bounded_compositions(length - 1, total - first, cap):
            yield (first,) + tail


def kind(beta, q):
    if 0 in beta:
        return "T1"
    if any(sum(beta[i] for i in triangle) == q for triangle in dimension.TRIANGLES6):
        return "T2"
    return "G"


def exact_controls():
    controls = []
    for alpha in [(3,) + (5,) * 9, (4, 4) + (5,) * 8,
                  (4, 5, 5, 5, 5, 5, 5, 4, 5, 5)]:
        terms = dimension.lower_momentum(alpha, 8)
        result = dimension.scaled_sum(terms, lambda a: read_momentum(7)[canon5(a)])
        require(result == read_momentum(8)[canon5(alpha)], "independent q=8 table control failed")
        controls.append({"normalization": "I", "q": 8, "weights": list(alpha),
                         "lower_terms": len(terms), "value": encode_value(result)})
    for q, alpha in [(4, (2, 2, 2, 3, 3, 3, 2, 3, 2, 2)), (5, (3,) * 10)]:
        beta = dimension.complete_momentum(alpha, q)
        terms = dimension.lower_completed(beta, q)
        result = dimension.scaled_sum(terms, lambda b: read_completed(q - 1)[canon6(b)])
        normalization = F(factorial(q - 1))
        for a in alpha:
            normalization /= factorial(q - a - 1)
        expected = {label: value * normalization
                    for label, value in read_momentum(q)[canon5(alpha)].items()}
        require(result == expected, "complete-J Gram control failed")
        controls.append({"normalization": "J", "q": q, "weights": list(beta),
                         "lower_terms": len(terms), "value": encode_value(result)})
    return controls


def q9_values():
    labeled = [tuple(6 - a for a in extra) for extra in bounded_compositions(10, 6, 4)]
    require(len(labeled) == comb(15, 9) - 10 * comb(10, 9) == 4905,
            "wrong complete q=9 interior enumeration")
    classes = sorted({canon5(alpha) for alpha in labeled})
    values = []
    for alpha in classes:
        beta = dimension.complete_momentum(alpha, 9)
        terms = dimension.lower_momentum(alpha, 9)
        result = dimension.scaled_sum(terms, lambda a: read_momentum(8)[canon5(a)])
        require(set(result) <= {"1", "z3", "z5"}, "unexpected coefficient label")
        values.append({"alpha": list(alpha), "kind": kind(beta, 9),
                       "coefficients": encode_value(result), "lower_terms": len(terms),
                       "complete_lower_row_sha256": digest([[list(a), str(c)] for a, c in terms])})
    return {"q": 9, "region": "3 <= beta_e <= 7; sum beta_e = 36; alpha=9-beta",
            "labelled_inputs": len(labeled), "orbits": len(classes),
            "kind_counts": dict(sorted(Counter(item["kind"] for item in values).items())),
            "complete_values": values}


def boundary_scan():
    shifts = np.array([powers for powers, value in dimension.gram_terms()], dtype=np.int16)
    results = []
    for q in range(3, 9):
        inputs = np.array(list(read_momentum(q)), dtype=np.int16)
        nonzero = divergent = all_zero = safe = completed_safe = 0
        for start in range(0, len(inputs), 256):
            alpha = inputs[start:start + 256]
            beta = q - alpha
            keep = (beta[:, None, :] > shifts[None, :, :]).all(axis=2)
            lower = alpha[:, None, :] - 1 + shifts[None, :, :]
            converges = (lower > 0).all(axis=2) & (lower < q - 1).all(axis=2)
            for triangle in dimension.TRIANGLES5:
                converges &= lower[:, :, triangle].sum(axis=2) > q - 1
            for clique in dimension.CLIQUES5:
                converges &= lower[:, :, clique].sum(axis=2) > 3 * (q - 1)
            nonzero += int(keep.sum())
            divergent += int((keep & ~converges).sum())
            all_zero += int((~keep.any(axis=1)).sum())
            safe += int((keep.all(axis=1) & converges.all(axis=1)).sum())
            for row in alpha:
                completed = dimension.complete_momentum(tuple(map(int, row)), q)
                completed_safe += (max(completed) <= q - 2 and all(
                    sum(completed[i] for i in triangle) <= 2 * q - 3
                    for triangle in dimension.TRIANGLES6))
        results.append({"q": q, "orbits": len(inputs), "nonzero_formal_terms": nonzero,
                        "divergent_nonzero_formal_terms": divergent,
                        "all_formal_terms_zero_orbits": all_zero,
                        "safe_scalar_recurrences": safe,
                        "safe_completed_recurrences": completed_safe})
    return results


def explicit_failures():
    alpha3 = (1, 1) + (2,) * 8
    require(dimension.momentum_convergent(alpha3, 3), "first failure input changed")
    require(all(any(3 - a <= shift for a, shift in zip(alpha3, powers))
                for powers, value in dimension.gram_terms()), "formal derivative is not zero")
    expected = {"z3": F(5, 3), "z5": F(-5, 3)}
    require(read_momentum(3)[canon5(alpha3)] == expected, "nonzero base control changed")
    try:
        dimension.lower_momentum(alpha3, 3)
    except dimension.BoundaryContributionRequired:
        pass
    else:
        raise AssertionError("a boundary term was silently omitted")

    alpha5 = (2, 2, 3, 4, 2, 3, 3, 4, 3, 4)
    shift = (0, 0, 1, 0, 0, 1, 1, 0, 1, 0)
    lower = tuple(a - 1 + b for a, b in zip(alpha5, shift))
    require(dimension.momentum_convergent(alpha5, 5), "second failure input is not convergent")
    require(dict(dimension.gram_terms())[shift] == F(1, 8), "wrong exact offending coefficient")
    require(all(5 - a > b for a, b in zip(alpha5, shift)), "offending falling factor vanishes")
    require(all(0 < a < 4 for a in lower) and not dimension.momentum_convergent(lower, 4),
            "the positive-weight lower triangle is not divergent")
    require(sum(lower[i] for i in (0, 1, 4)) == 3, "wrong lower triangle sum")
    try:
        dimension.lower_completed(dimension.complete_momentum(alpha5, 5), 5)
    except dimension.DivergentLowerTerm:
        pass
    else:
        raise AssertionError("a divergent lower Gram term was accepted")
    return [{"q": 3, "alpha": list(alpha3), "formal_interior_value": {},
             "actual_value": encode_value(expected), "missing": "coordinate boundary"},
            {"q": 5, "alpha": list(alpha5), "shift": list(shift),
             "gram_coefficient": "1/8", "falling_factor": "1",
             "lower_q": 4, "lower_alpha": list(lower), "triangle_012_sum": 3}]


def x_rows(beta, q, classes):
    index = {item: i for i, item in enumerate(classes)}
    for x, x1, c in permutations(range(6), 3):
        edge = lambda a, b: dimension.I6[tuple(sorted((a, b)))]
        items = [(beta[edge(x, c)] - beta[edge(x, x1)], beta)]
        okay = True
        for i in range(6):
            if i in (x, x1, c) or not beta[edge(x, i)]:
                continue
            coefficient = beta[edge(x, i)]
            for sign, changes in [(-1, [(x, x1, -1), (x, i, 1), (x1, c, 1), (i, c, -1)]),
                                  (1, [(x, i, 1), (x1, i, -1), (x1, c, 1), (x, c, -1)])]:
                shifted = list(beta)
                for first, second, change in changes:
                    shifted[edge(first, second)] += change
                shifted = tuple(shifted)
                if not dimension.completed_convergent(shifted, q):
                    okay = False
                    break
                items.append((sign * coefficient, shifted))
            if not okay:
                break
        if okay:
            row = defaultdict(F)
            for coefficient, weights in items:
                row[index[canon6(weights)]] += coefficient
            yield (x, x1, c), {j: coefficient for j, coefficient in row.items() if coefficient}


def pruefer_trees(vertices):
    """Every labeled tree on five vertices, independently of historical generators."""
    vertices = tuple(vertices)
    for code in product(vertices, repeat=len(vertices) - 2):
        degree = {v: 1 + code.count(v) for v in vertices}
        tree = []
        for vertex in code:
            leaf = min(v for v in vertices if degree[v] == 1)
            tree.append(tuple(sorted((leaf, vertex))))
            degree[leaf] -= 1
            degree[vertex] -= 1
        last = [v for v in vertices if degree[v] == 1]
        tree.append(tuple(last))
        yield tuple(tree)


def missing_input_relation():
    classes = [tuple(map(int, line.split()))
               for line in (K5 / "data/classes_q3.txt").read_text().splitlines()]
    low = tuple(map(int, (K5 / "data/classes_q2.txt").read_text().split()))
    require(len(classes) == 15, "the existing level-3 indexing changed")
    row_specs = [
        (F(22), 10, {13: 1, 9: -1}),
        (F(1, 2), 4, {0: 1, 1: -1}),
        (F(-1, 2), 2, {2: 2, 4: -2, 3: 2, 1: -1}),
        (F(-1), 6, {6: 1, 8: -1, 7: 1, 2: -1}),
        (F(3), 8, {3: 1, 11: -1}),
        (F(3), 8, {8: -2, 10: -2, 4: 1, 12: 1, 11: 1}),
        (F(-1), 4, {7: -1, 12: 1}),
        (F(-1, 2), 7, {6: -2, 8: 2}),
        (F(-10), 13, {13: 1, 10: -2, 14: 2}),
        (F(6), 10, {10: -1, 13: -2, 8: 1, 11: 1}),
        (F(-2), 3, {3: 1, 7: -1, 9: 1, 8: -2, 10: 2}),
    ]
    witness = []
    total = defaultdict(F)
    generated = {}
    for coefficient, source, expected in row_specs:
        if source not in generated:
            generated[source] = list(x_rows(classes[source], 3, classes))
        match = next((triple for triple, row in generated[source] if row == expected), None)
        require(match is not None, "an asserted vertex row was not independently generated")
        for index, value in expected.items():
            total[index] += coefficient * value
        witness.append({"kind": "X", "coefficient": str(coefficient), "source_class": source,
                        "ordered_vertices": list(match),
                        "row": {str(index): str(value) for index, value in sorted(expected.items())}})
    index = {item: i for i, item in enumerate(classes)}
    dimension_row = defaultdict(F)
    dimension_row[15] = F(1)
    trees = list(pruefer_trees(range(1, 6)))
    require(len(trees) == len(set(trees)) == 125, "wrong labeled spanning trees")
    for tree in trees:
        coefficient = F(1, 2)
        shifted = list(low)
        for edge in tree:
            position = dimension.I6[tuple(sorted(edge))]
            coefficient *= low[position]
            shifted[position] += 1
        if not coefficient:
            continue
        for v in range(1, 6):
            position = dimension.I6[(0, v)]
            shifted[position] += 6 - sum(shifted[i] for i in dimension.STARS6[v])
        shifted = tuple(shifted)
        require(dimension.completed_convergent(shifted, 3), "dimension row crosses a boundary")
        dimension_row[index[canon6(shifted)]] -= coefficient
    dimension_row = {j: a for j, a in dimension_row.items() if a}
    for index, value in dimension_row.items():
        total[index] += value
    total = {index: value for index, value in total.items() if value}
    require(total == {15: F(1), 9: F(-24), 14: F(-20)}, "complete row witness failed")
    witness.append({"kind": "DS", "coefficient": "1", "decompletion": 0,
                    "row": {str(index): str(value) for index, value in sorted(dimension_row.items())}})

    # Exact complete graph maps to the two radial wheels.
    for q, expected in [(2, low), (3, classes[9])]:
        beta = [0] * 15
        for u in (2, 3):
            for v in (0, 1, 4, 5):
                beta[dimension.I6[tuple(sorted((u, v)))]] = 1
        for edge in [(0, 4), (1, 4), (1, 5), (0, 5)]:
            beta[dimension.I6[edge]] = q - 1
        beta[dimension.I6[(2, 3)]] = 2 * q - 4
        require(tuple(beta) == expected, "a radial-wheel graph map failed")

    # The two remaining graph maps keep every edge and the finite source frame.
    require(all(classes[5][i] == (2 if (u < 3) != (v < 3) else 0)
                for i, (u, v) in enumerate(dimension.E6)), "the three-bubble graph map failed")
    reduced = list(classes[6])
    require([reduced[dimension.I6[(0, v)]] for v in range(1, 6)] == [0, 0, 2, 2, 2],
            "the unique-star weights changed")
    for edge in combinations((3, 4, 5), 2):
        reduced[dimension.I6[edge]] += 1
    four_edges = tuple(combinations(range(4), 2))
    four_index = {edge: i for i, edge in enumerate(four_edges)}
    obtained = tuple(3 - reduced[dimension.I6[edge]] for edge in combinations((2, 3, 4, 5), 2))
    expected = (2, 1, 1, 1, 2, 2)
    matches = [permutation for permutation in permutations(range(4)) if tuple(
        obtained[four_index[tuple(sorted((permutation[u], permutation[v])))]]
        for u, v in four_edges) == expected]
    require(matches, "the K4 unique-star/decompletion map failed")

    # Check the uniform spectral polynomial identity coefficientwise, not by zeta fitting.
    radial_checks = 0
    for n in range(4, 21):
        for q in range(2, n - 1):
            numerator = {0: F(1)}
            for j in range(1, q - 1):
                changed = defaultdict(F)
                for degree, value in numerator.items():
                    changed[degree + 2] += value
                    changed[degree] -= j * j * value
                numerator = dict(changed)
            lhs = defaultdict(F)
            for degree, value in numerator.items():
                lhs[degree + 2] += value
                lhs[degree] -= (q - 1) ** 2 * value
                lhs[degree] += (q - 1) ** 2 * value
            require({degree: value for degree, value in lhs.items() if value}
                    == {degree + 2: value for degree, value in numerator.items() if value},
                    "radial multiplicity identity failed")
            require(F(comb(2 * n - 2, n - 1), comb(2 * n - 4, n - 2))
                    == F(2 * (2 * n - 3), n - 1), "radial normalization ratio failed")
            radial_checks += 1

    def spectral(s):
        return (-F(5, 12 * (s + 1)) - F(1, 12 * (s + 1) ** 2)
                + F(1, 2 * (s + 1) ** 3) - F(5, 12 * (s - 1))
                + F(1, 12 * (s - 1) ** 2) + F(1, 2 * (s - 1) ** 3)
                + F(5, 6 * s) - F(1, 2 * s ** 3))

    def boundary(s):
        return (-F(5, 12) * (F(1, s - 1) - F(1, s))
                + F(1, 12) * (F(1, (s - 1) ** 2) + F(1, s ** 2))
                + F(1, 2) * (F(1, (s - 1) ** 3) - F(1, s ** 3)))

    principal = {(0, 2): -1, (0, 1): -3, (1, 3): 1, (1, 2): -2, (1, 1): 3}
    # The full residue expression has denominator s^5*(s-1)^4*(s+1)^4.
    # After multiplication by the degree-four multiplicity, a difference
    # numerator has degree at most sixteen. Twenty exact points prove the
    # rational identity; these are polynomial identity checks, not quadrature.
    for s in range(2, 22):
        convolution = sum(((-1) ** (5 + t) * comb(k + t - 2, k - 1)
                           * F(first * second, (s + l - j) ** (k + t - 1))
                           for (j, k), first in principal.items()
                           for (l, t), second in principal.items()), F())
        require(convolution * F(s * s * (s * s - 1), 12) == spectral(s),
                "the complete K4 convolution residue expression failed")
        require(spectral(s) - F(1, 2 * s ** 3) == boundary(s) - boundary(s + 1),
                "complete rational telescoping identity failed")
    require(boundary(2) == F(1, 3), "wrong retained finite boundary")

    # Existing values are a separate control; they were not used to obtain the row witness.
    actual = defaultdict(F)
    for label, value in read_completed(2)[low].items():
        actual[label] += value
    for j, coefficient in [(14, -4), (6, -16), (5, F(-16, 3))]:
        for label, value in read_completed(3)[classes[j]].items():
            actual[label] += coefficient * value
    require(not any(actual.values()), "the independently derived input relation fails the stored values")
    return {"level_3_classes": [list(beta) for beta in classes], "row_witness": witness,
            "proved_intermediate": {"J2": "1", "J3_9": "-24", "J3_14": "-20"},
            "proved_extra_input_relation": {"J2": "1", "J3_14": "-4", "J3_6": "-16",
                                             "J3_5": "-16/3"},
            "class_6_k4_relabeling": list(matches[0]),
            "radial_coefficient_controls": radial_checks, "telescoping_rational_controls": 20,
            "exact_telescoping_boundary_at_2": "1/3"}


def build():
    determinant = determinant_checks()
    controls = exact_controls()
    values = q9_values()
    scan = boundary_scan()
    failures = explicit_failures()
    relation = missing_input_relation()
    inputs = [Path(__file__), HERE / "k5_dimension.py", HERE / "DERIVATION.md", HERE / "README.md"]
    inputs += [K5 / f"tables/K5_q{q}.json" for q in range(3, 9)]
    inputs += [K5 / f"data/J_q{q}_tower.json" for q in range(2, 5)]
    inputs += [K5 / f"data/classes_q{q}.txt" for q in (2, 3)]
    inputs += [ROOT / "notes/methods-paper/sections/periods.tex"]
    bindings = {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    certificate = {"scope": "guarded dimension lowering and origin of the extra input relation; all-q K5 closure remains open",
                   "determinant": determinant, "independent_table_controls": controls,
                   "new_level": {key: value for key, value in values.items() if key != "complete_values"},
                   "q9_values_sha256": digest(values), "complete_boundary_scan": scan,
                   "explicit_failures": failures, "extra_input_relation": relation,
                   "source_sha256": bindings}
    return values, certificate


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="explicitly refresh both complete certificates")
    arguments = parser.parse_args()
    values, certificate = build()
    if arguments.write:
        write_json(NEW_VALUES, values)
        write_json(CERTIFICATE, certificate)
        print("Wrote every exact q=9 value and the complete guarded K5 certificate.")
    else:
        require(values == json.loads(NEW_VALUES.read_text()), "a complete q=9 value or row digest changed")
        require(certificate == json.loads(CERTIFICATE.read_text()), "certificate or bound input digest changed")
        print("Replayed every exact value, guarded relation, failure, and source digest.")
    print("q=9: 4905 labeled inputs; 73 orbits; 13 generic orbits; 9490 complete lower terms.")
    print("Independent existing-table controls: 5; complete boundary scan: 54391 K5 orbits.")
    print("Extra input relation: 11 vertex rows + 1 dimension row, radial identity, retained boundary 1/3.")
    print("All-q generic K5 numerical closure remains open.")
    print("canonical_certificate_sha256:", digest(certificate))


if __name__ == "__main__":
    main()
