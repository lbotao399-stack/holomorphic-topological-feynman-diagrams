"""Independent exact controls for the fixed-core stabilization theorem.

Only the Python standard library is used.  No production graph, period,
sector, or integration module is imported.  The unbounded theorem uses
algebraic de Rham comparison; these checks audit its finite algebraic
identities and boundary examples, not a numerical extrapolation of ranks.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial, prod
from pathlib import Path
import argparse
import json


HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "checks.json"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value):
    return sha256(canonical(value)).hexdigest()


def clean(poly):
    require(all(isinstance(value, (int, Fraction)) for value in poly.values()),
            "a polynomial coefficient is not exact rational arithmetic")
    return {tuple(power): Fraction(value) for power, value in poly.items() if value}


def add(*polys):
    result = {}
    for poly in polys:
        for power, value in poly.items():
            result[power] = result.get(power, Fraction()) + value
    return clean(result)


def scale(poly, value):
    return clean({power: value * coefficient for power, coefficient in poly.items()})


def multiply(left, right):
    result = {}
    for first, coefficient in left.items():
        for second, other in right.items():
            power = tuple(a + b for a, b in zip(first, second))
            result[power] = result.get(power, Fraction()) + coefficient * other
    return clean(result)


def shift(poly, powers):
    return clean({tuple(a + b for a, b in zip(power, powers)): value
                  for power, value in poly.items()})


def derivative(poly, coordinate):
    result = {}
    for power, value in poly.items():
        if power[coordinate]:
            changed = list(power)
            changed[coordinate] -= 1
            result[tuple(changed)] = value * power[coordinate]
    return clean(result)


def evaluate(poly, point):
    return sum((value * prod(x ** exponent for x, exponent in zip(point, power))
                for power, value in poly.items()), Fraction())


def encode_poly(poly):
    return [{"powers": list(power), "coefficient": str(value)}
            for power, value in sorted(poly.items())]


def determinant(matrix, variables):
    """Leibniz determinant of a matrix of sparse exact polynomials."""
    size = len(matrix)
    result = {}
    for order in permutations(range(size)):
        inversions = sum(order[i] > order[j]
                         for i in range(size) for j in range(i + 1, size))
        term = {(0,) * variables: Fraction((-1) ** inversions)}
        for row, column in enumerate(order):
            term = multiply(term, matrix[row][column])
        result = add(result, term)
    return result


def matrix_rank(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    if not matrix:
        return 0
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((row for row in range(rank, len(matrix))
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        divisor = matrix[rank][column]
        matrix[rank] = [value / divisor for value in matrix[rank]]
        for row in range(rank + 1, len(matrix)):
            multiple = matrix[row][column]
            matrix[row] = [a - multiple * b
                           for a, b in zip(matrix[row], matrix[rank])]
        rank += 1
    return rank


def incidence_rank(vertices, edges, subset):
    matrix = [[0] * len(subset) for _ in range(vertices)]
    for column, edge in enumerate(subset):
        left, right = edges[edge]
        matrix[left][column] += 1
        matrix[right][column] -= 1
    return matrix_rank(matrix)


def greedy_tree(vertices, edges, order):
    parent = list(range(vertices))

    def root(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex

    chosen = []
    for edge in order:
        left, right = map(root, edges[edge])
        if left != right:
            parent[left] = right
            chosen.append(edge)
    return tuple(sorted(chosen))


def spanning_trees(vertices, edges):
    return tuple(subset for subset in combinations(range(len(edges)), vertices - 1)
                 if len(greedy_tree(vertices, edges, subset)) == vertices - 1)


def kirchhoff_from_determinant(vertices, edges):
    """Compute prod(t) det L(1/t) through conductance monomials."""
    count = len(edges)
    laplacian = [[{} for _ in range(vertices - 1)] for _ in range(vertices - 1)]
    for edge, (left, right) in enumerate(edges):
        if left == right:
            continue
        powers = tuple(int(index == edge) for index in range(count))
        variable = {powers: Fraction(1)}
        for first, second, sign in ((left, left, 1), (right, right, 1),
                                    (left, right, -1), (right, left, -1)):
            if first < vertices - 1 and second < vertices - 1:
                laplacian[first][second] = add(
                    laplacian[first][second], scale(variable, sign))
    tree_polynomial = determinant(laplacian, count)
    require(all(all(value in (0, 1) for value in powers)
                for powers in tree_polynomial), "tree determinant is not square-free")
    return {tuple(1 - value for value in powers): coefficient
            for powers, coefficient in tree_polynomial.items()}


def sector_image(poly, order):
    result = {}
    for powers, value in poly.items():
        image = tuple(sum(powers[edge] for edge in order[:length])
                      for length in range(1, len(order)))
        result[image] = result.get(image, Fraction()) + value
    return clean(result)


def sector_numerator(poly, pole, order, prefix_loops):
    return shift(sector_image(poly, order),
                 tuple(index - pole * loops
                       for index, loops in enumerate(prefix_loops)))


def compositions(total, length):
    if total < length:
        return
    if length == 1:
        yield (total,)
    else:
        for first in range(1, total - length + 2):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail


GRAPHS = (
    ("two_edge_cycle", 2, ((0, 1), (0, 1)), 5),
    ("triangle", 3, ((0, 1), (0, 2), (1, 2)), 5),
    ("square", 4, ((0, 1), (1, 2), (2, 3), (0, 3)), 6),
    ("diamond", 4, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3)), 4),
    ("complete_four", 4, tuple(combinations(range(4), 2)), 4),
    ("triple_parallel", 2, ((0, 1), (0, 1), (0, 1)), 4),
    ("triangle_with_bridge", 4, ((0, 1), (0, 2), (1, 2), (2, 3)), 5),
    ("parallel_with_self_loop", 2, ((0, 1), (0, 1), (0, 0)), 4),
)


def graph_control(name, vertices, edges, last_pole):
    count = len(edges)
    trees = spanning_trees(vertices, edges)
    require(trees, "graph must be connected")
    cotrees = tuple(tuple(edge for edge in range(count) if edge not in tree)
                    for tree in trees)
    from_trees = {tuple(int(edge in cotree) for edge in range(count)): Fraction(1)
                  for cotree in cotrees}
    from_determinant = kirchhoff_from_determinant(vertices, edges)
    require(from_trees == from_determinant, "independent Kirchhoff constructions disagree")

    loops = {}
    for mask in range(1 << count):
        subset = tuple(edge for edge in range(count) if mask & (1 << edge))
        loops[mask] = len(subset) - incidence_rank(vertices, edges, subset)
    require(loops[(1 << count) - 1] == count - vertices + 1, "wrong full loop number")

    sector_trace = []
    flags = []
    cotree_checks = 0
    for order in permutations(range(count)):
        prefixes = tuple(sum(1 << edge for edge in order[:length])
                         for length in range(1, count))
        prefix_loops = tuple(loops[mask] for mask in prefixes)
        greed = greedy_tree(vertices, edges, order)
        require(greed in trees, "greedy result is not a spanning tree")
        greedy_cotree = set(range(count)) - set(greed)
        require(all(len(greedy_cotree.intersection(order[:length])) == prefix_loops[length - 1]
                    for length in range(1, count)), "greedy prefix equality failed")
        normalized = shift(sector_image(from_determinant, order),
                           tuple(-value for value in prefix_loops))
        require(all(all(value >= 0 for value in power) for power in normalized),
                "a cotree violates greedy monomial divisibility")
        require(normalized.get((0,) * (count - 1)) == 1, "constant term is not exactly one")
        require(all(value > 0 for value in normalized.values()), "negative normalized coefficient")
        require(sum(normalized.values()) == len(trees), "wrong value at the unit corner")
        require(shift(normalized, prefix_loops) == sector_image(from_determinant, order),
                "exact sector factorization failed")
        cotree_checks += len(cotrees)
        sector_trace.append({"order": order, "prefix_loops": prefix_loops,
                             "greedy_tree": greed, "polynomial": encode_poly(normalized)})
        flags.append((order, prefixes, prefix_loops))

    proper_subsets = tuple((mask, tuple(edge for edge in range(count) if mask & (1 << edge)))
                           for mask in range(1, (1 << count) - 1))

    def converges(weights, pole):
        return all(sum(weights[edge] for edge in subset) > pole * loops[mask]
                   for mask, subset in proper_subsets)

    level_records = []
    convergence_checks = 0
    nested_terms = 0
    coordinate_checks = 0
    rejected_witnesses = 0
    for pole in range(1, last_pole + 1):
        accepted = []
        tested = 0
        for weights in compositions(pole * (count - vertices + 1), count):
            tested += 1
            convergence_checks += 1
            if not converges(weights, pole):
                mask, subset = next((mask, subset) for mask, subset in proper_subsets
                                    if sum(weights[edge] for edge in subset) <= pole * loops[mask])
                order = subset + tuple(edge for edge in range(count) if edge not in subset)
                prefix_loops = tuple(loops[sum(1 << edge for edge in order[:length])]
                                     for length in range(1, count))
                numerator = {tuple(value - 1 for value in weights): Fraction(1)}
                pulled = sector_numerator(numerator, pole, order, prefix_loops)
                require(next(iter(pulled))[len(subset) - 1] <= -1,
                        "rejected weight has no negative exponent at its witness")
                rejected_witnesses += 1
                continue
            accepted.append(list(weights))
            for order, prefixes, prefix_loops in flags:
                numerator = {tuple(value - 1 for value in weights): Fraction(1)}
                pulled = sector_numerator(numerator, pole, order, prefix_loops)
                powers = next(iter(pulled))
                expected = tuple(sum(weights[edge] for edge in order[:length])
                                 - pole * prefix_loops[length - 1] - 1
                                 for length in range(1, count))
                require(powers == expected and all(value >= 0 for value in powers),
                        "convergent weight has a coordinate pole")
                coordinate_checks += count - 1
            numerator = {tuple(value - 1 for value in weights): Fraction(1)}
            raised = {}
            for cotree in cotrees:
                changed = tuple(value + int(edge in cotree)
                                for edge, value in enumerate(weights))
                require(sum(changed) == (pole + 1) * (count - vertices + 1),
                        "pole shift changed the homogeneous degree")
                require(converges(changed, pole + 1), "pole shift lost absolute convergence")
                raised = add(raised, {tuple(value - 1 for value in changed): Fraction(1)})
                nested_terms += 1
            require(raised == multiply(numerator, from_determinant),
                    "the complete cotree shift is not multiplication by U")
        level_records.append({"pole": pole, "tested_weights": tested,
                              "convergent_weights": accepted})
    return {
        "name": name, "vertices": vertices, "edges": [list(edge) for edge in edges],
        "loop_number": count - vertices + 1,
        "spanning_trees": [list(tree) for tree in trees],
        "kirchhoff_polynomial": encode_poly(from_determinant),
        "sector_count": len(sector_trace), "cotree_factorization_checks": cotree_checks,
        "all_sector_records_sha256": digest(sector_trace),
        "first_sector": sector_trace[0],
        "last_checked_pole": last_pole, "levels": level_records,
        "tested_weight_count": convergence_checks,
        "rejected_weight_witness_count": rejected_witnesses,
        "convergent_sector_coordinate_checks": coordinate_checks,
        "convergence_preserving_shift_terms": nested_terms,
    }


def jacobian_controls():
    records = []
    for variables in range(1, 6):
        parameters = [{tuple(int(column >= row) for column in range(variables)): Fraction(1)}
                      for row in range(variables)]
        jacobian = [[derivative(poly, column) for column in range(variables)]
                    for poly in parameters]
        calculated = determinant(jacobian, variables)
        expected = {tuple(range(variables)): Fraction(1)}
        require(calculated == expected, "ratio-coordinate Jacobian failed")
        records.append({"variables": variables, "jacobian": encode_poly(calculated)})
    return records


def signed_numerator_controls():
    _, vertices, edges, _ = next(row for row in GRAPHS if row[0] == "complete_four")
    count = len(edges)
    trees = spanning_trees(vertices, edges)
    cotrees = [tuple(edge for edge in range(count) if edge not in tree) for tree in trees]
    numerator = {tuple(int(edge in cotree) for edge in range(count)): coefficient
                 for cotree, coefficient in zip(cotrees[:3], (Fraction(1), Fraction(-2), Fraction(1, 3)))}
    trace = []
    for order in permutations(range(count)):
        prefix_loops = tuple(length - incidence_rank(vertices, edges, order[:length])
                             for length in range(1, count))
        pulled = sector_numerator(numerator, 3, order, prefix_loops)
        require(len(pulled) == len(numerator), "distinct homogeneous monomials were identified")
        require(all(all(value >= 0 for value in powers) for powers in pulled),
                "a convergent complete signed numerator has a coordinate pole")
        trace.append({"order": order, "numerator": encode_poly(pulled)})

    divergent = {(3, 0, 0, 0, 0, 0): Fraction(1),
                 (0, 3, 0, 0, 0, 0): Fraction(-1)}
    order = (1, 2, 3, 4, 5, 0)
    prefix_loops = tuple(length - incidence_rank(vertices, edges, order[:length])
                         for length in range(1, count))
    pulled = sector_numerator(divergent, 3, order, prefix_loops)
    least = min(powers[-1] for powers in pulled)
    require(least == -2, "wrong divergent signed-numerator exponent")
    leading = {powers[:-1]: value for powers, value in pulled.items() if powers[-1] == least}
    require(leading == {(0, 1, 2, 0): Fraction(1)}, "wrong nonzero leading boundary coefficient")
    return {
        "graph": "complete_four", "pole": 3,
        "convergent_signed_numerator": encode_poly(numerator),
        "checked_sector_count": len(trace), "all_pullbacks_sha256": digest(trace),
        "first_pullback": trace[0],
        "divergent_signed_numerator": encode_poly(divergent),
        "divergent_witness_order": list(order),
        "divergent_pullback": encode_poly(pulled),
        "least_last_coordinate_exponent": least,
        "leading_boundary_coefficient": encode_poly(leading),
    }


def beta_by_polynomial_integral(first, second):
    return sum((Fraction((-1) ** power * comb(second - 1, power), first + power)
                for power in range(second)), Fraction())


def simplex_by_beta_recursion(weights):
    return prod(beta_by_polynomial_integral(weight, sum(weights[index + 1:]))
                for index, weight in enumerate(weights[:-1]))


def cycle_controls():
    records = []
    value_count = 0
    nesting_count = 0
    trace = []
    for edges in range(2, 9):
        require(all(not list(compositions(pole, edges)) for pole in range(1, edges)),
                "a cycle has positive integer weights before its first level")
        levels = []
        for pole in range(edges, edges + 4):
            values = []
            for weights in compositions(pole, edges):
                calculated = simplex_by_beta_recursion(weights)
                expected = Fraction(prod(factorial(value - 1) for value in weights), factorial(pole - 1))
                require(calculated == expected, "simplex moment disagrees with factorial formula")
                raised = Fraction()
                for edge in range(edges):
                    changed = tuple(value + int(index == edge) for index, value in enumerate(weights))
                    raised += simplex_by_beta_recursion(changed)
                    nesting_count += 1
                require(raised == calculated, "cycle nesting failed under independent integration")
                values.append({"weights": list(weights), "value": str(calculated)})
                value_count += 1
            levels.append({"pole": pole, "checked_weight_count": len(values)})
            trace.append({"edges": edges, "pole": pole, "values": values})
        first = simplex_by_beta_recursion((1,) * edges)
        require(first == Fraction(1, factorial(edges - 1)), "wrong first nonzero cycle period")
        records.append({"edge_count": edges, "zero_levels": list(range(1, edges)),
                        "first_nonzero_level": edges, "first_period": str(first),
                        "checked_levels": levels})
    return {"cycles": records, "integrated_weight_count": value_count,
            "integrated_shift_term_count": nesting_count, "all_integrals_sha256": digest(trace)}


def divide_univariate_polynomial(numerator, denominator):
    quotient = {}
    remainder = dict(numerator)
    denominator_degree = max(power[0] for power in denominator)
    while remainder and max(power[0] for power in remainder) >= denominator_degree:
        degree = max(power[0] for power in remainder)
        coefficient = remainder[(degree,)] / denominator[(denominator_degree,)]
        term = {(degree - denominator_degree,): coefficient}
        quotient = add(quotient, term)
        remainder = add(remainder, scale(multiply(term, denominator), -1))
    return quotient, remainder


def univariate_controls():
    """Use u=1+x, so every rational function here is a Laurent polynomial."""
    trace = []
    selected = []
    endpoint_factor = {(2,): Fraction(-1), (1,): Fraction(3), (0,): Fraction(-2)}
    for power in range(13):
        for pole in range(13):
            integrand = {(index - pole,): Fraction(comb(power, index) * (-1) ** (power - index))
                         for index in range(power + 1)}
            log_coefficient = integrand.get((-1,), Fraction())
            primitive = {(exponent + 1,): value / (exponent + 1)
                         for (exponent,), value in integrand.items() if exponent != -1}
            constant = evaluate(primitive, (Fraction(2),)) - evaluate(primitive, (Fraction(1),))
            vanishing = add(primitive,
                            {(0,): -evaluate(primitive, (Fraction(1),)) + constant,
                             (1,): -constant})
            require(evaluate(vanishing, (Fraction(1),)) == 0
                    and evaluate(vanishing, (Fraction(2),)) == 0,
                    "a relative primitive has a nonzero endpoint")
            reduced = add({(0,): constant}, {(-1,): log_coefficient}, derivative(vanishing, 0))
            require(reduced == integrand, "complete rational one-variable reduction failed")
            minimum = min((exponent for (exponent,) in vanishing), default=0)
            clearing_power = max(0, -minimum)
            polynomial = shift(vanishing, (clearing_power,))
            quotient, remainder = divide_univariate_polynomial(polynomial, endpoint_factor)
            require(not remainder, "the primitive lacks the complete endpoint factor")
            rational_quotient = shift(quotient, (-clearing_power,))
            require(multiply(endpoint_factor, rational_quotient) == vanishing,
                    "endpoint-factor division was not exact")
            row = {"numerator_x_power": power, "denominator_1_plus_x_power": pole,
                   "rational_value": str(constant), "log_2_coefficient": str(log_coefficient),
                   "primitive_in_u_equals_1_plus_x": encode_poly(vanishing),
                   "primitive_divided_by_x_times_1_minus_x": encode_poly(rational_quotient)}
            trace.append(row)
            if (power, pole) in {(0, 0), (0, 1), (0, 2), (1, 1), (7, 5), (12, 12)}:
                selected.append(row)
    return {"numerator_powers": [0, 12], "pole_orders": [0, 12],
            "complete_identity_count": len(trace), "all_identities_sha256": digest(trace),
            "selected_complete_identities": selected,
            "proved_relative_basis": ["1", "1/(1+x)"],
            "basis_periods": ["1", "log(2)"],
            "basis_scope": "Only the one-variable relative example; not an original-core master basis."}


def polynomial_cube_controls():
    trace = []
    for variables in (2, 3):
        for powers in product(range(5), repeat=variables):
            integral = prod(Fraction(1, power + 1) for power in powers)
            divergence = {}
            primitives = []
            for coordinate in range(variables):
                coefficient = Fraction(1, powers[coordinate] + 1) * prod(
                    Fraction(1, powers[index] + 1) for index in range(coordinate))
                first = tuple(0 if index < coordinate else
                              powers[index] + 1 if index == coordinate else powers[index]
                              for index in range(variables))
                second = tuple(0 if index < coordinate else
                               1 if index == coordinate else powers[index]
                               for index in range(variables))
                primitive = add({first: coefficient}, {second: -coefficient})
                for endpoint in (Fraction(0), Fraction(1)):
                    restricted = {}
                    for exponents, value in primitive.items():
                        remaining = exponents[:coordinate] + exponents[coordinate + 1:]
                        restricted[remaining] = restricted.get(remaining, Fraction()) + value * endpoint ** exponents[coordinate]
                    require(not clean(restricted), "polynomial cube primitive has a boundary value")
                divergence = add(divergence, derivative(primitive, coordinate))
                primitives.append(encode_poly(primitive))
            expected = add({powers: Fraction(1)}, {(0,) * variables: -integral})
            require(divergence == expected, "polynomial cube reduction is not exact")
            trace.append({"powers": powers, "integral": str(integral), "primitives": primitives})
    return {"dimensions": [2, 3], "each_coordinate_power": [0, 4],
            "complete_identity_count": len(trace), "all_identities_sha256": digest(trace),
            "last_complete_identity": trace[-1]}


def build_certificate():
    graphs = [graph_control(*graph) for graph in GRAPHS]
    certificate = {
        "schema": "fixed-core-period-stabilization-v1",
        "claim": "Finite exact controls of the written all-order theorem; no general core master basis computed.",
        "arithmetic": "Python integers and fractions.Fraction; no floating-point evaluation.",
        "production_module_imports": [],
        "jacobian_controls": jacobian_controls(),
        "graph_controls": graphs,
        "signed_numerator_controls": signed_numerator_controls(),
        "cycle_controls": cycle_controls(),
        "one_variable_relative_controls": univariate_controls(),
        "polynomial_cube_controls": polynomial_cube_controls(),
        "summary": {
            "graphs": len(graphs),
            "complete_hepp_sectors": sum(row["sector_count"] for row in graphs),
            "cotree_factorization_checks": sum(row["cotree_factorization_checks"] for row in graphs),
            "tested_scalar_weights": sum(row["tested_weight_count"] for row in graphs),
            "rejected_scalar_weight_witnesses": sum(row["rejected_weight_witness_count"] for row in graphs),
            "convergent_sector_coordinate_checks": sum(row["convergent_sector_coordinate_checks"] for row in graphs),
            "convergence_preserving_cotree_shift_terms": sum(row["convergence_preserving_shift_terms"] for row in graphs),
        },
        "source_sha256": {name: sha256((HERE / name).read_bytes()).hexdigest()
                          for name in ("README.md", "DERIVATION.md", "verify.py", "run_checks.py")},
    }
    # Normalize tuples before equality with the persisted JSON representation.
    return json.loads(canonical(certificate))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate the exact certificate")
    arguments = parser.parse_args()
    calculated = build_certificate()
    if arguments.write:
        temporary = CERTIFICATE.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(calculated, indent=2, sort_keys=True) + "\n")
        temporary.replace(CERTIFICATE)
        print("Wrote exact stabilization controls.")
    else:
        require(CERTIFICATE.exists(), "missing checks.json; run verify.py --write first")
        recorded = json.loads(CERTIFICATE.read_text())
        require(recorded == calculated, "certificate changed; inspect the difference before regenerating")
        print("Replayed every exact stabilization control and source digest.")
    for name, count in calculated["summary"].items():
        print(f"{name}: {count}")
    print("signed_numerator_sectors:", calculated["signed_numerator_controls"]["checked_sector_count"])
    print("independently_integrated_cycle_weights:", calculated["cycle_controls"]["integrated_weight_count"])
    print("one_variable_complete_relative_identities:", calculated["one_variable_relative_controls"]["complete_identity_count"])
    print("polynomial_cube_complete_relative_identities:", calculated["polynomial_cube_controls"]["complete_identity_count"])
    print("canonical_certificate_sha256:", digest(calculated))
    print("No general core master basis, numerical independence, or stabilization threshold was computed.")


if __name__ == "__main__":
    main()
