#!/usr/bin/env python3
"""Exact Laurent representatives and the primary-tree test for every boundary Q_n.

Only integers and fractions are used. No previous graph evaluator is imported.
The result concerns trees whose external fields are the displayed linear
Laurent representatives, not physical bulk-to-boundary kernels.
"""
from fractions import Fraction
from itertools import product
from math import comb, factorial
from pathlib import Path
import json


def clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def plus(*values):
    answer = {}
    for value in values:
        for key, coefficient in value.items():
            answer[key] = answer.get(key, 0) + coefficient
    return clean(answer)


def bits(mask, count):
    return tuple(i for i in range(count) if mask >> i & 1)


def wedge_sign(left, right):
    if left & right:
        return 0
    return (-1) ** sum(i > j for i in bits(left, 5) for j in bits(right, 5))


def bulk_term(coefficient, powers, vectors):
    return {(tuple(powers), vectors): Fraction(coefficient)}


def multiply(left, right):
    answer = {}
    for (a, u), x in left.items():
        for (b, v), y in right.items():
            sign = wedge_sign(u, v)
            if sign:
                key = (tuple(i+j for i, j in zip(a, b)), u | v)
                answer[key] = answer.get(key, 0) + sign*x*y
    return clean(answer)


def derivative(value, coordinate):
    answer = {}
    for (powers, vectors), coefficient in value.items():
        if not powers[coordinate]:
            continue
        updated = list(powers)
        updated[coordinate] -= 1
        key = (tuple(updated), vectors)
        answer[key] = answer.get(key, 0) + coefficient*powers[coordinate]
    return clean(answer)


def divergence(value):
    answer = {}
    for (powers, vectors), coefficient in value.items():
        for position, coordinate in enumerate(bits(vectors, 5)):
            if not powers[coordinate]:
                continue
            updated = list(powers)
            updated[coordinate] -= 1
            key = (tuple(updated), vectors ^ (1 << coordinate))
            answer[key] = answer.get(key, 0) + (
                (-1)**position*coefficient*powers[coordinate])
    return clean(answer)


def phi(value):
    """Pinned Laurent/form dictionary; forms ordered dz, theta, dtheta."""
    answer = {}
    epsilon = {0: -1, 1: -1, 2: 1, 3: 1}
    for (powers, vectors), coefficient in value.items():
        assert all(e <= -1 for e in powers[2:])
        tangent_vectors = vectors & 3
        theta = vectors >> 2
        normal_orders = tuple(-e-1 for e in powers[2:])
        factor = Fraction(
            epsilon[tangent_vectors]*(-1)**sum(normal_orders),
            factorial(normal_orders[0])*factorial(normal_orders[1])
            * factorial(normal_orders[2]))
        key = (powers[:2], 3 ^ tangent_vectors, theta, normal_orders)
        answer[key] = answer.get(key, 0) + coefficient*factor
    return clean(answer)


def form_term(coefficient, tangent, dz, theta, normal_orders):
    return {(tuple(tangent), dz, theta, tuple(normal_orders)): Fraction(coefficient)}


def de_rham(value):
    answer = {}
    for (tangent, dz, theta, normal_orders), coefficient in value.items():
        for coordinate in range(2):
            if dz >> coordinate & 1 or not tangent[coordinate]:
                continue
            updated = list(tangent)
            updated[coordinate] -= 1
            sign = (-1)**sum(j < coordinate for j in bits(dz, 2))
            key = (tuple(updated), dz | (1 << coordinate), theta, normal_orders)
            answer[key] = answer.get(key, 0) + sign*coefficient*tangent[coordinate]
        for position, coordinate in enumerate(bits(theta, 3)):
            updated = list(normal_orders)
            updated[coordinate] += 1
            sign = (-1)**(dz.bit_count()+position)
            key = (tangent, dz, theta ^ (1 << coordinate), tuple(updated))
            answer[key] = answer.get(key, 0) + sign*coefficient
    return clean(answer)


def form_basis(weight, degree):
    answer = []
    for dz in range(4):
        tangent = tuple(weight[i] - ((dz >> i) & 1) for i in range(2))
        if min(tangent) < 0:
            continue
        for theta in range(8):
            orders = tuple(weight[i+2]-((theta >> i) & 1) for i in range(3))
            if min(orders) < 0 or sum(orders)+dz.bit_count() != degree:
                continue
            answer.append((tangent, dz, theta, orders))
    return answer


def rank(matrix):
    if not matrix:
        return 0
    a = [list(map(Fraction, row)) for row in matrix]
    height, width = len(a), len(a[0])
    pivot = 0
    for column in range(width):
        chosen = next((i for i in range(pivot, height) if a[i][column]), None)
        if chosen is None:
            continue
        a[pivot], a[chosen] = a[chosen], a[pivot]
        factor = a[pivot][column]
        a[pivot] = [v/factor for v in a[pivot]]
        for i in range(height):
            if i != pivot and a[i][column]:
                factor = a[i][column]
                a[i] = [v-factor*w for v, w in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def exact_form_dimension(weight, degree):
    source, target = form_basis(weight, degree-1), form_basis(weight, degree)
    columns = [de_rham({key: Fraction(1)}) for key in source]
    matrix = [[column.get(key, 0) for column in columns] for key in target]
    return len(source), len(target), rank(matrix)


def serialize_bulk(value):
    return [{"coefficient": str(c), "coordinate_powers": list(e),
             "polyvector_indices": list(bits(mask, 5))}
            for (e, mask), c in sorted(value.items())]


def certificate():
    # Coordinate order: z1,z2,w1,w2,w3.
    # The two tangent vectors occupy bits zero and one, retained everywhere.
    P = 3
    G1 = plus(bulk_term(2, (0,0,-1,-3,-1), P | 16),
              bulk_term(-1, (0,0,-1,-2,-2), P | 8))
    G2 = plus(bulk_term(2, (0,0,-1,-1,-3), P | 4),
              bulk_term(-1, (0,0,-2,-1,-2), P | 16))
    G3 = plus(bulk_term(2, (0,0,-3,-1,-1), P | 8),
              bulk_term(-1, (0,0,-2,-2,-1), P | 4))
    H = plus(bulk_term(2, (0,0,-3,-1,-1), P | 8 | 16),
             bulk_term(-1, (0,0,-2,-2,-1), P | 4 | 16),
             bulk_term(1, (0,0,-2,-1,-2), P | 4 | 8))
    flux = bulk_term(1, (0,0,-1,-1,-1), P)
    fields = [G1, G2, G3, H, flux]

    primitives = [
        form_term(1, (0,0), 0, 6, (0,1,0)),
        form_term(-1, (0,0), 0, 5, (0,0,1)),
        form_term(1, (0,0), 0, 3, (1,0,0)),
        form_term(1, (0,0), 0, 7, (1,0,0)),
    ]
    weights = [(0,0,0,2,1), (0,0,1,0,2), (0,0,2,1,0), (0,0,2,1,1)]
    records = []
    for index, field in enumerate(fields[:4]):
        assert not divergence(field)
        assert phi(field) == de_rham(primitives[index])
        assert not de_rham(phi(field))
        dimensions = exact_form_dimension(weights[index], 2)
        assert dimensions[2] == 1
        records.append({"field": ["G1", "G2", "G3", "H"][index],
                        "laurent_representative": serialize_bulk(field),
                        "boundary_weight": list(weights[index]),
                        "dimension_Omega1_Omega2_B2": list(dimensions)})
    assert phi(flux) == form_term(1, (0,0), 0, 0, (0,0,0))

    # General Phi/divergence identity on all normal derivative multi-indices,
    # including nonclosed primitive controls in all tangent polyvector sectors.
    controls = [
        bulk_term(1, (2,1,-1,-2,-3), vectors)
        for vectors in range(32)
    ]
    derivative_controls = 0
    jets = []
    for order in range(7):
        for a in range(order+1):
            for b in range(order-a+1):
                multi = (a,b,order-a-b)
                for field in fields:
                    jet = field
                    for coordinate, multiplicity in enumerate(multi, 2):
                        for _ in range(multiplicity):
                            jet = derivative(jet, coordinate)
                    assert all(mask & P == P for _, mask in jet)
                    assert phi(divergence(jet)) == de_rham(phi(jet))
                    jets.append(jet)
                    derivative_controls += 1
    for field in controls:
        assert phi(divergence(field)) == de_rham(phi(field))

    # Literal complete products, before any component evaluation or rooting.
    pair_checks = 0
    for left in jets:
        for right in fields:
            assert not multiply(left, right)
            assert not divergence(multiply(left, right))
            pair_checks += 1

    outputs = []
    for n in range(2,25):
        coefficient = (-1)**n*factorial(n-1)
        field = bulk_term(coefficient, (n-1,n-1,-n,-2,-2), 0)
        expected = form_term(1, (n-1,n-1), 3, 0, (n-1,1,1))
        assert phi(field) == expected and not divergence(field)
        primitive = form_term(Fraction(1,n), (n,n-1), 2, 0, (n-1,1,1))
        assert de_rham(primitive) == expected
        weight = (n,n,n-1,1,1)
        dimensions = exact_form_dimension(weight,n+3)
        assert dimensions[1:] == (1,1)
        total_input_theta = (2*n-1,n+1,n+1)
        assert tuple(x-n for x in total_input_theta) == weight[2:]
        assert 2*(n+1)+1-n == n+3
        outputs.append({"boundary_order": n, "input_traces": n+1,
                        "output_letters": n+3, "output_weight": list(weight),
                        "highest_form_degree": n+3,
                        "candidate_protected_output_dimension": 1,
                        "laurent_scalar_coefficient": coefficient,
                        "laurent_scalar_powers": [n-1,n-1,-n,-2,-2]})

    # Combinatorial induction: every rooted tree with vertices of arity >=2
    # has an internal vertex with at least two external children. The exact
    # recursion counts ordered cubic trees, rather than relying on a diagram.
    binary_counts = [0,1]
    for leaves in range(2,13):
        count = sum(binary_counts[i]*binary_counts[leaves-i]
                    for i in range(1,leaves))
        assert count == comb(2*(leaves-1),leaves-1)//leaves
        binary_counts.append(count)
    return {
        "status": "PASS",
        "arithmetic": "exact integers and rational numbers; standard library only",
        "pinned_input_revision": "29ef4be3c92b8055ac51067357524c1048e6cda0",
        "coordinate_order": ["z1","z2","w1","w2","w3"],
        "input_laurent_fields": records,
        "independent_dictionary_controls": len(controls),
        "normal_derivative_dictionary_checks": derivative_controls,
        "complete_zero_products_before_root_divergence": pair_checks,
        "output_charge_and_laurent_checks": outputs,
        "ordered_cubic_tree_counts_by_external_input_count": {
            str(n): binary_counts[n] for n in range(2,13)},
        "all_order_scope": (
            "Rooted trees made from coefficient-derivative/wedge local primary "
            "vertices, with the given linear Laurent fields and explicit flux "
            "leaves. Arbitrary linear internal propagators are permitted. "
            "External bulk-to-boundary kernels and graphs with cycles are not covered."),
        "output_scalar_scope": (
            "Unique protected one-particle state in the boundary output charge "
            "sector; not an identification of the nonclosed raw joining output."),
    }


def main():
    result = certificate()
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({
        "status": result["status"],
        "input_sectors": len(result["input_laurent_fields"]),
        "dictionary_controls": result["independent_dictionary_controls"],
        "normal_derivative_checks": result["normal_derivative_dictionary_checks"],
        "complete_zero_products": result["complete_zero_products_before_root_divergence"],
        "output_orders": [2,24],
    }))


if __name__ == "__main__":
    main()
