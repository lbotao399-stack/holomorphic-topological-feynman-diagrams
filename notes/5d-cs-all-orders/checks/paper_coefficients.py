#!/usr/bin/env python3
"""Independent exact arithmetic assembly for the five-dimensional paper.

Reads the four core-class period rows and marked-theta rows from the
research notes, evaluates the half-integer Gamma factors independently,
and verifies the census and every automorphism order from the chain data.
It does not import the original assembly or graph-enumeration programs,
and does not claim to reintegrate the four-loop periods.
"""

import ast
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import factorial, prod
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
HIGHER = ROOT / 'deepen/o3/higher_orders'


class Polynomial:
    """One-variable polynomial with rational coefficients, in ascending order."""

    def __init__(self, values):
        values = [Fraction(x) for x in values]
        while len(values) > 1 and values[-1] == 0:
            values.pop()
        self.values = tuple(values or [Fraction(0)])

    @staticmethod
    def number(value):
        return Polynomial([value])

    def __add__(self, other):
        return Polynomial([self.coefficient(i) + other.coefficient(i)
                           for i in range(max(len(self.values), len(other.values)))])

    def __neg__(self):
        return Polynomial([-x for x in self.values])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        out = [Fraction(0)] * (len(self.values) + len(other.values) - 1)
        for i, a in enumerate(self.values):
            for j, b in enumerate(other.values):
                out[i+j] += a*b
        return Polynomial(out)

    def scale(self, value):
        return Polynomial([x * value for x in self.values])

    def coefficient(self, degree):
        return self.values[degree] if degree < len(self.values) else Fraction(0)

    def evaluate(self, value):
        result = Fraction(0)
        for coefficient in reversed(self.values):
            result = result * value + coefficient
        return result

    def __eq__(self, other):
        return self.values == other.values


ZERO = Polynomial.number(0)
ONE = Polynomial.number(1)
VARIABLE = Polynomial([0, 1])


def parse_cell(cell, variable):
    """A restricted arithmetic parser; no eval, functions, or arbitrary names."""
    text = cell.strip().strip('$').replace(variable, '*x')
    if text.startswith('*'):
        text = text[1:]
    text = text.replace('(*', '(1*').replace('+*', '+1*').replace('-*', '-1*')
    expression = ast.parse(text, mode='eval')

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return Polynomial.number(node.value)
        if isinstance(node, ast.Name) and node.id == 'x':
            return VARIABLE
        if isinstance(node, ast.UnaryOp):
            result = visit(node.operand)
            if isinstance(node.op, ast.USub):
                return -result
            if isinstance(node.op, ast.UAdd):
                return result
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                assert len(right.values) == 1 and right.values[0]
                return left.scale(1 / right.values[0])
        raise ValueError(f'Unrecognized arithmetic in table cell: {cell}')

    return visit(expression)


def gamma_half(twice_argument):
    """Coefficient of sqrt(pi) in Gamma(twice_argument/2), exactly."""
    assert twice_argument > 0 and twice_argument % 2 == 1
    integer = (twice_argument - 1) // 2
    return Fraction(factorial(2*integer), 4**integer * factorial(integer))


def theta_coefficients():
    # H=2, D=3. The Dirichlet/Gamma formula leaves one factor pi.
    rows = {}
    for shape in ((1, 3, 3), (2, 2, 3)):
        coefficient = (Fraction(1, 2)
                       * prod(gamma_half(7-2*n) for n in shape)
                       / gamma_half(7)
                       / prod(factorial(n-1) for n in shape))
        positive_cycle_sign = (-1)**(3*(shape[1]+1))
        signed_coefficient = -positive_cycle_sign * coefficient
        automorphisms = 2 * prod(factorial(n) for n in Counter(shape).values())
        rows[''.join(map(str, shape))] = (signed_coefficient, automorphisms)
    assert rows == {'133': (Fraction(-1, 20), 4), '223': (Fraction(1, 30), 4)}
    assembled = {key: value / automorphisms for key, (value, automorphisms) in rows.items()}
    # Substituting F133=beta/2-F223 in Killing normalization.
    f223 = assembled['223'] - assembled['133']
    beta = assembled['133'] / 2
    assert f223 == Fraction(1, 48) and beta == Fraction(-1, 160)
    assert -f223 / 64 == Fraction(-1, 3072)
    return {'signed_kernel_over_pi': {k: str(v[0]) for k, v in rows.items()},
            'automorphism_orders': {k: v[1] for k, v in rows.items()},
            'assembly_in_units_of_pi_nu': {k: str(v) for k, v in assembled.items()},
            'killing_reduction': {'F223': str(f223), 'beta': str(beta)},
            'anomaly_class_in_units_of_U_over_pi4': '-1/3072'}


def period_rows(text):
    rows = {}
    for line in text.splitlines():
        match = re.match(r'^\|\s*(\d+)[\u2013-](\d+)\s*\|', line)
        if not match:
            continue
        start, finish = map(int, match.groups())
        cells = [x.strip() for x in line.strip().strip('|').split('|')]
        assert len(cells) == 4
        rows[start, finish] = tuple(parse_cell(x, r'\pi^2') for x in cells[1:])
    assert set(rows) == {(0, 2), (3, 10), (11, 38), (39, 49)}
    assert all(len(value.values) <= 2 for row in rows.values() for value in row)
    return rows


def marked_rows(text):
    rows = {}
    for line in text.splitlines():
        if not (line.startswith('| $(1,3,3)$ |') or line.startswith('| $(2,2,3)$ |')):
            continue
        cells = [x.strip() for x in line.strip().strip('|').split('|')]
        assert len(cells) == 3
        rows[cells[0].strip('$')] = tuple(parse_cell(x, 'N^2') for x in cells[1:])
    assert set(rows) == {'(1,3,3)', '(2,2,3)'}
    return rows


def fourth_order(periods, marked):
    bare = tuple(sum((row[rank] for row in periods.values()), ZERO) for rank in range(3))
    expected_bare = (Polynomial([Fraction(-172, 945)]),
                     Polynomial([Fraction(139229, 48600), Fraction(-35701, 120960)]),
                     Polynomial([Fraction(760981, 3402000), Fraction(-130727, 5806080)]))
    assert bare == expected_bare
    # Couplings from the second- and third-order equations, divided by
    # epsilon*pi*nu. Here VARIABLE=N^2; no original assembly code is used.
    a_squared = VARIABLE.scale(Fraction(-1, 4))
    second = Polynomial([10, 3]).scale(Fraction(-1, 20))
    third = Polynomial([5, 2]).scale(Fraction(-1, 30))
    tree = (a_squared * third).scale(10) + (second * second).scale(-4)
    one_mark = sum((row[0] for row in marked.values()), ZERO)
    two_marks = sum((row[1] for row in marked.values()), ZERO)
    theta = second * one_mark + a_squared * two_marks
    tree_and_marked_primitive = (tree + theta).scale(Fraction(1, 16))
    assert tree_and_marked_primitive == Polynomial([192480, 85304, 21221]).scale(Fraction(1, 3628800))
    coefficients = tuple(Polynomial.number(tree_and_marked_primitive.coefficient(rank))
                         + bare[rank].scale(Fraction(1, 16)) for rank in range(3))
    expected = (Polynomial([Fraction(1, 24)]),
                Polynomial([Fraction(1102559, 5443200), Fraction(-35701, 1935360)]),
                Polynomial([Fraction(4216, 212625), Fraction(-130727, 92897280)]))
    assert coefficients == expected
    # Each nonconstant coefficient decreases with pi²; pi²<10 gives
    # strictly positive, rational lower bounds. No decimal evaluation.
    lower = [value.evaluate(Fraction(10)) for value in coefficients]
    assert all(value > 0 for value in lower)
    assert coefficients[0].coefficient(1) == 0
    assert coefficients[1].coefficient(1) < 0 and coefficients[2].coefficient(1) < 0

    def show(rows):
        return {str(2*degree): {'rational': str(value.coefficient(0)),
                               'pi_squared': str(value.coefficient(1))}
                for degree, value in enumerate(rows)}

    return {'bare_core_sum_over_pi_squared': show(bare),
            'trees_by_rank_power': {str(2*i): str(value) for i, value in enumerate(tree.values)},
            'marked_thetas_by_rank_power': {str(2*i): str(value) for i, value in enumerate(theta.values)},
            'trees_and_marked_after_division_by_16':
                {str(2*i): str(value) for i, value in enumerate(tree_and_marked_primitive.values)},
            'primitive_coefficients_in_units_of_epsilon_pi_nu_squared': show(coefficients),
            'strict_rational_lower_bounds_from_pi_squared_less_than_10':
                {str(2*i): str(value) for i, value in enumerate(lower)},
            'pi_squared_column_nonzero_for_every_positive_rank': True}


def low_order_coefficients():
    # VARIABLE=N^2. Each term comes separately from the bracket and
    # anomaly boundary (order two), then the tree and marked graph (three).
    second_tree = (VARIABLE + Polynomial.number(4)).scale(Fraction(1, 512))
    second_graph = VARIABLE.scale(Fraction(1, 2560))
    second = second_tree + second_graph
    assert second == Polynomial([10, 3]).scale(Fraction(1, 1280))
    third_tree = second.scale(Fraction(1, 3))
    third_marked = VARIABLE.scale(Fraction(1, 3840))
    third = third_tree + third_marked
    assert third == Polynomial([5, 2]).scale(Fraction(1, 1920))
    assert second.evaluate(4) == Fraction(11, 640)
    return {'i2_over_epsilon_U_pi_minus4': {'constant': str(second.coefficient(0)),
                                           'N_squared': str(second.coefficient(1))},
            'i3_over_a_epsilon_U_pi_minus4': {'constant': str(third.coefficient(0)),
                                             'N_squared': str(third.coefficient(1))},
            'rank_two_i2_coefficient': '11/640',
            'each_tree_and_graph_term_assembled_separately': True}


def vogel_checks():
    n = VARIABLE
    def number(value):
        return Polynomial.number(value)
    def product(rows):
        out = ONE
        for row in rows:
            out = out * row
        return out
    families = {
        'sl_N': ((number(-2), number(2), n), n*n-ONE,
                 ((n*n-ONE)*(n*n-number(4))).scale(-3), n*n),
        'so_N': ((number(-2), number(4), n-number(4)), (n*(n-ONE)).scale(Fraction(1, 2)),
                 product((n,n-ONE,n-number(3),n-number(4),n+number(2))).scale(-3),
                 product((n-number(2),)*3).scale(2)),
        'sp_N': ((number(-2), ONE, (n+number(4)).scale(Fraction(1, 2))),
                 (n*(n+ONE)).scale(Fraction(1, 2)),
                 product((n,n-number(2),n+ONE,n+number(3),n+number(4))).scale(-3),
                 product((n+number(2),)*3).scale(2)),
    }
    verified = {}
    for name, (parameters, dimension, expected_numerator, expected_denominator) in families.items():
        total = sum(parameters, ZERO)
        denominator = product((total, total, total)).scale(2)
        numerator = (dimension*product(tuple(total+v for v in parameters))).scale(-3)
        assert numerator*expected_denominator == expected_numerator*denominator
        # Check the universal dimension and the unfactored Casimir formula.
        assert dimension*product(parameters) == product(tuple(v-total.scale(2) for v in parameters))
        casimir_numerator = dimension*(
            (total*sum((v*v for v in parameters), ZERO)).scale(3)
            -sum((v*v*v for v in parameters), ZERO)-product((total, total, total)).scale(8))
        assert casimir_numerator == numerator
        verified[name] = 'polynomial identities after denominator clearing'
    exceptional = {
        'G2': ((Fraction(-2), Fraction(10, 3), Fraction(8, 3)), 14, Fraction(-385, 12)),
        'F4': ((Fraction(-2), Fraction(5), Fraction(6)), 52, Fraction(-12740, 81)),
        'E6': ((Fraction(-2), Fraction(6), Fraction(8)), 78, Fraction(-975, 4)),
        'E7': ((Fraction(-2), Fraction(8), Fraction(12)), 133, Fraction(-34580, 81)),
        'E8': ((Fraction(-2), Fraction(12), Fraction(20)), 248, Fraction(-12152, 15)),
    }
    for name, (parameters, dimension, expected) in exceptional.items():
        total = sum(parameters)
        assert prod(v-2*total for v in parameters)/prod(parameters) == dimension
        detector = -3*dimension*prod(total+v for v in parameters)/(2*total**3)
        assert detector == expected
        verified[name] = str(detector)
    return {'detector_8n4_minus_3dim': verified,
            'sp_N_convention': 'N is the even defining dimension; parameters (-2,1,(N+4)/2)'}


def graph_census(records):
    assert len(records) == 50
    all_core_permutations = tuple(permutations(range(6)))
    core_groups = defaultdict(list)
    weighted_keys = set()
    checked_subsets = 0
    for index, record in enumerate(records):
        chains = record['chains']
        assert len(chains) == 9
        assert sorted(v for chain in chains for v in chain[1:-1]) == [6, 7, 8, 9]
        core_edges = [(min(chain[0], chain[-1]), max(chain[0], chain[-1])) for chain in chains]
        degrees = Counter(v for edge in core_edges for v in edge)
        assert degrees == Counter({v: 3 for v in range(6)})
        full_edges = [tuple(sorted((a, b))) for chain in chains for a, b in zip(chain, chain[1:])]
        assert len(full_edges) == 13 and len(set(full_edges)) == 13

        def signature(order, include_lengths):
            entries = []
            for chain in chains:
                ends = sorted((order[chain[0]], order[chain[-1]]))
                entries.append(tuple(ends) + ((len(chain)-1,) if include_lengths else ()))
            return tuple(sorted(entries))

        core_key = min(signature(order, False) for order in all_core_permutations)
        weighted_key = min(signature(order, True) for order in all_core_permutations)
        assert weighted_key not in weighted_keys
        weighted_keys.add(weighted_key)
        core_groups[core_key].append(index)
        initial = signature(tuple(range(6)), True)
        core_automorphisms = sum(signature(order, True) == initial for order in all_core_permutations)
        identical_path_maps = prod(factorial(multiplicity) for multiplicity in Counter(initial).values())
        assert core_automorphisms * identical_path_maps == record['aut']
        for size in range(2, 11):
            for subset in combinations(range(10), size):
                vertices = set(subset)
                internal_edges = sum(a in vertices and b in vertices for a, b in full_edges)
                assert 2*internal_edges <= 3*size-4
                checked_subsets += 1
    groups = sorted(core_groups.values(), key=lambda values: min(values))
    assert [len(group) for group in groups] == [3, 8, 28, 11]
    assert groups == [list(range(0, 3)), list(range(3, 11)), list(range(11, 39)), list(range(39, 50))]
    return {'graphs': len(records), 'core_classes': len(groups),
            'class_counts': [len(group) for group in groups],
            'all_graphs_distinct': True, 'all_automorphism_orders_verified': True,
            'induced_subset_inequalities_checked': checked_subsets,
            'automorphism_order_histogram': dict(sorted(Counter(row['aut'] for row in records).items()))}


def main():
    inputs = {'higher_order_notes': HIGHER / 'README.md',
              'period_rows': HIGHER / 'four_loop_reduction/README.md',
              'graph_census': HIGHER / 'four-loop-graphs.json'}
    result = {'arithmetic': 'exact rational arithmetic; no floating point',
              'scope': 'Independent assembly of supplied period rows, not reintegration of four-loop periods.',
              'theta': theta_coefficients(),
              'low_order': low_order_coefficients(),
              'vogel': vogel_checks(),
              'fourth_order': fourth_order(period_rows(inputs['period_rows'].read_text()),
                                           marked_rows(inputs['higher_order_notes'].read_text())),
              'census': graph_census(json.loads(inputs['graph_census'].read_text())),
              'source_sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                                for path in inputs.values()}}
    destination = HERE / 'verification/paper_coefficients.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
