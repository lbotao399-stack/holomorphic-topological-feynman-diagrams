"""Complete reductions, checked by independent symbolic differentiation."""
from fractions import Fraction
import hashlib
from itertools import permutations
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
from htk.graph import Graph
from htk.core_sector_periods import CoreSector, CubePeriodRelations


def build():
    cases = []
    cube = CubePeriodRelations({(0,): 1, (1,): 1})
    result = cube.search(cube.monomial((7,), 5),
                        [cube.monomial((0,)), cube.monomial((0,), 1)],
                        max_primitive_power=4, max_primitive_degree=8)
    cases.append({'name': 'one_variable_nonrational', 'witness': result.record()})
    cube = CubePeriodRelations({(0, 0): 1, (1, 1): 1})
    result = cube.search(cube.monomial((2, 2), 3),
                        [cube.monomial((0, 0)), cube.monomial((0, 0), 1),
                         cube.monomial((1, 0), 1)],
                        max_primitive_power=2, max_primitive_degree=4)
    cases.append({'name': 'two_variable_nonrational', 'witness': result.record()})
    graph = Graph(2, [(0, 1)] * 3)
    for order in permutations(range(3)):
        sector = CoreSector(graph, order)
        result = sector.relations.search(sector.pullback((2, 2, 2), 3),
                        [sector.relations.monomial((0, 0))],
                        max_primitive_power=2, max_primitive_degree=4)
        cases.append({'name': 'three_edge_banana', 'edge_order': list(order),
                      'witness': result.record()})
    return {'schema': 'htk.core-sector-controls.v1', 'cases': cases,
            'scope': 'Complete specified target identities, not a universal numerical basis.'}


def _expression(rows, variables):
    # Construct exact expressions from coefficient/exponent data; do not execute
    # certificate text or ask SymPy to integrate the production generator.
    return sum((sp.Rational(Fraction(c).numerator, Fraction(c).denominator)
                * sp.prod(x ** a for x, a in zip(variables, powers))
                for powers, c in rows), sp.Integer(0))


def verify(data):
    if data['schema'] != 'htk.core-sector-controls.v1':
        raise ArithmeticError('unexpected controls schema')
    counts = {'complete_certificates': 0, 'independent_rational_identities': 0,
              'independent_primitive_faces': 0, 'full_graph_sectors': 0,
              'intended_target_bindings': 0}
    banana = Fraction()
    orders = set()
    names = []
    for case in data['cases']:
        record = case['witness']
        denominator = {tuple(a): Fraction(c) for a, c in record['denominator']}
        cube = CubePeriodRelations(denominator)
        cube.read_certificate(record)
        variables = sp.symbols('x:' + str(cube.variables))
        U = _expression(record['denominator'], variables)
        def form(value):
            return _expression(value['numerator'], variables) / U ** value['power']
        # Bind the algebraic witness to the specified integral, not merely to
        # whichever target and master forms a certificate happens to contain.
        # These densities are written independently of the production sector
        # pullback and boundary-derivative generators.
        if case['name'] == 'one_variable_nonrational':
            if len(variables) != 1:
                raise ArithmeticError('wrong dimension for the one-variable control')
            x, = variables
            intended_U = 1 + x
            intended_target = x ** 7 / intended_U ** 5
            intended_masters = [sp.Integer(1), 1 / intended_U]
        elif case['name'] == 'two_variable_nonrational':
            if len(variables) != 2:
                raise ArithmeticError('wrong dimension for the two-variable control')
            x, y = variables
            intended_U = 1 + x * y
            intended_target = x ** 2 * y ** 2 / intended_U ** 3
            intended_masters = [sp.Integer(1), 1 / intended_U, x / intended_U]
        elif case['name'] == 'three_edge_banana':
            if len(variables) != 2:
                raise ArithmeticError('wrong dimension for a banana sector')
            x, y = variables
            intended_U = 1 + x + x * y
            intended_target = x / intended_U ** 3
            intended_masters = [sp.Integer(1)]
        else:
            raise ArithmeticError('unrecognized control')
        if (sp.cancel(U - intended_U) != 0
                or sp.cancel(form(record['target']) - intended_target) != 0
                or len(record['masters']) != len(intended_masters)
                or any(sp.cancel(form(actual) - expected) != 0
                       for actual, expected in zip(record['masters'], intended_masters))):
            raise ArithmeticError('certificate does not bind the intended target and masters')
        counts['intended_target_bindings'] += 1
        remainder = form(record['target'])
        for coefficient, master in zip(record['coefficients'], record['masters']):
            remainder -= sp.Rational(coefficient) * form(master)
        for coefficient, coordinate, powers, q in record['primitives']:
            x = variables[coordinate]
            primitive = sp.Rational(coefficient) * x * (1 - x)
            primitive *= sp.prod(z ** a for z, a in zip(variables, powers)) / U ** q
            remainder -= sp.diff(primitive, x)
            for endpoint in (0, 1):
                if sp.cancel(primitive.subs(x, endpoint)) != 0:
                    raise ArithmeticError('primitive leaves a nonzero cube-face contribution')
                counts['independent_primitive_faces'] += 1
        if sp.cancel(remainder) != 0:
            raise ArithmeticError('independent full rational identity failed')
        counts['complete_certificates'] += 1
        counts['independent_rational_identities'] += 1
        coefficients = tuple(Fraction(c) for c in record['coefficients'])
        names.append(case['name'])
        if case['name'] == 'one_variable_nonrational':
            if coefficients != (Fraction(1553, 64), Fraction(-35)):
                raise ArithmeticError('one-variable exact value changed')
        elif case['name'] == 'two_variable_nonrational':
            if coefficients != (Fraction(-1, 2), Fraction(1), Fraction(-3, 4)):
                raise ArithmeticError('two-variable exact value changed')
        elif case['name'] == 'three_edge_banana':
            if coefficients != (Fraction(1, 12),):
                raise ArithmeticError('banana sector is not complete')
            order = tuple(case['edge_order'])
            if order in orders:
                raise ArithmeticError('duplicate graph sector')
            orders.add(order)
            banana += coefficients[0]
            counts['full_graph_sectors'] += 1
        else:
            raise ArithmeticError('unrecognized control')
    if names.count('one_variable_nonrational') != 1 or names.count('two_variable_nonrational') != 1:
        raise ArithmeticError('a full nonrational control is missing')
    if orders != set(permutations(range(3))) or banana != Fraction(1, 2):
        raise ArithmeticError('complete original graph failed sector reconstruction')
    counts['three_edge_banana_value'] = str(banana)
    return counts


if __name__ == '__main__':
    if '--write' in sys.argv:
        data = build()
        verify(data)
        (HERE / 'reductions.json').write_text(json.dumps(data, indent=2) + '\n')
    data = json.loads((HERE / 'reductions.json').read_text())
    result = verify(data)
    result['certificate_sha256'] = hashlib.sha256((HERE / 'reductions.json').read_bytes()).hexdigest()
    print(json.dumps(result, indent=2))
