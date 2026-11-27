"""Independently parse every printed coefficient and complete boundary witness.

No table generator or production period/current implementation is imported.
The displayed TeX, not invisible annotations, is the checked input.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path
import re

import sympy as sp

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
GENERATED = HERE / 'generated'


def rational(text):
    text = text.strip()
    fraction = re.fullmatch(r'(-?)\\frac\{(\d+)\}\{(\d+)\}', text)
    if fraction:
        sign, numerator, denominator = fraction.groups()
        return F(int(numerator) * (-1 if sign else 1), int(denominator))
    if not re.fullmatch(r'-?\d+', text):
        raise ArithmeticError('nonrational printed coefficient: ' + text)
    return F(text)


def fields(line):
    if not line.endswith(r'\\'):
        raise ArithmeticError('unterminated printed row')
    answer = []
    for value in line[:-2].split('&'):
        value = value.strip()
        answer.append(value[1:-1] if value.startswith('$') and value.endswith('$') else value)
    return answer


def full_coefficients(data):
    expected = {}
    reconstruction = data['reconstruction']
    for (i, j), values in zip(reconstruction['frame_monomials'], reconstruction['bare_frame_coefficients']):
        if values:
            expected[None, None, 5-i-j, i, j] = tuple(F(values.get(k, 0)) for k in ('1', 'z3', 'z5'))
    for edge, directions in enumerate(reconstruction['denominator_free_original_source_covariants']):
        for source, rows in zip((1, 2, 5), directions):
            for (i, j), values in rows:
                expected[edge, source, 5-i-j, i, j] = tuple(F(values.get(k, 0)) for k in ('1', 'z3', 'z5'))
    found, edge, sections = {}, None, []
    bracket = re.compile(r'\\langle\\lambda_(\d),\\lambda_(\d)\\rangle\^\{(\d+)\}')
    for line in (GENERATED / 'all_currents.tex').read_text().splitlines():
        heading = re.fullmatch(r'\\section\{Original edge (\d+): \$(\d+)\\to (\d+)\$\}', line)
        if heading:
            edge, left, right = map(int, heading.groups())
            if [left, right] != data['edges'][edge]:
                raise ArithmeticError('printed original edge orientation changed')
            sections.append(edge)
        if not line.startswith('$'):
            continue
        monomial, *values = fields(line)
        if len(values) != 3:
            raise ArithmeticError('a printed period channel is missing')
        powers = {(1, 2): 0, (5, 2): 0, (1, 5): 0}
        seen = set()
        for left, right, power in bracket.findall(monomial):
            pair = int(left), int(right)
            if pair not in powers or pair in seen or int(power) <= 0:
                raise ArithmeticError('unexpected source bracket')
            powers[pair] = int(power)
            seen.add(pair)
        remainder = bracket.sub('', monomial)
        if edge is None:
            if remainder:
                raise ArithmeticError('the bare row has an extra factor')
            source = None
        else:
            match = re.fullmatch(r'\\lambda_\{([125]),\+\\dot a\}', remainder)
            if not match:
                raise ArithmeticError('a full free spinor index or original source was lost')
            source = int(match.group(1))
        if sum(powers.values()) != 5:
            raise ArithmeticError('a complete bracket degree changed')
        key = edge, source, powers[1, 2], powers[5, 2], powers[1, 5]
        if key in found:
            raise ArithmeticError('a printed monomial is duplicated')
        found[key] = tuple(map(rational, values))
    if sections != list(range(11)) or found != expected:
        raise ArithmeticError('the complete printed source polynomial differs from its exact certificate')
    bare = sum(key[0] is None for key in found)
    current = len(found)-bare
    channels = sum(sum(bool(x) for x in value) for key, value in found.items() if key[0] is not None)
    if (bare, current, channels) != (16, 468, 1334):
        raise ArithmeticError('the full printed counts changed')
    return {'bare_monomials': bare, 'current_monomials': current,
            'nonzero_current_rational_coefficients': channels,
            'all_printed_rational_entries': len(found)*3}


def period(text):
    if text == '0':
        return F(), F(), F()
    token = re.compile(r'([+-]?)(\\frac\{\d+\}\{\d+\}|\d+)(?:\\zeta\(([35])\))?')
    result = [F(), F(), F()]
    offset = 0
    for match in token.finditer(text):
        if match.start() != offset:
            raise ArithmeticError('unparsed scalar expression in printed specialization')
        sign, value, weight = match.groups()
        result[0 if weight is None else 1 if weight == '3' else 2] += rational(value) * (-1 if sign == '-' else 1)
        offset = match.end()
    if offset != len(text):
        raise ArithmeticError('unparsed tail in printed scalar expression')
    return tuple(result)


def zero_marked(data):
    found = []
    for line in (GENERATED / 'zero_marked_currents.tex').read_text().splitlines():
        if not line.startswith('$'):
            continue
        edge, first, second = fields(line)
        match = re.fullmatch(r'(\d+)\\to (\d+)', edge)
        if not match or list(map(int, match.groups())) != data['edges'][len(found)]:
            raise ArithmeticError('specialized current has wrong original orientation')
        found.append((period(first), period(second)))
    expected = []
    for directions in data['reconstruction']['denominator_free_original_source_covariants']:
        values = []
        for direction in directions[:2]:
            coefficient = next((value for power, value in direction if power == [0, 0]), {})
            values.append(tuple(F(coefficient.get(k, 0)) for k in ('1', 'z3', 'z5')))
        expected.append(tuple(values))
    if found != expected or len(found) != 11:
        raise ArithmeticError('complete vanishing-bare specialization differs')
    for vertex in range(7):
        for source in range(2):
            for channel in range(3):
                divergence = sum((int(left == vertex)-int(right == vertex))*found[e][source][channel]
                                 for e, (left, right) in enumerate(data['edges']))
                if divergence:
                    raise ArithmeticError('printed specialization violates source conservation')
    lower = tuple(value[0] + F(5, 4)*value[1] for value in found[0])
    if lower != (F(653, 1828915200), F(2, 1828915200)):
        raise ArithmeticError('the displayed exact positivity calculation changed')
    return {'complete_specialized_spinors': len(found),
            'coefficientwise_specialized_conservation_checks': 42,
            'strict_positive_lower_bounds': list(map(str, lower))}


def wheel_bases():
    stored = json.loads((ROOT / 'explore/circuit-general/core-period-closure/certificate.json').read_text())
    rows = []
    for line in (GENERATED / 'wheel_bases.tex').read_text().splitlines():
        if not re.match(r'^\d+ &', line):
            continue
        q, first, second = fields(line)
        q, first, second = int(q), rational(first), rational(second)
        if second != F(12*factorial(q-2)**2, factorial(2*q-2)):
            raise ArithmeticError('the printed uniform zeta coefficient differs')
        if (q == 3 and first != -F(1, 6)) or (q >= 4 and first <= 0):
            raise ArithmeticError('a finite wheel rational direction vanishes')
        rows.append({'q': q, 'value': {'0': str(first), '3': str(second)}})
    if rows != stored['finite_k4_base_range']:
        raise ArithmeticError('the complete printed wheel base range differs')
    return {'complete_wheel_base_rows': len(rows)}


def group(text, start=0):
    if text[start] != '{':
        raise ArithmeticError('expected a TeX brace group')
    depth = 0
    for end in range(start, len(text)):
        if text[end] == '{': depth += 1
        if text[end] == '}': depth -= 1
        if depth == 0: return text[start+1:end], end+1
    raise ArithmeticError('unbalanced printed primitive')


def boundary_witnesses():
    x, y = sp.symbols('x y')
    cases = {
        'One-variable target.': (1+x, x**7/(1+x)**5, F(1553, 64)-35/(1+x), 6),
        'Two-variable nonrational target.': (1+x*y, x*x*y*y/(1+x*y)**3, -F(1, 2)+1/(1+x*y)-F(3, 4)*x/(1+x*y), 5),
        'Each complete banana sector.': (1+x+x*y, x/(1+x+x*y)**3, F(1, 12), 5)}
    denominator_text = {'1+x':1+x, '1+xy':1+x*y, '1+x+xy':1+x+x*y}
    checked, faces, count, current, residual = [], 0, 0, None, None
    def finish():
        if current is not None:
            if count != cases[current][3] or sp.cancel(residual) != 0:
                raise ArithmeticError('a complete printed derivative identity failed')
            checked.append(current)
    for line in (GENERATED / 'boundary_witnesses.tex').read_text().splitlines():
        heading = re.fullmatch(r'\\paragraph\{(.+)\}', line)
        if heading:
            finish()
            current = heading.group(1)
            if current not in cases: raise ArithmeticError('unknown scalar target')
            _, target, masters, _ = cases[current]
            residual, count = target-masters, 0
        if not line.startswith('$'):
            continue
        coefficient, coordinate, primitive = fields(line)
        if coordinate not in ('x', 'y'):
            raise ArithmeticError('unrecognized differentiation coordinate')
        variable = x if coordinate == 'x' else y
        prefix = coordinate+'(1-'+coordinate+')'
        if not primitive.startswith(prefix):
            raise ArithmeticError('a complete vanishing-face factor is absent')
        remaining = primitive[len(prefix):]
        expression = variable*(1-variable)
        pole = 0
        if remaining.startswith(r'\frac'):
            numerator, end = group(remaining, 5)
            denominator, finish_at = group(remaining, end)
            if finish_at != len(remaining): raise ArithmeticError('extra printed primitive factor')
            match = re.fullmatch(r'\((1\+x|1\+xy|1\+x\+xy)\)\^\{(\d+)\}', denominator)
            if not match: raise ArithmeticError('unexpected primitive denominator')
            base, pole = match.groups()
            if denominator_text[base] != cases[current][0]: raise ArithmeticError('wrong target denominator')
            expression /= denominator_text[base]**int(pole)
        else:
            numerator = remaining or '1'
        if numerator != '1':
            tokens = re.findall(r'([xy])\^\{(\d+)\}', numerator)
            if ''.join(a+'^{'+b+'}' for a,b in tokens) != numerator:
                raise ArithmeticError('unparsed primitive monomial')
            for name, exponent in tokens:
                expression *= (x if name == 'x' else y)**int(exponent)
        for endpoint in (0, 1):
            if sp.cancel(expression.subs(variable, endpoint)) != 0:
                raise ArithmeticError('a printed primitive leaves a face term')
            faces += 1
        residual -= sp.Rational(rational(coefficient))*sp.diff(expression, variable)
        count += 1
    finish()
    if checked != list(cases): raise ArithmeticError('a complete printed scalar target is missing')
    return {'complete_printed_boundary_identities': len(checked),
            'all_printed_primitives': faces//2, 'printed_primitive_faces': faces}


def verify():
    source_path = ROOT / 'explore/circuit-general/forest-source-currents/certificate.json'
    data = json.loads(source_path.read_text())
    result = {'schema':'methods-paper-printed-data.v1', **full_coefficients(data),
              **zero_marked(data), **wheel_bases(), **boundary_witnesses()}
    result['generated_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(GENERATED.glob('*.tex'))}
    result['source_certificate_sha256'] = hashlib.sha256(source_path.read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = verify()
    path = HERE / 'verification/printed_tables.json'
    if args.write: path.write_text(json.dumps(result,indent=2)+'\n')
    elif json.loads(path.read_text()) != result: raise ArithmeticError('printed-table verification record changed')
    print(json.dumps(result,indent=2))
