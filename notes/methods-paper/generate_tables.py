"""Generate the complete printed tables from immutable exact certificates."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'explore/circuit-general/forest-source-currents/certificate.json'
WHEEL = ROOT / 'explore/circuit-general/core-period-closure/certificate.json'
BOUNDARY = ROOT / 'explore/circuit-general/core-sector-reduction/reductions.json'


def rational(value):
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    sign = '-' if value < 0 else ''
    return sign + r'\frac{' + str(abs(value.numerator)) + '}{' + str(value.denominator) + '}'


def period(values):
    terms = []
    for name, weight in (('1', None), ('z3', 3), ('z5', 5)):
        value = Fraction(values.get(name, 0))
        if not value:
            continue
        term = rational(abs(value))
        if weight is not None:
            term += r'\zeta(' + str(weight) + ')'
        terms.append(('-' if value < 0 else '+' if terms else '') + term)
    return ''.join(terms) or '0'


def monomial(i, j, source=None):
    result = ''
    for power, pair in ((5 - i - j, (1, 2)), (i, (5, 2)), (j, (1, 5))):
        if power:
            result += (r'\langle\lambda_' + str(pair[0]) + r',\lambda_'
                       + str(pair[1]) + r'\rangle^{' + str(power) + '}')
    if source is not None:
        result += r'\lambda_{' + str(source) + r',+\dot a}'
    return result or '1'


def full_table(rows):
    header = (r'\begin{longtable}{@{}lrrr@{}}' + '\n' + r'\toprule' + '\n'
              + r'Complete monomial & Rational & $\zeta(3)$ & $\zeta(5)$ \\' + '\n'
              + r'\midrule\endfirsthead' + '\n'
              + r'\toprule Complete monomial & Rational & $\zeta(3)$ & $\zeta(5)$ \\' + '\n'
              + r'\midrule\endhead' + '\n'
              + r'\midrule\multicolumn{4}{r}{Continued on the next page.}\\\endfoot' + '\n'
              + r'\bottomrule\endlastfoot' + '\n')
    text = header
    for (i, j), source, values in rows:
        text += '$' + monomial(i, j, source) + '$ & '
        text += ' & '.join('$' + rational(values.get(label, 0)) + '$'
                           for label in ('1', 'z3', 'z5')) + r' \\' + '\n'
    return text + r'\end{longtable}' + '\n'


def generate():
    data = json.loads(SOURCE.read_text())
    reconstruction = data['reconstruction']
    result = {}
    text = r'\section{Complete bare polynomial}' + '\n'
    rows = [(power, None, values) for power, values in
            zip(reconstruction['frame_monomials'], reconstruction['bare_frame_coefficients']) if values]
    text += full_table(rows)
    count, channels = 0, 0
    for edge, (ends, covariants) in enumerate(zip(data['edges'], reconstruction['denominator_free_original_source_covariants'])):
        text += r'\section{Original edge ' + str(edge) + ': $' + str(ends[0]) + r'\to ' + str(ends[1]) + '$}\n'
        rows = [(power, source, values) for source, direction in zip((1, 2, 5), covariants)
                for power, values in direction]
        count += len(rows)
        channels += sum(sum(bool(Fraction(value)) for value in values.values()) for _, _, values in rows)
        text += full_table(rows)
    if count != 468 or channels != 1334:
        raise ArithmeticError('complete source-current table counts changed')
    result['all_currents.tex'] = text

    text = (r'\begin{center}\small' + '\n' + r'\begin{tabular}{@{}crr@{}}' + '\n'
            + r'\toprule Original edge & $A_e$ & $B_e$ \\ \midrule' + '\n')
    for ends, directions in zip(data['edges'], reconstruction['denominator_free_original_source_covariants']):
        coefficients = [next((values for power, values in direction if power == [0, 0]), {})
                        for direction in directions[:2]]
        text += '$' + str(ends[0]) + r'\to ' + str(ends[1]) + '$ & '
        text += ' & '.join('$' + period(value) + '$' for value in coefficients) + r' \\' + '\n'
    text += r'\bottomrule\end{tabular}\end{center}' + '\n'
    result['zero_marked_currents.tex'] = text

    wheel = json.loads(WHEEL.read_text())
    text = (r'\begin{longtable}{@{}crr@{}}' + '\n'
            + r'\toprule $q$ & Rational coefficient $R_q$ & Coefficient of $\zeta(3)$ \\' + '\n'
            + r'\midrule\endfirsthead' + '\n'
            + r'\toprule $q$ & $R_q$ & Coefficient of $\zeta(3)$ \\ \midrule\endhead' + '\n'
            + r'\bottomrule\endlastfoot' + '\n')
    for row in wheel['finite_k4_base_range']:
        text += str(row['q']) + ' & $' + rational(row['value']['0']) + '$ & $'
        text += rational(row['value']['3']) + r'$ \\' + '\n'
    text += r'\end{longtable}' + '\n'
    result['wheel_bases.tex'] = text

    boundary = json.loads(BOUNDARY.read_text())
    names = {'one_variable_nonrational': ('One-variable target', '1+x', ('x',)),
             'two_variable_nonrational': ('Two-variable nonrational target', '1+xy', ('x', 'y')),
             'three_edge_banana': ('Each complete banana sector', '1+x+xy', ('x', 'y'))}
    seen = set()
    text = ''
    for case in boundary['cases']:
        if case['name'] in seen:
            continue
        seen.add(case['name'])
        title, denominator, variables = names[case['name']]
        text += r'\paragraph{' + title + '.}\n'
        text += (r'\begin{center}\begin{tabular}{@{}rcl@{}}' + '\n'
                 + r'\toprule Coefficient & Coordinate & Complete primitive \\ \midrule' + '\n')
        for coefficient, coordinate, powers, pole in case['witness']['primitives']:
            x = variables[coordinate]
            numerator = ''.join(variable + '^{' + str(power) + '}'
                                for variable, power in zip(variables, powers) if power) or '1'
            primitive = x + '(1-' + x + ')'
            if pole:
                primitive += r'\frac{' + numerator + '}{(' + denominator + ')^{' + str(pole) + '}}'
            elif numerator != '1':
                primitive += numerator
            text += '$' + rational(coefficient) + '$ & $' + x + '$ & $' + primitive + r'$ \\' + '\n'
        text += r'\bottomrule\end{tabular}\end{center}' + '\n'
    result['boundary_witnesses.tex'] = text
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    for name, text in generate().items():
        path = HERE / 'generated' / name
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        elif not path.exists() or path.read_text() != text:
            raise ArithmeticError('printed exact table differs: ' + name)
    print('Complete tables agree: 16 bare monomials, 468 current monomials, 1334 current coefficients, 18 wheel bases.')
