"""Render every exact charge coefficient with the original free spinor index."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import json

import quartic_b_squared as calculation

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'FULL_SOURCE_TENSOR.md'
VECTORS = (r'\lambda_{\gamma_4,+\dot a}', r'\lambda_{c,+\dot a}',
           r'\lambda_{\gamma_5,+\dot a}')
BRACKETS = (r'\langle\lambda_{\gamma_4},\lambda_c\rangle',
            r'\langle\lambda_{\gamma_5},\lambda_c\rangle',
            r'\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle')


def monomial(powers):
    return ''.join(value if power == 1 else value + '^{%d}' % power
                   for value, power in zip(BRACKETS, powers) if power)


def scalar(value):
    value = abs(Fraction(value))
    if value == 1:
        return ''
    if value.denominator == 1:
        return str(value.numerator)
    return r'\frac{%d}{%d}' % (value.numerator, value.denominator)


def spinor(terms):
    grouped = defaultdict(list)
    for vector, powers, coefficient in terms:
        if Fraction(coefficient):
            grouped[vector].append((tuple(powers), Fraction(coefficient)))
    if not grouped:
        return '\\[0.\\]\n'
    result = ['\\[', '\\begin{aligned}']
    for index, (vector, polynomial) in enumerate(sorted(grouped.items())):
        result.append(('&' if index == 0 else r'{}+{}&') + VECTORS[vector]
                      + r'\left[\begin{aligned}')
        for offset in range(0, len(polynomial), 3):
            line = '&'
            for position, (powers, coefficient) in enumerate(polynomial[offset:offset + 3], offset):
                sign = '-' if coefficient < 0 else ('' if position == 0 else '+')
                line += sign + scalar(coefficient) + monomial(powers)
            if offset + 3 < len(polynomial):
                line += r'\\'
            result.append(line)
        result.append(r'\end{aligned}\right]' + (r'\\' if index + 1 < len(grouped) else '.'))
    result += ['\\end{aligned}', '\\]', '']
    return '\n'.join(result)


def charge_word(word):
    return ''.join('q_%d' % flavor + ('^{%d}' % count if count > 1 else '')
                   for flavor, count in sorted(Counter(word).items()))


def render(data):
    result = [r'''# Complete four-loop source tensor

The interaction is

\[
-igc\sum_{i=1}^{4}q_i\gamma_i\beta^i
+\gamma_1\gamma_2\gamma_3\gamma_4,
\]

and the inputs are

\[
\beta^1(\partial_{+\dot a}\beta^2)\gamma_5,
\qquad b^2\beta^3.
\]

Here \(q_i\) is the charge of the indicated interacting matter field.
The field \(\gamma_5\) is a separate free neutral spectator, absent from
the interaction. The output is \(c\gamma_4\gamma_5\), with derivatives
specified by the source polynomial below. Each source \(\lambda\)
belongs to the field in its subscript; the root balances the three
sources. The bracket is the manuscript's antisymmetric spinor pairing,
and \(\dot a\) is the original free input index.

Every displayed coefficient is the complete signed geometric sum.
To obtain the coefficient of \(\hbar^4\), multiply the charge monomial
by its displayed spinor, sum the charge monomials, replace each source
by the derivative on its indicated output field, and multiply by

\[
\frac{ig^3}{1536(8\pi^2)^4}.
\]

This factor already contains the manuscript substitution
\(\lambda=-\partial/2\), the differentiated head edge, and the four
superpotential placements divided by \(4!\). The input is the actual
monomial \(b^2\), without division by \(2!\).

## Complete physical specialization

For the gauge-invariant, anomaly-free charge assignment
\((q_1,q_2,q_3,q_4)=(1,-1,0,0)\), the complete spinor is
''']
    result.append(spinor(data['complete_physical_covariant_spinor']['1']))
    result.append(r'''
Setting the free spectator to one sets \(\lambda_{\gamma_5}=0\), and
the whole fourth-order coefficient vanishes. The source restriction
\(\lambda_{\gamma_5}=\lambda_{\gamma_4}\), which does not identify the
two physical fields, gives the exact shorter expression

\[
\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^4}{9600}
\bigl(10\lambda_{\gamma_4,+\dot a}-3\lambda_{c,+\dot a}\bigr).
\]

## Every symbolic charge coefficient

The following coefficients are computed before imposing
\(q_1+q_2+q_3+q_4=0\). Each is rational separately.
''')
    coefficients = {tuple(word): channels for word, channels in data['complete_covariant_charge_polynomial']}
    for word in combinations_with_replacement(range(1, 5), 3):
        result.append('\n### \\(' + charge_word(word) + '\\)\n\n')
        result.append(spinor(coefficients.get(word, {}).get('1', ())))
    result.append(r'''
Every \(\zeta(3)\) coefficient is zero. The source reconstruction also
checks the \(\zeta(5)\) channel and every other period label present in
the exact graph evaluations. No omitted scalar integral contributes.

Generated directly from [the complete certificate](four_loop_certificate.json)
by [render_tensor.py](render_tensor.py); `python render_tensor.py --check`
checks every printed coefficient against the stored exact tensor.
''')
    return ''.join(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    options = parser.parse_args()
    data = json.loads(calculation.CERTIFICATE.read_text())
    value = render(data)
    if options.check:
        if OUTPUT.read_text() != value:
            raise ArithmeticError('a printed full-spinor coefficient differs')
        print('Every printed physical and symbolic-charge spinor coefficient agrees.')
    else:
        OUTPUT.write_text(value, encoding='utf-8')
        print('Wrote every complete full-source tensor coefficient.')


if __name__ == '__main__':
    main()
