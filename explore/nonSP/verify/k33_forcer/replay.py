"""From the K3,3 Kirchhoff density to the exact four-source coefficient.

Requires the dependencies of ../k33_positive/check.py, FORM and FORCER.
All FORM jobs run serially in one temporary directory.
"""
import argparse
import re
import shutil
import subprocess
import sys
from fractions import Fraction
from itertools import permutations
from math import factorial, prod
from pathlib import Path
from tempfile import TemporaryDirectory

import sympy as sp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'k33_positive'))
from check import AUTOMORPHISMS, INJECTIONS, K33, grid_input, numerator
from vk import G, density_exact

EP = sp.Symbol('ep')
FIELD = sp.QQ.old_frac_field(EP)
# Ordered so that the five ordinary graph integrals reduce triangularly.
MISSING = ((), (5, 8), (5, 6), (1, 3), (1, 5))
MASTERS = ('BasicNOIntegral', 'GschemeConstants(0,0)*BasicT1Integral',
           'GschemeConstants(0,0)^2*GschemeConstants(1,0)',
           'GschemeConstants(0,0)^2*GschemeConstants(2,0)',
           'GschemeConstants(0,0)*GschemeConstants(1,0)*GschemeConstants(2,0)')


def vertices():
    incident = [[] for _ in range(6)]
    for i, (v, w) in enumerate(K33):
        name = f'p{i}' if i else 'Q'
        incident[v].append(name)
        incident[w].append('-' + name)
    return '*'.join('vx(' + ','.join(edges) + ')' for edges in incident)


def gram():
    # Edges 1,...,5 form a spanning tree after deleting edge 0.
    momenta = ['Q', 'p6', 'p7', 'p8']
    terms = []
    for p in permutations(range(4)):
        sign = sum(p[i] > p[j] for j in range(4) for i in range(j))
        monomial = '*'.join(f'{momenta[i]}.{momenta[j]}' for i, j in enumerate(p))
        terms.append(('-' if sign % 2 else '+') + monomial)
    return '(' + '\n'.join(terms) + ')'


def propagators(powers):
    return '*'.join(f'1/(p{i}.p{i})^{power}'
                    for i, power in enumerate(powers, 1) if power)


def graph_terms(injections=INJECTIONS, source_factor=6):
    trees, records = grid_input(injections)
    uniform = len(trees)**6 * density_exact(G(6, K33), [1]*9, injections[:-1])
    raw = numerator(records, uniform, symmetrize=False)
    weights = sum(raw.transpose(p) for p in AUTOMORPHISMS)
    subsets = [tuple(j for j in range(9) if mask >> j & 1) for mask in range(1, 511)]
    cycles = [min(sum(tree[j] for j in edges) for tree in trees) for edges in subsets]
    seen = set()
    terms = [[], []]
    for row in zip(*weights.nonzero()):
        exponent = tuple(map(int, row))
        if exponent in seen:
            continue
        orbit = {tuple(exponent[e] for e in p) for p in AUTOMORPHISMS}
        seen.update(orbit)
        coefficient = Fraction(len(orbit) * int(weights[row]) * prod(map(factorial, exponent)),
                               len(AUTOMORPHISMS) * source_factor * factorial(5))
        assert all(sum(exponent[e] + 1 for e in edges) > 6*loops
                   for edges, loops in zip(subsets, cycles))
        for out, order in zip(terms, (min, max)):
            cut = order(orbit)
            out.append(f'+({coefficient})*' + propagators([k + 1 for k in cut[1:]]))
    if injections == INJECTIONS:
        assert len(terms[0]) == 175 and len(seen) == 10980
    # Integrability of the original density also controls the dimension-raising
    # subdivision identity. Include orbits cancelled by the automorphism sum.
    for row in zip(*raw.nonzero()):
        exponent = tuple(map(int, row))
        if exponent not in seen:
            assert all(sum(exponent[e] + 1 for e in edges) > 6*loops
                       for edges, loops in zip(subsets, cycles))
            seen.update(tuple(exponent[e] for e in p) for p in AUTOMORPHISMS)
    return ['(' + '\n'.join(column) + ')' for column in terms]


def run_form(directory, executable, name, expression, options=''):
    job = directory / (name + '.frm')
    job.write_text('#-\n#include- forcer.h\nL F = ' + vertices() + '*'
                   + expression + ';\n#call Forcer(' + options + ')\n'
                   + f'#write <{name}.out> "%E",F\n.end\n')
    log = directory / (name + '.log')
    with log.open('w') as stream:
        status = subprocess.run([executable, '-q', job.name], cwd=directory,
                                stdout=stream, stderr=subprocess.STDOUT)
    if status.returncode:
        raise RuntimeError(log.read_text()[-2500:])
    return ''.join((directory / (name + '.out')).read_text().split()).replace('\\', '')


def read_masters(result, shift):
    # FORCER writes one coefficient rat(integer polynomial,integer polynomial) per master.
    for i, master in enumerate(MASTERS):
        result = result.replace(master, f'M{i}')
    parts = re.fullmatch(r'(M[0-4]\*rat\([-+ep^*0-9]+,[-+ep^*0-9]+\))'
                         r'(\+M[0-4]\*rat\([-+ep^*0-9]+,[-+ep^*0-9]+\))*', result)
    if not parts:
        raise ValueError(result)
    vector = [FIELD.zero] * 5
    for column, top, bottom in re.findall(r'M([0-4])\*rat\(([^,]+),([^()]+)\)', result):
        def coefficients(polynomial):
            p = sp.Poly(sp.sympify(polynomial.replace('^', '**'), locals={'ep': EP}),
                        EP, domain=sp.QQ)
            return p.shift(-shift).all_coeffs()
        vector[int(column)] += FIELD((coefficients(top), coefficients(bottom)))
    return vector


def ordinary_coordinates(vector, rows):
    out = []
    for k, row in enumerate(rows):
        assert all(not coefficient for coefficient in row[:k]) and row[k]
        scale = vector[k] / row[k]
        out.append(scale)
        vector = [entry - scale * basis for entry, basis in zip(vector, row)]
    assert not any(vector)
    return out


def laurent(coefficient, last):
    if not coefficient:
        return {}
    top = coefficient.numer().to_list()[::-1]
    bottom = coefficient.denom().to_list()[::-1]
    first_top = next(i for i, a in enumerate(top) if a)
    first_bottom = next(i for i, a in enumerate(bottom) if a)
    first = first_top - first_bottom
    top, bottom = top[first_top:], bottom[first_bottom:]
    series = []
    for i in range(last - first + 1):
        n = top[i] if i < len(top) else 0
        n -= sum(bottom[j] * series[i-j] for j in range(1, min(i+1, len(bottom))))
        series.append(n / bottom[0])
    return dict(zip(range(first, last+1), series))


def master_tables(directory, executable):
    expressions = [propagators([int(i not in missing) for i in range(1, 9)]) for missing in MISSING]
    ordinary = [run_form(directory, executable, f'basis{i}', expression)
                for i, expression in enumerate(expressions)]
    inserted = [run_form(directory, executable, f'gram{i}', expression + '*' + gram())
                for i, expression in enumerate(expressions)]
    for i, expression in enumerate(expressions):
        expanded = run_form(directory, executable, f'poles{i}', expression, 'msbarexpand=5,keepzeta2')
        series = sp.sympify(expanded.replace('^', '**'), locals={'ep': EP})
        assert sp.expand(EP**3 * series).is_polynomial(EP)
    bases = [[read_masters(result, j) for result in ordinary] for j in range(5)]
    moments = [[read_masters(result, j) for result in inserted] for j in range(4)]

    return expressions, bases, moments


def integrate(directory, executable, tables, expression, cut):
    expressions, bases, moments = tables
    vector = read_masters(run_form(directory, executable, f'cut{cut}', expression), 4)
    for j in (4, 3, 2, 1):
        coords = ordinary_coordinates(vector, bases[j])
        d = FIELD([-2, 4 + 2*(j-1)])
        factor = FIELD(8) / prod((d-FIELD(k) for k in (1, 2, 3)), start=FIELD.one)
        vector = [factor * sum((a * row[k] for a, row in zip(coords, moments[j-1])), FIELD.zero)
                  for k in range(5)]
    coords = ordinary_coordinates(vector, bases[0])
    # A three-loop basis integral has at most a triple pole. Expand one
    # extra power to distinguish the two continuations away from d=12.
    terms = []
    for a, graph in zip(coords, expressions):
        series = laurent(a, 4)
        terms.append('+(' + '+'.join(f'({v})*ep^({k})' for k,v in series.items()) + ')*' + graph)
    output = run_form(directory, executable, f'value{cut}', '(' + '\n'.join(terms) + ')',
                      'msbarexpand=7,keepzeta2')
    # No special marker from the master table, and no unexpanded functions.
    assert re.fullmatch(r'[0-9+\-*/^epz]+', output)
    result = sp.sympify(output.replace('^', '**'), locals={'ep': EP})
    assert result.is_polynomial(EP), result
    return result


def solve(directory, executable):
    tables = master_tables(directory, executable)
    z3, z5 = sp.symbols('z3 z5')
    expected = (-539 - 1414*z3 + 2160*z5) / 14400
    continuations = []
    for cut, expression in enumerate(graph_terms()):
        result = integrate(directory, executable, tables, expression, cut)
        assert sp.expand(result.subs(EP, 0) - expected) == 0
        print(('minimum', 'maximum')[cut], 'cut: poles cancel, scalar =', expected, flush=True)
        continuations.append(result.coeff(EP, 1))
    assert sp.expand(continuations[0] - continuations[1]) != 0


def run(calculate):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--form', required=True, help='FORM executable or command name')
    parser.add_argument('--forcer', required=True, type=Path, help='directory containing forcer.h and forcer/')
    args = parser.parse_args()
    executable = shutil.which(args.form)
    if not executable:
        parser.error('FORM executable missing')
    with TemporaryDirectory() as scratch:
        directory = Path(scratch)
        for name in ('forcer.h', 'forcer'):
            source = (args.forcer / name).resolve(strict=True)
            (directory / name).symlink_to(source, target_is_directory=source.is_dir())
        calculate(directory, str(Path(executable).resolve()))


if __name__ == '__main__':
    run(solve)
