"""Three-loop prism cuts: transport the five ordinary graph integrals to d=12."""
from collections import defaultdict
from functools import lru_cache
from itertools import permutations
from math import prod
from pathlib import Path
import re
import sys

import networkx as nx
import sympy as sp
from flint import fmpq_poly

from periods import PRISM, form, read_vectors, tagged

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'k33_forcer'))
from replay import EP, FIELD, K33, laurent, ordinary_coordinates, propagators
import replay as k33

G = [f'GschemeConstants({i},0)' for i in range(3)]
MASTERS = (G[0]+'*BasicT1Integral', G[0]+'*'+G[1]+'*'+G[2], G[0]+'^2*'+G[2],
           G[0]+'^2*'+G[1], G[0]+'^3')
MISSING = ((6,7), (1,3), (1,4), (1,8), (5,8))
CHORDS = (5,7,8)


def gram(edges):
    chords = CHORDS if edges == PRISM else (6,7,8)
    assert nx.is_tree(nx.Graph(edges[e] for e in range(1,9) if e not in chords))
    momenta = ['Q'] + [f'p{i}' for i in chords]
    terms = []
    for p in permutations(range(4)):
        inversions = sum(x > y for i,x in enumerate(p) for y in p[i+1:])
        terms.append(('-' if inversions%2 else '+') + '*'.join(
            f'{a}.{momenta[j]}' for a,j in zip(momenta,p)))
    return '('+''.join(terms)+')'


@lru_cache(maxsize=2048)
def shifted_polynomial(source, shift):
    assert re.fullmatch(r'[+\-0-9*ep^]+', source)
    pattern = r'([+-]?)(?:(\d+)\*)?ep(?:\^(\d+))?|([+-]?\d+)'
    matches = list(re.finditer(pattern,source))
    assert ''.join(m.group() for m in matches)==source
    terms = [m.groups() for m in matches]
    coefficients = {}
    for sign, c, exponent, constant in terms:
        k = 0 if constant else int(exponent or 1)
        assert k not in coefficients
        coefficients[k] = int(constant) if constant else (-1 if sign=='-' else 1)*int(c or 1)
    polynomial = fmpq_poly([coefficients.get(i,0) for i in range(max(coefficients)+1)])
    return polynomial(fmpq_poly([-shift, 1]))


def rational(pair, shift):
    top, bottom = (shifted_polynomial(s,shift) for s in pair)
    return FIELD((list(reversed(top.coeffs())) or [0], list(reversed(bottom.coeffs()))))


def vector(row, shift, edges):
    masters = MASTERS if edges==PRISM else k33.MASTERS
    assert row.keys() <= set(masters)
    return [rational(row[k],shift) if k in row else FIELD.zero for k in masters]


def dimension_shift(folder, executable, shifts=4, edges=PRISM):
    assert edges==PRISM or edges==K33
    missing = MISSING if edges==PRISM else k33.MISSING
    ordinary = [propagators([int(i not in m) for i in range(1,9)]) for m in missing]
    basis = read_vectors(form(folder, executable, 'basis', tagged(enumerate(ordinary)),edges=edges))
    inserted = [a+'*'+gram(edges) for a in ordinary]
    moments = read_vectors(form(folder, executable, 'gram', tagged(enumerate(inserted)),edges=edges))
    # Transfer complete master vectors, with no expansion between dimension shifts.
    transfer = [[FIELD(int(i==j)) for j in range(5)] for i in range(5)]
    for j in range(shifts,0,-1):
        basis_j = [vector(basis[i],j,edges) for i in range(5)]
        gram_j = [vector(moments[i],j-1,edges) for i in range(5)]
        d = FIELD([-2, 4+2*(j-1)])
        factor = FIELD(8) / prod((d-FIELD(i) for i in (1,2,3)),start=FIELD.one)
        for i, row in enumerate(transfer):
            coords = ordinary_coordinates(row, basis_j)
            transfer[i] = [factor*sum((v*b[k] for v,b in zip(coords,gram_j)),FIELD.zero)
                           for k in range(5)]
    transfer = [ordinary_coordinates(row,[vector(basis[i],0,edges) for i in range(5)])
                for row in transfer]
    return ordinary, transfer


def master_series(folder, executable, last, shifts=4, edges=PRISM):
    ordinary, transfer = dimension_shift(folder, executable, shifts, edges)
    # FORM's msbarexpand cutoff is the number of powers retained above the
    # default three-loop pole. The ordinary basis has poles of order at most three.
    def valuation(polynomial):
        return next(i for i,c in enumerate(reversed(polynomial.to_list())) if c)
    lower = min(valuation(c.numer())-valuation(c.denom()) for row in transfer for c in row if c)
    need = max(0, last-lower)
    series = []
    for i, expr in enumerate(ordinary):
        text = ''.join(form(folder,executable,f'ordinary{i}',expr,
                            option=f'msbarexpand={need+3},keepzeta2',edges=edges).split()).replace('\\','')
        assert re.fullmatch(r'[0-9+\-*/^epz]+',text),text
        expanded = sp.sympify(text.replace('^','**'), locals={'ep':EP}).expand()
        assert (expanded*EP**3).expand().is_polynomial(EP)
        assert expanded.coeff(EP,need) != 0
        series.append({k:expanded.coeff(EP,k) for k in range(-3,need+1)})
    output = []
    for row in transfer:
        result = defaultdict(int)
        for c, basis in zip(row,series):
            for j,coef in laurent(c,last+3).items():
                for k,value in basis.items():
                    if j+k <= last:
                        result[j+k] += sp.Rational(coef)*value
        output.append({k:sp.expand(v) for k,v in result.items() if sp.expand(v)})
    return output


def period_evaluator(folder, executable, last, shifts=4, edges=PRISM):
    master = master_series(folder, executable, last, shifts, edges)
    symbols = sorted(set.union(*(set().union(*(v.free_symbols for v in row.values()))
                                  for row in master)),key=str)
    data = [{k:sp.Poly(v,*symbols).terms() for k,v in row.items()} for row in master]

    def integrate(row):
        result = defaultdict(int)
        for i, c in enumerate(vector(row,shifts,edges)):
            if not c:
                continue
            basis = data[i]
            if not basis:
                continue
            series = laurent(c, -min(basis))
            assert not series or min(series)>=-last
            for j,coef in series.items():
                for k, terms in basis.items():
                    if j+k<=0:
                        for ex,v in terms:
                            result[j+k,ex] += coef*v
        assert all(power==0 or coeff==0 for (power,ex),coeff in result.items())
        return sp.Add(*(sp.Rational(v)*prod(z**k for z,k in zip(symbols,ex))
                        for (power,ex),v in result.items() if v))
    return integrate
