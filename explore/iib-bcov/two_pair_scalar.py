"""Exact restriction of the nine-primary trace to two disjoint collinear pairs."""
from fractions import Fraction
from itertools import permutations
from math import lcm, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from flint import fmpq_mat, fmpz_mat, fmpz_mpoly_ctx

from elementary_scalar import ROOT
from three_loop_soft_jet import DOUBLE_SOFT, SUBSETS, orbit_integers


def chart(s, t, w):
    frame = [tuple(int(i == j) for j in range(5)) for i in range(5)]
    first = [frame[0], (s,0,0,0,0), frame[1], (0,t,0,0,0), *frame[2:], w]
    return [*first, tuple(-sum(column) for column in zip(*first))]


def determinant(rows, orders):
    result = 0
    for order, sign in orders:
        entries = [row[i] for row, i in zip(rows,order)]
        if all(entries):
            result += sign * prod(entries)
    return result


def symbolic_traces(exe):
    ctx = fmpz_mpoly_ctx.get(('s','t','w1','w2','w3','w4','w5'))
    s, t, *w = ctx.gens()
    points = chart(s,t,w)
    orders = [(order, (-1)**sum(order[i]>order[j] for j in range(5) for i in range(j)))
              for order in permutations(range(5))]
    factors = [ctx.constant(0)]
    by_terms = {}
    labels = []
    for subset in SUBSETS:
        minor = ctx.constant(0) + determinant([points[i] for i in subset], orders)
        if not minor:
            labels.append(0)
            continue
        sign = 1 if minor.leading_coefficient()>0 else -1
        minor *= sign
        key = tuple(minor.terms())
        if key not in by_terms:
            by_terms[key] = len(factors)
            factors.append(minor)
        labels.append(sign * by_terms[key])
    data = [(len(DOUBLE_SOFT),1,-1), *DOUBLE_SOFT, labels]
    stdin = '\n'.join(' '.join(map(str,row)) for row in data)
    result = run([exe],input=stdin,text=True,capture_output=True,check=True)
    traces, polynomial = [], ctx.constant(0)
    for line in result.stdout.splitlines():
        if line == 'end':
            traces.append(polynomial)
            polynomial = ctx.constant(0)
        else:
            multiplicity, *indices = map(int,line.split())
            polynomial += multiplicity*prod(factors[i] for i in indices)
    assert len(traces) == len(DOUBLE_SOFT)
    return traces


def restricted_trace(traces):
    equations, periods = [], []
    for name in ('three_loop_constraints.json','three_loop_four_pair.json'):
        data = json.loads((ROOT/'operator'/name).read_text())
        equations.extend(data['rows'])
        periods.extend(data['reduced_periods'])
    reduced, rank = fmpq_mat(equations).rref()
    pivots = [next(i for i in range(18) if reduced[j,i]) for j in range(rank)]
    assert rank == len(equations) == 15
    square = fmpq_mat([[r[i] for i in pivots] for r in equations])
    # Every column identity holds over the full seven-variable polynomial ring.
    projection = square.solve(fmpq_mat(equations))
    periods = square.solve(fmpq_mat(periods))
    def combine(column):
        den = lcm(*(int(x.q) for x in column))
        return sum(int(x*den)*traces[i] for x,i in zip(column,pivots)), den
    for j,trace in enumerate(traces):
        polynomial, den = combine([projection[i,j] for i in range(rank)])
        assert trace*den == polynomial, j
    return [combine([periods[i,j] for i in range(rank)]) for j in range(2)]


def certificate():
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        traces = symbolic_traces(exe)
        channels = restricted_trace(traces)
        fixtures = [(1,1,1,0,1,-1,1), (-2,3,1,-1,2,-3,4), (3,-4,-2,1,3,2,-1)]
        rows = orbit_integers(exe,[chart(s,t,w) for s,t,*w in fixtures])
        for point,row in zip(fixtures,rows):
            assert [int(poly(*point)) for poly in traces] == row
    terms = [p.to_dict() for p in traces]
    powers = sorted(set().union(*terms))
    assert all(e[0]>=2 and e[1]>=2 for e in powers)
    matrix = fmpz_mat([[int(poly.get(e,0)) for e in powers] for poly in terms])
    values = [[str(Fraction(int(poly(*point)),den)) for poly,den in channels]
              for point in fixtures]
    assert len(powers) == 1676 and matrix.rank() == 12
    assert values[0] == ['-4582336574567/5402002814208000','1355079563/2084106024000']
    # A smooth complex zero away from the common minor zero set:
    # in this chart set s=t=1 and w=(x,1,0,0,0).
    complex_slice = [{i:sum(Fraction(int(a),den) for e,a in poly.terms()
                            if e[2]==i and not any(e[4:]))
                      for i in range(int(poly.degrees()[2])+1)} for poly,den in channels]
    assert complex_slice == [
        {0:Fraction(-49188779,2715939072000),1:0,2:Fraction(-7896562031,211843247616000)},
        {0:Fraction(7783,523908000),1:0,2:Fraction(2500061,81729648000)}]
    for i in (0,2):
        assert complex_slice[0][i]+Fraction(1202057,10**6)*complex_slice[1][i]<0
    print('Exact two-pair restriction rank 12, monomials',len(powers))
    print('Rational and zeta(3) coefficients',values)
    return channels


if __name__=='__main__':
    certificate()
