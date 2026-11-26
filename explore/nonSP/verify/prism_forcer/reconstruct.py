"""Reconstruct and verify the full four-loop prism kernel for conserved sources."""
import argparse
from fractions import Fraction
from math import factorial, prod
import json
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

import sympy as sp

from periods import BATCH, CUTS, HERE, PRISM, monomial_jobs, read_vectors, support
from shift import period_evaluator
from sources import PAIRS, brackets, column, sampling, weights

Z3, Z5 = sp.symbols('z3 z5')


def period_bank(folder, executable, output):
    evaluate = period_evaluator(folder, executable, 4)
    result = []
    for start in range(0, len(support()), BATCH):
        rows = [read_vectors((folder/f'period{start}_{cut}.out').read_text()) for cut in CUTS]
        indices = range(start, min(start+BATCH,len(support())))
        assert all(set(row)==set(indices) for row in rows)
        for i in indices:
            a,b = (evaluate(row[i]) for row in rows)
            assert a==b, (i,a,b)
            coeffs = [a.subs({Z3:0,Z5:0}),a.coeff(Z3),a.coeff(Z5)]
            assert all(c.is_Rational for c in coeffs) and a==coeffs[0]+coeffs[1]*Z3+coeffs[2]*Z5
            normalization = Fraction(prod(map(factorial,support()[i])),factorial(5))
            result.append([str(c*normalization) for c in coeffs])
        print('normalized periods',len(result),flush=True)
    output.write_text(json.dumps({'exponents':support(),'coefficients':result},separators=(',',':'))+'\n')


def source_data(folder):
    binary = folder/'grid'
    subprocess.run([os.getenv('CXX','c++'),'-O3','-std=c++17',str(HERE.parent/'k33_positive/grid.cpp'),
                    '-o',str(binary)],check=True)
    _,probes,_ = sampling()
    return [weights(source,binary) for source in probes]


def restricted_checks(seeds, coefficients):
    from fractions import Fraction as F
    import sys
    sys.path.insert(0,str(HERE.parent))
    from prism_rung import projection
    ref = json.loads((HERE.parents[1]/'results_nonSP_C2.json').read_text())['prism_4src']['table']
    def evaluate(probes):
        tables = brackets(probes)
        return sp.Matrix.hstack(*(sp.Matrix(column(tables,seed).tolist()) for seed in seeds))*coefficients
    # A quartic in the three brackets at four marks is fixed on [01]=1 by
    # fifteen triangular nodes (and hence on the dense chart [01] != 0).
    probes = []
    expected = []
    for x in range(5):
        for y in range(5-x):
            probes.append([(1,0),(0,1),(x,y),(-1-x,-1-y),(0,0),(0,0)])
            v = [0,0,0]
            for name,pair in ref.items():
                i,j,k = map(int,name[-3:])
                for u,c in enumerate(pair):
                    v[u] += F(c)*y**j*(-x)**k
            expected.append(v)
    assert evaluate(probes)==sp.Matrix(expected)
    pairings = (((0,1),(2,3)),((0,2),(1,3)),((0,1),(2,4)),((0,5),(1,2)))
    answer = []
    for pair in pairings:
        terms = list(projection(pair))
        matrix = evaluate([[source.get(i,(0,0)) for i in range(6)] for source,_ in terms])
        answer.append(list(sp.Matrix([[weight for _,weight in terms]])*matrix))
    assert answer == [[-F(13,28800),F(1,3200),0],[-F(13,57600),F(1,6400),0],
                      [-F(1,28800),F(1,9600),0],[-F(237,6400),-F(71,720),F(3,20)]]


def reconstruct(weight_rows, bank, output):
    seeds, probes, matrix = sampling()
    assert bank['exponents']==list(map(list,support()))
    periods = sp.Matrix([[sp.Rational(c) for c in row] for row in bank['coefficients']])
    values = sp.Matrix(weight_rows)*periods
    coefficients = matrix.inv()*values[:40,:]
    extra = brackets([probes[-1]])
    test = sp.Matrix([[int(column(extra,seed)[0]) for seed in seeds]])
    assert all(test*matrix.inv()) and test*coefficients==values[40:,:]
    restricted_checks(seeds,coefficients)
    print('held-out sources',probes[-1],'coefficients',list(values[40,:]),flush=True)
    record = {'edges':PRISM, 'basis':'edge-sign weighted sum over all twelve prism automorphisms',
              'periods':['1','zeta(3)','zeta(5)'],
              'rows':[[[PAIRS[i] for i in seed],list(map(str,coefficients[j,:]))]
                      for j,seed in enumerate(seeds)],
              'holdout':[probes[-1],list(map(str,values[40,:]))]}
    output.write_text(json.dumps(record,indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--form',required=True)
    parser.add_argument('--forcer',required=True,type=Path)
    args = parser.parse_args()
    executable = shutil.which(args.form)
    if not executable:
        parser.error('FORM executable missing')
    with TemporaryDirectory(prefix='prism-c2-') as temp:
        folder = Path(temp)
        for name in ('forcer','forcer.h'):
            (folder/name).symlink_to((args.forcer/name).resolve(strict=True))
        monomial_jobs(folder,executable)
        period_bank(folder,executable,folder/'periods.json')
        reconstruct(source_data(folder),json.loads((folder/'periods.json').read_text()),folder/'kernel.json')
        for name in ('periods.json','kernel.json'):
            assert (folder/name).read_bytes()==(HERE/name).read_bytes(),name
    print('Both stored tables reproduced exactly.')


if __name__=='__main__':
    main()
