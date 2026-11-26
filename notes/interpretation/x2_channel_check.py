"""Killing-normalized X2 trace in the (01)|(23) color channel.

For an alternating cochain X, exchanging slots 1 and 2 swaps u and v and
negates w. Its uv color tensor is therefore antisymmetric in slots 1,2.
The probe of flagship.tex is varpi=d(<X,kk0>-2<X,ff0>).
In this channel kk0=d*P_singlet, ff0=-P_ad and
P_X2=(kk1-kk2)/2+ff0. The crossing antisymmetry proves
varpi(X)=-2d*tr(P_X2 X^uv). This is the trace of the X2 block.
The complete F223^uv has a nonzero singlet component.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.dont_write_bytecode = True
source = Path(__file__).resolve().parents[2] / 'deepen' / 'simple-all'
sys.path.insert(0, str(source / 'code'))
import lie
import t4


def exact_matrix(matrix):
    coo = matrix.tocoo()
    entries = {(int(i), int(j)): int(v) for i, j, v in zip(coo.row, coo.col, coo.data)}
    return sp.SparseMatrix(*coo.shape, entries)


def matrix_example(key):
    algebra = lie.build(key)
    d = algebra.dim
    ratio, valid = algebra.killing_ratio()
    assert valid
    ratio = sp.Rational(ratio)
    metric = ratio * exact_matrix(algebra.kap)
    dual_pair = sp.kronecker_product(metric.inv(),metric.inv())
    bracket = exact_matrix(algebra.FU).T
    adjoint = dual_pair * bracket.T * metric
    assert bracket * adjoint == -sp.eye(d)
    swap = sp.SparseMatrix(d*d,d*d,{(i*d+j,j*d+i):1 for i in range(d) for j in range(d)})
    projector_ad = -adjoint*bracket
    projector_x = (sp.eye(d*d)-swap)/2-projector_ad
    assert projector_x**2 == projector_x and bracket*projector_x == sp.zeros(d,d*d)
    assert projector_x.trace() == d*(d-3)//2

    engine = t4.Engine(algebra)
    def operator(tensor, scaling):
        # Convert the stored four-lower tensor to Killing normalization and raise two slots.
        return (dual_pair*exact_matrix(tensor.mat())*ratio**scaling
                / sp.Integer(engine.s)**tensor.sp)

    singlet = operator(engine.kk,2)/d
    assert singlet**2 == singlet and singlet.trace() == 1
    assert operator(engine.ff,1) == -projector_ad
    uv = t4.QB.index('uv')
    _, w223 = engine.weights()
    F = operator(t4.alt_bb(w223)[uv],-1)
    beta = operator(t4.beta_t(engine.ff)[uv],1)
    Gkk = operator(t4.alt_bb(engine.kk)[uv],2)
    assert (singlet*F).trace() == -6
    assert projector_x*beta*projector_x == sp.zeros(d*d)
    assert projector_x*Gkk*projector_x == -8*projector_x
    varpi = d*((d*singlet+2*projector_ad)*F).trace()
    assert varpi == -2*d*(projector_x*F).trace()
    print(key,'dim X2',projector_x.trace(),'singlet trace of F', (singlet*F).trace(),
          'X2 trace of F', (projector_x*F).trace(),'varpi',varpi)


def archived_contractions():
    records = json.loads((source/'logs'/'b_data.json').read_text())
    for row in records:
        d = row['dim']
        n4 = Fraction(row['moments']['n4'])
        for name, vector in row['v'].items():
            # For each of six momentum coefficients: kk0,kk1,kk2,ff0,ff1,ff2.
            offset = 6*t4.QB.index('uv')
            kk0,kk1,kk2,ff0 = map(Fraction,vector[offset:offset+4])
            assert (kk1,kk2)==(-kk0,0)
            trace = (kk1-kk2)/2+ff0
            assert d*(kk0-2*ff0) == -2*d*trace
            if name=='F223': assert trace==3*d-8*n4 and kk0==-6*d
            if name=='beta': assert trace==0
            if name=='Gkk': assert trace==-4*d*(d-3)
    print(len(records),'archived simple-algebra contraction tables agree.')


if __name__ == '__main__':
    for key in ('A1','A2','B2'):
        matrix_example(key)
    archived_contractions()
