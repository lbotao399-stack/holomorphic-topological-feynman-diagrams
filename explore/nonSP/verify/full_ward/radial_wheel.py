"""Evaluate wheel source periods using its massless radial integrals in D=4 and D=6."""
from collections import Counter
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import numpy as np
from flint import fmpq_mat, fmpz_mat

from contractions import HERE, types
from sampling import pivots
from wheel_faces import DUAL, face_data, star_value
from wheel_ibp import WheelMoments, compile_grid, kernel, wheel_coefficients


PROBES=(
    ((0,0),(-1,1),(0,0),(0,1),(-1,0),(2,-2)),
    ((1,1),(1,1),(0,0),(-1,-1),(0,-1),(-1,0)),
)


def radial(q):
    # m+q-1 is the shifted angular harmonic. In 4D its spectral multiplicity
    # factor is m+1, in 6D it is ((m+2)^3-(m+2))/6.
    if q==2:
        return F(0),F(0),F(20)
    if q==3:
        return F(0),F(5,3),-F(5,3)
    raise ValueError('The radial reconstruction uses q=2 and q=3.')


def lift(moments,exponents,q):
    # Multiply numerator and denominator by U until the denominator is U^6.
    terms=Counter({tuple(exponents):1})
    for _ in range(q,6):
        new=Counter()
        for ex,v in terms.items():
            for c in moments.cotrees:
                new[tuple(a+b for a,b in zip(ex,c))]+=v
        terms=new
    row=Counter()
    for ex,c in terms.items():
        row[moments.index[moments.orbit(ex)]]+=c
    return row


def equations(moments):
    rows=[]
    rhs=[]
    for edge in (0,4):
        for _,row,value in face_data(moments.index,edge):
            rows.append(row)
            rhs.append(value+(F(),))
    assert len(rows)==751
    for ex,j in moments.index.items():
        dual=tuple(4-ex[k] for k in DUAL)
        values=[a for e in (ex,dual) if (a:=star_value(e)) is not None]
        if values:
            assert len(set(values))==1
            rows.append({j:1})
            rhs.append(values[0]+(F(),))
    assert len(rows)==1299
    for q,ex in [(2,(0,)*8),(3,(1,)*4+(0,)*4)]:
        rows.append(lift(moments,ex,q))
        rhs.append(radial(q))
    return rows,rhs


def base_matrix(rows,basis):
    # Clear denominators in the homogeneous wheel IBP at the moment level.
    table=np.array(basis.tolist(),dtype=np.int64)
    bound=max(sum(abs(int(v)) for v in row.values()) for row in rows)*int(np.abs(table).max())
    assert bound<2**63
    return np.array([sum(c*table[k] for k,c in row.items()) for row in rows])


def recover(probes):
    moments=WheelMoments()
    homogeneous=kernel(moments.relations())
    rows,rhs=equations(moments)
    matrix=base_matrix(rows,homogeneous)
    chosen=pivots(matrix.T)
    assert len(chosen)==219 and chosen[-2:]==[1299,1300]
    independent=matrix[chosen]
    cols=pivots(independent)
    left=fmpz_mat(independent[:,cols].tolist())
    values=fmpq_mat([[str(y) for y in rhs[i]] for i in chosen])
    numerator,denominator=values.numer_denom()
    solution=left.solve(numerator)/denominator
    assert fmpq_mat(fmpz_mat(matrix[:,cols].tolist()))*solution==fmpq_mat([[str(x) for x in r] for r in rhs])
    with TemporaryDirectory() as tmp:
        folder=Path(tmp)
        weights=wheel_coefficients(moments,types(6)[10],probes,folder,compile_grid(folder))*homogeneous
    target=weights.tolist()
    restricted=fmpz_mat([[row[j] for j in cols] for row in target])
    coeff=left.transpose().solve(restricted.transpose()).transpose()
    assert coeff*fmpz_mat(independent.tolist())==weights
    # wheel_coefficients clears all serial beta moments by 7!.
    integral=(coeff*values)/factorial(7)
    return [tuple(F(str(x)) for x in row) for row in integral.tolist()]


if __name__=='__main__':
    data=json.loads((HERE/'wheel_inputs.json').read_text())
    ids=(10,11)
    assert list(PROBES)==[tuple(map(tuple,data['probes'][i])) for i in ids]
    values=recover(PROBES)
    assert values==[tuple(map(F,data['coefficients'][i])) for i in ids]
    print('Wheel source coefficients from angular modes, beta integrals and lower-loop cuts:',values)
