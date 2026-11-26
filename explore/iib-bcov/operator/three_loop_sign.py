"""Integer replay of the global real sign of the nine-primary BCOV trace."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod
from pathlib import Path
import json

from flint import fmpq_mat, fmpz_mat
import numpy as np

from three_loop_soft_jet import DOUBLE_SOFT, SUBSETS

ROOT = Path(__file__).resolve().parent


def covers(word):
    return len(set().union(*(set(SUBSETS[i]) for i in word))) == 9


def incidence_type(word):
    # Relabeling the nine slots permutes their three incidence bits.
    return min(tuple(sorted(sum(int(i in SUBSETS[j]) << k for k,j in enumerate(order))
                            for i in range(9))) for order in permutations(word))


def verify_squares(reps):
    # The first five-row minor can always be carried to subset zero.
    shapes = {incidence_type((0,j,k)) for j,k in combinations_with_replacement(range(126),2)
              if covers((0,j,k))}
    assert len(shapes)==len(reps)==10 and shapes=={incidence_type(w) for w in reps}


def minor_images(group):
    images = np.empty((126,len(group)),dtype=np.int8)
    subset_id = np.zeros(512,dtype=np.int16)
    for j,subset in enumerate(SUBSETS):
        subset_id[sum(1<<i for i in subset)] = j+1
    for j,subset in enumerate(SUBSETS):
        labels = group[:,subset]
        masks = sum(1<<labels[:,i].astype(np.int16) for i in range(5))
        signed = subset_id[masks]
        for left,right in combinations(range(5),2):
            signed[labels[:,left]>labels[:,right]] *= -1
        images[j] = signed
    return images


def cosets(group,sets,alternating):
    canonical = group.copy()
    signs = np.ones(len(group),dtype=np.int8)
    order = prod(factorial(len(s)) for s in sets)
    for s in sets:
        canonical[:,s] = np.sort(group[:,s],axis=1)
        if alternating:
            for i,j in combinations(s,2):
                signs[group[:,i]>group[:,j]] *= -1
    # Base nine encodes all nine images without collision or overflow.
    keys = canonical @ np.array([9**i for i in range(9)],dtype=np.int64)
    _,ids,counts = np.unique(keys,return_inverse=True,return_counts=True)
    assert np.all(counts==order)
    return ids.astype(np.int32),order,signs


def young_cosets(group,shape):
    rows, offset = [],0
    for width in shape:
        rows.append(list(range(offset,offset+width)))
        offset += width
    assert offset==9 and list(shape)==sorted(shape,reverse=True)
    cols = [[row[j] for row in rows if j<len(row)] for j in range(shape[0])]
    r, nr, _ = cosets(group,rows,False)
    c, nc, signs = cosets(group,cols,True)
    col_order = np.argsort(c)
    return np.argsort(r).astype(np.int32), nr, r[col_order], nc, signs[col_order]


def projected(values,groups):
    row_order,nr,row_by_column,nc,signs = groups
    assert nr*nc*int(np.max(abs(values))) < 2**63
    row_sums = values[row_order].reshape(-1,nr).sum(axis=1)
    return (signs*row_sums[row_by_column]).reshape(-1,nc).sum(axis=1)


def evaluate(data):
    group = np.asarray(list(permutations(range(9))),dtype=np.int8)
    images = minor_images(group)
    reps = list(map(tuple,data['squares']))
    verify_squares(reps)
    blocks = data['blocks']
    young = [young_cosets(group,b['shape']) for b in blocks]
    words = set(reps).union(*(set(map(tuple,b['words'])) for b in blocks))
    assert all(len(w)==3 and covers(w) for w in words)
    design,trace,squares = [],[],[]
    for p in data['points']:
        assert len(p)==9 and all(len(q)==5 for q in p) and not any(map(sum,zip(*p)))
        minors = [int(fmpz_mat([p[j] for j in s]).det()) for s in SUBSETS]
        assert factorial(9)*max(map(abs,minors))**6 < 2**63
        table = np.array([0,*minors,*(-v for v in reversed(minors))],dtype=np.int64)

        def cubic(w):
            i,j,k = w
            return table[images[i]]*table[images[j]]*table[images[k]]

        values = {w:cubic(w) for w in words}
        squares.append(sum(int(values[w]@values[w]) for w in reps))
        trace.append([int(cubic(w[:3])@cubic(w[3:])) for w in DOUBLE_SOFT])
        row = []
        for b,groups in zip(blocks,young):
            vectors = [list(map(int,projected(values[tuple(w)],groups))) for w in b['words']]
            matrix = fmpz_mat(vectors)
            gram = matrix*matrix.transpose()
            for i,e in enumerate(b['scales']):
                for j in range(i,len(vectors)):
                    row.append(str(Fraction(int(gram[i,j])*(1+(i!=j)),1<<(e+b['scales'][j]))))
        design.append(row)
    assert int(fmpz_mat(trace).det())%1009 == 999
    return fmpq_mat(design),fmpq_mat(trace),squares


def verify():
    data = json.loads((ROOT/'three_loop_sign_certificate.json').read_text())
    periods = json.loads((ROOT/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,periods['basis_words']))==DOUBLE_SOFT
    low,high = map(Fraction,data['endpoints'])
    partial = sum(Fraction(1,j**3) for j in range(1,257))
    assert low < partial+Fraction(1,2*257**2) <= partial+Fraction(1,2*256**2) < high

    design,trace,squares = evaluate(data)
    targets = trace*fmpq_mat(periods['rational_and_zeta'])
    fixed = data['fixed']
    pivots = [i for i,x in enumerate(fixed) if x is None]
    assert design.nrows()==len(pivots)==18 and design.ncols()==len(fixed)==68
    independent = fmpq_mat([[design[i,j] for j in pivots] for i in range(18)])
    for endpoint in (low,high):
        target = fmpq_mat([[str(-Fraction(str(a))-endpoint*Fraction(str(b))
                                -Fraction(s,1<<data['square_power']))]
                           for (a,b),s in zip(targets.tolist(),squares)])
        coeffs = fmpq_mat([[x or '0'] for x in fixed])
        exact = independent.solve(target-design*coeffs)
        for j,i in enumerate(pivots):
            coeffs[i,0] = exact[j,0]
        assert design*coeffs==target
        pos = 0
        upper = Fraction(1,1<<data['square_power'])
        for block in data['blocks']:
            size = len(block['words'])
            gram = fmpq_mat(size,size)
            for i in range(size):
                for j in range(i,size):
                    gram[i,j] = gram[j,i] = coeffs[pos,0]
                    pos += 1
            for k in range(1,size+1):
                assert fmpq_mat([[gram[i,j] for j in range(k)] for i in range(k)]).det()>0
            shape = block['shape']
            rows = prod(factorial(width) for width in shape)
            columns = prod(factorial(sum(j<width for width in shape)) for j in range(shape[0]))
            upper += rows**2*columns*sum(
                abs(Fraction(str(gram[i,j])))/(1<<(e+f))
                for i,e in enumerate(block['scales']) for j,f in enumerate(block['scales']))
        assert pos==len(fixed)
        assert upper < Fraction(1,1<<19)
    print('18 unrestricted trace rows, determinant 999 modulo 1009.')
    print('Both endpoints: 16 positive definite rational Gram blocks and the positive ten-square term.')


if __name__=='__main__':
    verify()
