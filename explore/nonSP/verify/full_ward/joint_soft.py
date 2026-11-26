"""Order-two source zero in graph 9 and its remaining homogeneous Ward directions."""
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from flint import fmpz_mat

from contractions import basis,record,scaled_coefficients,substitution
from homogeneous import exact_kernel
from product_forests import placements


def remaining(kernels):
    words=basis(5)
    forbidden=[j for j,word in enumerate(words) if sum(i in (3,4) for col in word for i in col)<2]
    assert len(forbidden)==63
    _,physical=scaled_coefficients(record()[0])
    assert not np.any(physical[9][forbidden])
    restricted=fmpz_mat(kernels[9,forbidden].tolist())
    nspace,dimension=restricted.nullspace()
    assert restricted.rank()==5 and dimension==10
    kept=fmpz_mat([[nspace[i,j] for j in range(dimension)] for i in range(15)])
    assert restricted*kept==fmpz_mat(len(forbidden),dimension)
    return kept


def after_cut(kernels):
    kept=remaining(kernels)
    blocks=[fmpz_mat(k.tolist())*kept for k in kernels]
    def ranks(selected):
        return fmpz_mat([list(row) for i in selected for row in blocks[i].tolist()]).rank()
    print('individual graph restriction ranks:',[ranks((i,)) for i in range(13)],flush=True)
    for group in ((4,6,7,9),(8,10),(11,12),(8,10,11,12)):
        print('group',group,ranks(group),flush=True)
    _,graphs=record()
    forests=0
    for i,g in enumerate(graphs):
        for a,b,roots in placements(g):
            labels=[a]+roots
            images=tuple(tuple(int(v==j)-int(v==b) for j in labels) for v in range(5))
            mat=substitution(images)
            if mat.shape[0]:
                assert (fmpz_mat(mat.toarray().tolist())*blocks[i]).is_zero()
            forests+=1
    assert forests==23
    print('product-forest restriction rank after joint zero: 0 over 23 forests',flush=True)
    coeff=blocks[10]
    for n in range(3,7):
        rows=[]
        for supported in combinations(range(6),n):
            source=np.zeros((6,n-1),dtype=np.int64)
            for j,i in enumerate(supported[:-1]):
                source[i,j]=1;source[supported[-1],j]=-1
            mat=substitution(tuple(tuple(int(k) for k in row) for row in source))
            if mat.shape[0]:
                rows.extend((fmpz_mat(mat.toarray().tolist())*coeff).tolist())
        print('wheel 10 at most',n,'marked vertices:',fmpz_mat(rows).rank(),flush=True)


if __name__=='__main__':
    if len(sys.argv)>1:
        all_kernels=np.load(sys.argv[1])
    else:
        invariant,_,_,lift=exact_kernel()
        offsets=np.cumsum([0]+[g.shape[1] for g in invariant])
        all_kernels=np.stack([g@lift[a:b] for g,a,b in zip(invariant,offsets,offsets[1:])])
    after_cut(all_kernels)
