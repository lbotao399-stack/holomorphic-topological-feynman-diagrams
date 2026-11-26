"""Fifteen wheel periods determining all four-loop kernels given the lower faces."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
from random import Random
import numpy as np
from flint import fmpz_mat

from verify_identities import evaluate_basis
from contractions import basis, record, scaled_coefficients, substitution
from homogeneous import exact_kernel
from sampling import pivots


def marked_ranks(invariant,lift):
    offsets=np.cumsum([0]+[m.shape[1] for m in invariant])
    kernels=[m @ lift[x:y] for m,x,y in zip(invariant,offsets,offsets[1:])]
    ranks={}
    for count in (3,4,5):
        parts=[[] for _ in kernels]
        for marks in combinations(range(6),count):
            images=tuple(tuple(int(v==j)-int(v==marks[-1]) for j in marks[:-1]) for v in range(5))
            matrix=substitution(images)
            for group,kernel in zip(parts,kernels):
                group.append(matrix @ kernel)
        per_graph=[fmpz_mat(np.concatenate(p).tolist()).rank() for p in parts]
        total=fmpz_mat(np.concatenate([m for p in parts for m in p]).tolist()).rank()
        ranks[str(count)]={'per_graph':per_graph,'all_graphs':total}
    return ranks


def candidate_sources():
    rng=Random(721010)
    momenta=((1,0),(0,1),(-1,0),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1))
    rows=[]
    for count in range(3,7):
        for support in combinations(range(6),count):
            for _ in range(20):
                fixture=[(0,0)]*6
                p=rng.choices(momenta,k=count-1)
                last=tuple(-sum(v) for v in zip(*p))
                if last==(0,0):
                    continue
                for i,v in zip(support,p+[last]):
                    fixture[i]=v
                rows.append(fixture)
    return np.array(rows,dtype=np.int64)


def certificate():
    invariant,_,_,lift=exact_kernel()
    dims=[a.shape[1] for a in invariant]
    first=sum(dims[:10])
    wheel=invariant[10] @ lift[first:first+dims[10]]
    pts=candidate_sources()
    ev=evaluate_basis(pts,basis(5))
    # Use a small prime to find candidate rows; the determinant is checked over Z.
    prime=65521
    mod=(ev%prime) @ (wheel%prime) % prime
    chosen=pivots(mod.T)
    assert len(chosen)==15
    pts,ev=pts[chosen],ev[chosen]
    anchor=fmpz_mat(ev.tolist()) * fmpz_mat(wheel.tolist())
    assert anchor.det()!=0
    data,_=record()
    den,kernels=scaled_coefficients(data)
    coeff=(fmpz_mat(ev.tolist()) * fmpz_mat(kernels[10].tolist()))
    # Independent check of the homogeneous weight-five solution on all graphs.
    from flint import fmpq_mat
    combination=anchor.inv()*fmpz_mat([[coeff[i,2]] for i in range(15)])
    for i,b in enumerate(invariant):
        start=sum(dims[:i])
        block=b @ lift[start:start+dims[i]]
        assert fmpq_mat(fmpz_mat(block.tolist()))*combination==fmpq_mat(
            fmpz_mat([[int(x)] for x in kernels[i][:,2]]))
    return {'catalogue_index':10,'edges':data['edges'][10],
            'probes':pts.tolist(),'coefficients':[[str(F(int(coeff[i,j]),den)) for j in range(3)]
                                                for i in range(15)],
            'marked_ranks':marked_ranks(invariant,lift),
            'absolute_determinant_in_integer_kernel_basis':str(abs(anchor.det()))}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    data=json.dumps(certificate(),separators=(',',':'))+'\n'
    path=Path(__file__).with_name('wheel_inputs.json')
    if args.save:path.write_text(data)
    else:assert path.read_text()==data
