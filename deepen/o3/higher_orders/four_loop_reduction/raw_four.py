"""Position-space heat-form expansion at four loops, compared against the cycle Pfaffian."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
import json
import sys
from flint import fmpz_mat, fmpq_mat, fmpq

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[3]/'explore/ht-families/verify'))
import vcore


def reorder(groups):
    # Increasing keys mean: time generators first, then the remaining coordinate
    # generators ordered by coordinate type and edge. The input is edge-major.
    seen = inversions = 0
    for e in range(13):
        for j,mask in enumerate(groups):
            key = e if (mask>>e)&1 else 13*(j+1)+e
            inversions += (seen >> (key+1)).bit_count()
            seen |= 1 << key
    return (-1)**inversions


def direct(net,t,source,pivot):
    incidence = fmpq_mat(net.Bmat())
    scaled = fmpq_mat([[x/fmpq(t[i]) for x in row] for i,row in enumerate(net.Bmat())])
    laplace = incidence.transpose()*scaled
    project = scaled*laplace.inv()
    covariance = project*scaled.transpose()/2
    currents = project*fmpq_mat(source)
    U = laplace.det()*prod(t)
    full = (1<<13)-1
    wedge = {}
    for block in combinations([i for i in range(13) if i!=pivot],4):
        mask=sum(1<<i for i in block)
        det=int(fmpz_mat([row for i,row in enumerate(net.Bmat()) if not mask>>i&1]).det())
        if not det:
            continue
        i,j,k,l=block
        top=(covariance[i,j]*covariance[k,l]+covariance[i,k]*covariance[j,l]
             +covariance[i,l]*covariance[j,k])
        wedge[mask]=(det*top,det*prod(currents[e,0] for e in block),
                     det*prod(currents[e,1] for e in block))
    value=fmpq(0)
    for mask,(top,_,_) in wedge.items():
        rest=full^mask^(1<<pivot)
        for first,(_,y,_) in wedge.items():
            if mask&first:
                continue
            second=rest^first
            if second in wedge:
                value+=reorder((mask,first,second))*top*y*wedge[second][2]
    return F(str(value/((-1)**pivot*t[pivot]*U**2))),F(str(U))


def check(indices):
    records=json.loads((HERE.parent/'four-loop-graphs.json').read_text())
    t=[3,2,7,5,11,2,3,13,7,5,3,2,11]
    source=[[0,0]]*6+[[2,1],[1,4],[-2,3]]
    found=0
    for i in indices:
        edges=[edge for chain in records[i]['chains'] for edge in zip(chain,chain[1:])]
        net=vcore.Net(10,edges)
        pivot=i%13
        raw,U=direct(net,t,source,pivot)
        expected,U2,_=vcore.density_slice(net,net.cycle_basis(),t,source,e0=pivot)
        assert (raw,U)==(expected,U2),(i,raw,expected,U,U2)
        found+=bool(raw)
    print('raw heat-kernel facets',len(indices),'nonzero',found)
    assert found


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--all',action='store_true')
    args=p.parse_args()
    check(range(50) if args.all else (0,3,11,39))
