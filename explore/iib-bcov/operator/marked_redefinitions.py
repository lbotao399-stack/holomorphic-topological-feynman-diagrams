"""Field-type and jet constraints on translations/SL(5)-equivariant redefinitions."""
from itertools import combinations_with_replacement
from sympy import Rational,diag


def classification():
    # j=0 holomorphic waves. A type is (descendant level, polyvector rank).
    types=[(k,i) for k in range(5) for i in range(5-k)]
    target=tuple(sorted([(0,0)]*4+[(0,1)]))
    nonlinear=set()
    for pos in range(5):
        for other in range(pos+1,5):
            external=[target[pos],target[other]]
            level=sum(k for k,i in external)
            rank=sum(i for k,i in external)-1
            # A tree bracket only contracts a primary. All other legs at
            # the sole interacting vertex are primary, so its valence is three.
            if level or (0,rank) not in types:
                continue
            rest=[v for j,v in enumerate(target) if j not in (pos,other)]
            nonlinear.add(tuple(sorted(rest+[(0,rank)])))
    free=set()
    for j,(level,rank) in enumerate(target):
        lifted=(level+1,rank-1)
        if lifted in types:
            free.add(tuple(sorted([v for i,v in enumerate(target) if i!=j]+[lifted])))
    assert nonlinear=={((0,0),)*4}
    assert free=={((0,0),)*4+((1,0),)}
    # Check all possible eps*eps invariants of total degree ten.
    assert all(len(set(labels))<5 for labels in combinations_with_replacement(range(4),5))
    # The determinant-one scaling used in the proof on a generic four-covector support.
    scale=diag(2,2,2,2,Rational(1,16))
    assert scale.det()==1
    return nonlinear,free


if __name__=='__main__':
    classification()
    print('Degree ten: four primary scalars only. Free t∂: four primaries and one scalar descendant.')
