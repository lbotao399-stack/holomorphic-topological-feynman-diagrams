"""Independent eight-primary bulk fixture using the four-simplex cubature."""
from collections import Counter
from fractions import Fraction as F
from math import factorial,prod
from flint import fmpz_mat
from bulk_eight import SHAPES,potentials
from marked_theta import partitions,conductance
from integrated_cubature import green

MOMENTA=((1,0,0,0,0),(1,0,0,0,0),(0,1,0,0,0),(-1,-1,0,0,0),
         (0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1),(-1,0,-1,-1,-1))


def evaluate(shape):
    result=F()
    ds=[6-k for k in shape]
    base=conductance(ds)
    for groups in partitions(shape):
        path={i:h for h,g in enumerate(groups) for i in g}
        axes=[list(p[:2]) for p in MOMENTA[:4]]
        cols=[(a if path[i]!=2 else [0,0])+(a if path[i]!=0 else [0,0]) for i,a in enumerate(axes)]
        scale=int(fmpz_mat(cols).det())**2
        if scale:
            _,xy,bridge=green(tuple(path[i] for i in range(4,8)))
            for i,sg in ((2,1),(4,-1)):
                for j,p in enumerate(MOMENTA):
                    if p[0]:
                        term=xy[i][j]
                        if path[i]==path[j]:
                            term+=F(9,ds[path[i]])*bridge[i][j]
                        result+=base*scale*sg*p[0]*term
    return result/(2*prod(factorial(n) for n in Counter(shape).values()))


if __name__=='__main__':
    total=F()
    for shape in SHAPES:
        value=evaluate(shape)
        phi=potentials(MOMENTA,shapes=(shape,))
        assert value==phi[2][0]-phi[4][0]
        print(shape,value)
        total+=value
    assert total==F(1,19440)
    print('Transverse eight-primary potential pairing:',total)
