"""Second disk with the u-scalar and an arbitrary divergence-free descendant in KS coordinates."""
from collections import Counter
from fractions import Fraction as F
from descendant_disk import operation,hochschild


def scalar_bracket(k,ell):
    m=k+2*ell
    out=Counter()
    # Constant scalar contributions are a Hochschild boundary of the first disk,
    # hence vanish on the divergence-free input. Retain derivatives hitting phi.
    for word,value in operation(k,ell).items():
        if -1 in word:
            continue
        for j in word:
            remaining=tuple(i if i<j else i+1 for i in word if i!=j)
            sign=(-1)**(m+j+word.index(j))
            out[(-1,*remaining)]+=sign*value/2
    return {key:value for key,value in out.items() if value}


def verify(k,ell):
    m=k+2*ell
    actual=scalar_bracket(k,ell)
    # -1 now labels the scalar's momentum q. The other incoming tensor is divergence-free.
    projected={} if not k else {(-1,*key):value for key,value in operation(k-1,ell+1).items() if -1 not in key}
    assert actual=={key:(-1)**k*F(m+1,2)*value for key,value in projected.items()}
    closed={key:(-1)**k*value for key,value in projected.items()}
    defect={key:actual[key]-value for key,value in closed.items() if actual[key]!=value}
    factor=(-1)**k*F(m-1,2)
    replay={key:factor*value for key,value in hochschild(k,ell).items() if factor*value}
    assert defect==replay
    return factor


if __name__=='__main__':
    for ell in range(4):
        for k in range(6):
            verify(k,ell)
    print('24 scalar-descendant second disk Ward identities reproduced')
