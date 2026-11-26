"""Compare the Dirichlet bridge to Fernandez C.10 and verify higher longitudinal moments."""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, prod

from chain import gauge_box


def choose(n,k):
    return comb(n,k) if 0<=k<=n else 0


def published(hol,bar):
    r,s,m,n=hol
    k,h,l,j=bar
    if r+m!=k+l+1 or s+n!=h+j+1:
        return F(0),F(0)
    first=second=F(0)
    for a,c in product(range(m+1),range(n+1)):
        mult=(choose(l,a)*choose(j,c-1)-choose(l,a-1)*choose(j,c))*comb(m,a)*comb(n,c)
        if not mult:
            continue
        nums=(m+n-a-c,a+c-1,a,c,r+m-a,s+n-c,1+k+h,1+l+j-a-c)
        dens=(m+n,1+m+r+n+s-a-c,k,h,l,j)
        term=F(mult*prod(map(factorial,nums)),prod(map(factorial,dens)))
        first+=term
        second+=F(k+h-l-j+a+c,2+m+r+n+s-a-c)*term
    return -(first+second)/32,-(first-second)/32


def verify():
    count=0
    for total in range(2,9):
        for r in range(total+1):
            for s in range(total-r+1):
                for m in range(total-r-s+1):
                    n=total-r-s-m
                    for k in range(r+m):
                        for h in range(s+n):
                            hol=(r,s,m,n)
                            bar=(k,h,r+m-1-k,s+n-1-h)
                            assert gauge_box(hol,bar)==published(hol,bar),(hol,bar)
                            count+=1
    # With unit transverse determinant the simplex is uniform. The means
    # in x^0,y^0 are the two ordered positions along the unit segment.
    for p,q in product(range(5),repeat=2):
        expected=[F(0)]*(p+q+1)
        for i,j in product(range(p+1),range(q+1)):
            # Integrate x^i y^j over 0<x<y<1 and expand (z-w)^(i+j).
            weight=F(comb(p,i)*comb(q,j),16*(i+1)*(i+j+2))
            for h in range(i+j+1):
                expected[h]+=weight*comb(i+j,h)*(-1)**(i+j-h)
        assert gauge_box((0,1,1,0),(0,0,0,0),p,q)==tuple(expected)
    print('published transverse coefficients',count,'ordered longitudinal tests',25)


if __name__=='__main__':
    verify()
