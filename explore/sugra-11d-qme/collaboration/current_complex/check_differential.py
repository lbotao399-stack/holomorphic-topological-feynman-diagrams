"""Compare integration by parts in the original quotient jet algebra with the current differential."""
from pathlib import Path
from random import Random
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'cohomology'))
from jets import Alg, unit

def beta_shift(a):
    b=a.letter('b',0,a.zero)
    for i in range(a.n):
        b=a.padd(b,a.prod(a.letter('X',i,a.zero),a.letter('g',i,a.zero)))
    return b

def normal_form(a,p):
    """Convert from beta to beta' and put the unique w=1 letter first, undifferentiated."""
    result=[{} for _ in range(a.n+1)]
    for monomial,coefficient in a.reduce(p).items():
        slot,= [k for k,t in enumerate(monomial) if a.info[t][0] in ('g','b')]
        kind,i,alpha=a.info[monomial[slot]]
        rest=monomial[:slot]+monomial[slot+1:]
        other={rest:coefficient*(-1)**(slot*(kind=='b')+sum(alpha))}
        for j,count in enumerate(alpha):
            for _ in range(count):
                other=a.D(other,j)
        if kind=='g':
            result[i]=a.padd(result[i],other)
        else:
            result[-1]=a.padd(result[-1],other)
            for j in range(a.n):
                result[j]=a.padd(result[j],a.prod(a.letter('X',j,a.zero),other),-1)
    return [a.reduce(p) for p in result]

def convection(a,p):
    q=a.Q(p)
    for i in range(a.n):
        q=a.padd(q,a.prod(a.letter('X',i,a.zero),a.D(p,i)))
    return a.reduce(q)

def verify(n):
    a=Alg(n,quot=True)
    bp=beta_shift(a)
    assert not convection(a,bp)
    for i in range(a.n):
        gamma=a.letter('g',i,a.zero)
        rhs=a.D(bp,i)
        for j in range(a.n):
            rhs=a.padd(rhs,a.prod(a.letter('X',j,unit(n,i)),a.letter('g',j,a.zero)),-1)
        assert convection(a,gamma)==a.reduce(rhs)
    rng=Random(1024+n)
    for degree in range(4):
        for _ in range(6):
            coefficients=[]
            for s in range(n+1):
                polynomial={}
                for _ in range(2):
                    factors=[]
                    for _ in range(degree+(s<n)):
                        alpha=[0]*n
                        for _ in range(rng.randrange(3)):
                            alpha[rng.randrange(n)]+=1
                        factors.append(a.letter('X',rng.randrange(n),alpha))
                    polynomial=a.padd(polynomial,a.prod(*factors))
                coefficients.append(a.reduce(polynomial))
            current,scalar=coefficients[:-1],coefficients[-1]
            original=a.prod(bp,scalar)
            for i,p in enumerate(current):
                original=a.padd(original,a.prod(a.letter('g',i,a.zero),p))
            assert normal_form(a,original)==coefficients
            expected=[]
            for i,p in enumerate(current):
                dx=convection(a,p)
                for j,pj in enumerate(current):
                    dx=a.padd(dx,a.prod(a.letter('X',i,unit(n,j)),pj),-1)
                expected.append(a.reduce(dx))
            dk=convection(a,scalar)
            for i,p in enumerate(current):
                dk=a.padd(dk,a.D(p,i))
            expected.append(a.reduce(a.padd({},dk,-1)))
            assert normal_form(a,a.Q(original))==expected
    return 24

if __name__=='__main__':
    for n in (2,3,5):
        print(f'n={n}: {verify(n)} normal-form differentials')
