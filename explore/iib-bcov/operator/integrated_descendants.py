"""Scalar coefficient of a constant top-polyvector translation of the full BCOV potential."""
from fractions import Fraction as F
from itertools import product
from math import factorial,prod
from collections import Counter


def assignments(arity,level):
    if level<0:
        return
    for a4 in range(level//4+1):
        for a3 in range((level-4*a4)//3+1):
            for a2 in range((level-4*a4-3*a3)//2+1):
                a1=level-4*a4-3*a3-2*a2
                a0=arity-a1-a2-a3-a4
                if a0>=0:
                    yield (a0,a1,a2,a3,a4)


def source_factor(a):
    return prod(factorial(x)*factorial(j)**x for j,x in enumerate(a))


def translated(r,limit):
    return {a:F(factorial(m-2),factorial(r)*source_factor(a))
            for m in range(r+2,limit+1) for a in assignments(m,m-r-2)}


def derivative(poly,j):
    result={}
    for a,c in poly.items():
        if a[j]:
            b=list(a)
            b[j]-=1
            result[tuple(b)]=a[j]*c
    return result


def u_power(s,limit):
    # Lagrange's formula for u^s, u=sum phi_k u^k/k!.
    return {a:F(s*factorial(m-1),source_factor(a))
            for m in range(s,limit+1) for a in assignments(m,m-s)}


def multiply(left,right,limit):
    out=Counter()
    for a,c in left.items():
        for b,x in right.items():
            if sum(a)+sum(b)<=limit:
                out[tuple(v+w for v,w in zip(a,b))]+=c*x
    return dict(out)


def verify():
    limit=11
    for r in range(6):
        potential=translated(r,limit)
        for j in range(5):
            u=u_power(r+j+1,limit-1)
            divisor=(r+j+1)*factorial(r)*factorial(j)
            assert derivative(potential,j)=={a:c/divisor for a,c in u.items()}
        for j in range(4):
            f={a:v/factorial(j+1) for a,v in u_power(j+1,limit-2).items()}
            for k in range(5):
                lhs=derivative(derivative(potential,j+1),k)
                rhs=multiply(f,derivative(derivative(potential,0),k),limit-2)
                assert lhs==rhs
    for a,c in translated(1,limit).items():
        if a[2:]==(0,0,0):
            assert a[0]==3 and c==F(a[1]+1,6)
    for a,c in translated(0,limit).items():
        if a[2:]==(0,0,0):
            assert a[0]==2 and c==F(1,2)
    print('Six background levels through arity eleven: all five partial derivatives')
    print('and all twenty scalar/Beltrami integrability equations verified.')


if __name__=='__main__':
    verify()


def integrated_symbol():
    from sympy import Matrix
    from itertools import permutations
    rows=((1,2,0,0,1),(0,-1,2,0,0),(0,0,1,2,0),(0,0,0,1,1),(1,0,0,0,1))
    vector=tuple(-sum(col) for col in zip(*rows))
    pol=(1,2,3,-1,1)
    dot=lambda a,b:sum(v*w for v,w in zip(a,b))
    det=int(Matrix(rows).det())
    assert det
    raw=F()
    for s,row in enumerate(rows):
        active=(*rows[:s],*rows[s+1:],vector)
        d=int(Matrix(active).det())
        assert d*d==det*det
        raw+=F(d*d*(dot(vector,pol)-dot(row,pol)),24)
    assert raw==F(det*det*dot(vector,pol),4)
    # Direct coefficient of the five independent scalar waves in det(Hess phi).
    hessian=0
    for sigma in permutations(range(5)):
        matrix=[[rows[sigma[i]][i]*rows[sigma[i]][j] for j in range(5)] for i in range(5)]
        hessian+=int(Matrix(matrix).det())
    assert hessian==det**2
    for j in range(6):
        degree_G=5
        dual=-(-1)**degree_G
        preceding=(-1)**(5-j)
        divergence=(-1)**j
        assert raw+F(dual*preceding*divergence*hessian*dot(vector,pol),4)==0
    print('Five spectator choices and all six Dolbeault degree distributions cancel against the free descendant.')


if __name__=='__main__':
    integrated_symbol()
