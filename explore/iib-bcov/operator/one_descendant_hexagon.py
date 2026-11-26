"""Effective scalar and vector slots in the integrated-cubic one-loop hexagon."""
from fractions import Fraction as F
from math import factorial
from sympy import Matrix


def add(*vectors):
    return tuple(map(sum,zip(*vectors)))


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def determinant(rows):
    return int(Matrix(rows).det())


def bare(scalars,k,v):
    # The four scalar slots and one vector at the five bulk vertices.
    return F(determinant([*scalars,k])**2*dot(add(*scalars,k,k),v),24)


def original(p,k,v,q):
    parts=[F(),F(),F()]
    for i,scalar in enumerate(p):
        others=[s for j,s in enumerate(p) if j!=i]
        # F_1,00=(1+2b)phi at the root, u=(1+b)phi and w=(1+b)mu in the bulk.
        parts[0]+=2*bare(others,k,v)
        for j in range(4):
            replaced=others[:]
            replaced[j]=add(replaced[j],q)
            parts[1]+=bare(replaced,k,v)
        parts[2]+=bare(others,add(k,q),v)
    # On a transverse wave the nonlinear variation of b det Hess(phi) has Qphi=mu.dphi.
    primitive=sum(dot(s,v)*determinant([add(s,k) if j==i else x for j,x in enumerate(p)])**2
                  for i,s in enumerate(p))
    return parts,primitive


def transformed(p,k,v,q):
    root=sum(bare([s for j,s in enumerate(p) if j!=i],k,v) for i in range(5))
    # Qx=w.dx at transverse w; Qu=w.du.
    primitive=dot(q,v)*determinant(p)**2+sum(
        dot(s,v)*determinant([add(s,k) if j==i else x for j,x in enumerate(p)])**2
        for i,s in enumerate(p))
    return root,primitive


def check():
    p=[tuple(int(i==j) for j in range(5)) for i in range(5)]
    k,v=p[0],p[1]
    q=tuple(-x for x in add(*p,k))
    assert dot(k,v)==0
    parts,b=original(p,k,v,q)
    print('Original fields, placements root/scalar/vector',parts,'primitive',b)
    first,g=transformed(p,k,v,q)
    print('Variables x,u,w, hexagon and primitive',first,g)
    assert first==F(1,24) and g==0
    assert sum(parts)+F(1,4)*b==first
    for order in range(12):
        # Derivative of [(1-b)^(-6)-1]/24 compared with (1-b)^(-7)/4.
        from math import comb
        assert F((order+1)*comb(order+6,5),24)==F(comb(order+6,6),4)


if __name__=='__main__':
    check()
