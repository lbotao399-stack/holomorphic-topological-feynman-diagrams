"""Rooted primary BCOV hexagon on two scalar letters, four scalars and one vector."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
from itertools import permutations
from math import factorial,prod
from collections import Counter
from functools import cache

import sympy as sp


def sign(order):
    return (-1)**sum(x>y for i,x in enumerate(order) for y in order[i+1:])


def deltas():
    # theta[root]=0. Edge e runs from e to e+1 cyclically in the vertex labels 0..5.
    return [[int(e==v)-int((e+1)%6==v) for v in range(1,6)] for e in range(6)]


@cache
def determinant(omitted,at):
    rows=[row for i,row in enumerate(deltas()) if i not in omitted]
    if at:
        rows.append([int(v==at) for v in range(1,6)])
    return int(sp.Matrix(rows).det())


@cache
def tensor(at,axis):
    # Bulk traces are extracted in vertex-major order; five coordinate blocks
    # transpose five rank-five blocks with even parity.
    out={}
    for missing in set(permutations((*range(5),axis))):
        select=[]
        for j in range(5):
            skip=tuple(e for e,c in enumerate(missing) if c==j)
            select.append(determinant(skip,at if j==axis else 0))
        if not all(select):
            continue
        colors=[i for e in missing for i in range(5) if i!=e]+[axis]
        out[missing]=(-1)**sum(missing)*sign(colors)*prod(select)
    return out


def currents(momenta,t):
    total=sum(t)
    partial=[[sum(m[j] for m in momenta[:i]) for j in range(5)] for i in range(6)]
    average=[sum(F(w,total)*p[j] for w,p in zip(t,partial)) for j in range(5)]
    return [[x-c for x,c in zip(row,average)] for row in partial]


def coefficient_sign(e):
    others=[i for i in range(6) if i!=e]
    colors=[j for i in range(6) for j in range(5) if i==e or j!=others.index(i)]
    return (-1)**sum(range(5))*sign(colors)*prod(determinant((i,),0) for i in others)


def check_tensors():
    assert [coefficient_sign(e) for e in range(6)]==[-1,1,-1,1,-1,1]
    source=Path(__file__).resolve().parents[2]/'type-i'/'hexagon.py'
    spec=importlib.util.spec_from_file_location('type_i_heat_hexagon',source)
    heat=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(heat)
    # Polygon vertex labels start at p1, end at the operator. Rotate the edge list.
    g=heat.polygon(6)
    base=[[int(i==j) for i in range(5)] for j in range(5)]
    extra=[[2,-1,1,0,1],[1,2,0,0,-1],[0,-1,3,0,1],[-1,0,1,2,0],[0,2,1,1,1]]
    for momenta,t in [(base,(1,2,3,4,5,6)),(extra,(7,2,1,4,6,3))]:
        y=currents(momenta,t)
        determinant_p=int(sp.Matrix(momenta).det())
        heat_values=heat.edge_forms(g,y[1:]+y[:1])
        assert heat_values==[-F(w,sum(t))*determinant_p for w in t[1:]+t[:1]]
        for v in range(1,6):
            for j in range(5):
                x=sum(c*prod(y[e][i] for e,i in enumerate(miss))
                      for miss,c in tensor(v,j).items())
                s=sum(F(t[e],sum(t))*y[e][j] for e in range(v,6))
                assert x==-determinant_p*s


def compositions(total,size):
    if size==1:
        yield (total,)
    else:
        for j in range(total+1):
            for tail in compositions(total-j,size-1):
                yield (j,)+tail


def quadrature():
    # Check the rule on every homogeneous degree-seven monomial; sum(t)=1.
    nodes=[]
    for j in range(4):
        den=12-2*j
        w=F((-1)**j*den**7,64*factorial(j)*factorial(12-j))
        for beta in compositions(3-j,6):
            nodes.append((tuple(2*k+1 for k in beta),den,w))
    for ex in compositions(7,6):
        exact=F(prod(factorial(k) for k in ex),factorial(12))
        assert sum(w*F(prod(v**a for v,a in zip(t,ex)),den**7) for t,den,w in nodes)==exact
    return nodes


def ordered_moments():
    nodes=quadrature()
    momenta=[[int(i==j) for i in range(5)] for j in range(5)]
    table=[[F(0) for _ in range(5)] for _ in range(5)]
    for t,den,w in nodes:
        y=currents(momenta,t)
        values=[[int(x*den) for x in row] for row in y]
        assert all(F(v,den)==x for r,s in zip(values,y) for v,x in zip(r,s))
        for v in range(1,6):
            for j in range(5):
                # The geometric density at these momenta is -1.
                scalar=-sum(c*prod(values[e][k] for e,k in enumerate(miss))
                            for miss,c in tensor(v,j).items())
                table[v-1][j]+=w*F(scalar,den**6)
    expected=[[F(min(i,j)*(6-max(i,j)),42*factorial(5)) for i in range(1,6)] for j in range(1,6)]
    assert table==expected
    return table


def polarized(table):
    result=Counter()
    for assignments in permutations(range(5)):
        j=assignments.index(4)
        for pos,original in enumerate(assignments):
            result[original]+=table[j][pos]/2
    assert [result[i] for i in range(5)]==[F(1,24)]*4+[F(1,12)]
    return result


if __name__=='__main__':
    check_tensors()
    polarized(ordered_moments())
    print('Six heat facets, 50 ordered exterior tensors, 792 degree-seven simplex monomials,')
    print('25 directly integrated moments and 120 polarized orders. Scalar coefficient 1/24; vector 1/12.')
