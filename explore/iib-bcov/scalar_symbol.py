"""Maximal root order of the primary differential; its two-loop trilinear symbol."""
from fractions import Fraction as F
from collections import Counter
from math import factorial,prod
from flint import fmpq_mat,fmpz_mat
from marked_theta import partitions,conductance,green_moments
from point_graph import solve_graph

MOMENTA=tuple(tuple(int(i==j) for j in range(5)) for i in (0,0,1,1,2,2,3,4))


def symbol(shape):
    degrees=[6-k for k in shape]
    common=conductance(degrees)
    total=F()
    checked=F()
    active=0
    for groups in partitions(shape):
        path={i:h for h,g in enumerate(groups) for i in g}
        survives=all(path[i]!=path[i+1] for i in (0,2,4,6))
        # Independently form the entire determinant at unequal integer times.
        times=(2,3,5,7,11,13,17,19)
        v=[sum(s*p[i] for s,p in zip(times,MOMENTA)) for i in range(5)]
        zeros=[0]*5
        entries=[v+zeros,zeros+v]
        for h,group in enumerate(groups):
            for i in group:
                p=list(MOMENTA[i])
                entries.append((p if h!=2 else zeros)+(p if h!=0 else zeros))
        assert int(fmpz_mat(entries).det())**2==(times[6]*times[7])**2*survives
        if not survives:
            continue
        active+=1
        h=path[7]
        total+=common*(F(5,16)+F(9*(shape[h]-2),120*degrees[h]))
        for j in range(7):
            xy,local=green_moments((0,0,0,0,0,0,2,2),7,j)
            checked+=common*xy
            if path[j]==h:
                checked+=common*F(9,degrees[h])*local
    assert total==checked
    # The distinguished trivalent endpoint fixes the orientation of each path.
    return active,total/prod(factorial(m) for m in Counter(shape).values())


def extremal(n):
    edges=[(i,(i+1)%6) for i in range(6)]
    size=6
    tip=5
    for step in range(1,n):
        new=[0,*range(size,size+4),tip]
        edges+=list(zip(new,new[1:]))
        tip=size
        size+=4
    degrees=[sum(i in e for e in edges) for i in range(size)]
    assert degrees[0]==n+1
    src=[(0,)*5 for _ in degrees]
    for j,i in enumerate((i for i in range(2,size) if degrees[i]==2),1):
        src[i]=(1,j,j*j,j**3,j**4)
    src[1]=(0,1,2,3,4)
    src[0]=tuple(-sum(col) for col in zip(*src))
    t=[1+(i*i+7*i)%23 for i in range(len(edges))]
    flow,phi=solve_graph(edges,src,t,0)
    cycles=[[int(i<6) for i in range(len(edges))]]
    for step in range(1,n):
        start=6+5*(step-1)
        old=start-5 if step>1 else 5
        oldsign=-1 if step>1 else 1
        cycles.append([int(start<=i<start+5)+oldsign*int(i==old) for i in range(len(edges))])
    for c in cycles:
        assert all(sum(x*(int(i==v)-int(j==v)) for x,(i,j) in zip(c,edges))==0
                   for v in range(size))
    minor=fmpq_mat([[str(c[j]*flow[j][i]) for j in range(1,len(edges))]
                    for i in range(5) for c in cycles])
    assert minor.det() and phi[1][0]>0
    return degrees[0]


if __name__=='__main__':
    total=F()
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        active,value=symbol(shape)
        print(shape,active,value)
        total+=value
    print('two-loop trilinear symbol',total)
    print('nonzero root degrees at loops 1 to 8:',[extremal(i) for i in range(1,9)])
