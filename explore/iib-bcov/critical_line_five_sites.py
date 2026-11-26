"""Five hard sites for the mixed four-loop Ward potential on the critical line."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial,prod

from flint import fmpq_mat
from planar_cubic_primary import ROOT
from fourfold_mixed_soft_ward import hard_bipartite,integrated_soft_placements
from point_graph import solve_graph

AXES=tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
DIRECTION=(0,1,0,1)


def complementary_current(rays,cycles,lengths):
    """Remove the soft direction; the projected hard columns have one relation."""
    quotient=[[1]*len(lengths)]
    for coordinate in ([q[0] for q in rays],[q[2] for q in rays],[q[1]-q[3] for q in rays]):
        quotient.extend([[str(c*a) for c,a in zip(cycle,coordinate)] for cycle in cycles])
    assert len(quotient)+1==len(lengths)
    cofactors=[(-1)**j*F(str(fmpq_mat([row[:j]+row[j+1:] for row in quotient]).det()))
               for j in range(len(lengths))]
    current=[sum(c*q[3]*a for c,q,a in zip(cycle,rays,cofactors)) for cycle in cycles]
    M=fmpq_mat([[str(sum(l*a*b for l,a,b in zip(lengths,row,col))) for row in cycles]
                for col in cycles])
    z=fmpq_mat([[str(a)] for a in current])
    return F(str((z.transpose()*M.solve(z))[0,0]*M.det()))/sum(lengths)**2


def check(places):
    edges,cycles,lengths=hard_bipartite(places)
    nodes=(6,7,8,9,10)
    def solve(sites,charges):
        data=[(0,)*4 for _ in range(11)]
        for i,p in zip(sites,charges):
            data[i]=p
        data[sites[-1]]=tuple(-sum(col) for col in zip(*charges))
        return solve_graph(edges,data,lengths,sites[-1])

    rays,phi=solve(nodes,AXES)
    gram=integrated_soft_placements(rays,cycles,lengths,DIRECTION)/factorial(3)
    assert complementary_current(rays,cycles,lengths)==gram
    if not gram:
        return None
    # Direct sum over triples of segment placements, with every jump column retained.
    columns=list(zip(*cycles))
    def squared(triple):
        m=[[1]*len(edges)+[0]*3]
        for dim in range(4):
            for j,cycle in enumerate(cycles):
                m.append([str(c*a[dim]) for c,a in zip(cycle,rays)]+[
                    columns[k][j]*DIRECTION[dim] for k in triple])
        return F(str(fmpq_mat(m).det()))**2/sum(lengths)**2
    assert gram==sum(prod(lengths[k] for k in t)*squared(t) for t in combinations(range(len(edges)),3))

    def green(i,j):
        return phi[i][nodes.index(j)] if j!=nodes[-1] else F(0)
    def resistance(i,j):
        return green(i,i)+green(j,j)-2*green(i,j)
    excited=(7,9)
    other=(6,8,10)
    cross=sum(resistance(i,j) for i,j in product(excited,other))/6
    within=sum(resistance(i,j) for i,j in combinations(other,2))/3
    factor=F(3,4)*(2*cross-within-resistance(*excited))

    measured=F(0)
    for V,M in permutations(excited):
        for T,f,root in permutations(other):
            # s=1 after common reflection; the central second difference in t
            # extracts [s t²] from the homogeneous hard cubic in s,t.
            for t,weight in ((-1,F(1,2)),(0,-1),(1,F(1,2))):
                charges=[AXES[1],(1,0,1,1),AXES[0],(-1,0,-1,t)]
                physical,reading=solve((T,V,f,M,root),charges)
                square=integrated_soft_placements(physical,cycles,lengths,AXES[3])/factorial(3)
                assert square==(1+t)**2*gram
                # root momentum evaluated on d_1 is -1, and on -d_1+d_y is -t.
                voltage=-t*reading[T][0]-reading[V][0]+reading[V][3]
                measured+=weight*square*voltage/24  # twelve assignments, generator divisor two
    assert measured==factor*gram
    return gram,factor,measured


if __name__=='__main__':
    for places in ((0,1,3,4,8),(0,1,2,4,7),(0,0,4,5,8)):
        result=check(places)
        if result:
            print('Five sites',places,'squared cofactor, mean voltage, and product',*result,flush=True)
    assert check((0,1,3,4,8)) is not None
