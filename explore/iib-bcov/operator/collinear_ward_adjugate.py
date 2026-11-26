"""All zero-sum collinear insertions in the fourfold Ward density, at fixed edge lengths."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, factorial, prod

from flint import fmpq_mat
from fourfold_quartic_ward import NORMAL, collinear_square, det, fusion, restored_square
from quartic_seven_plane_bare import compressed_network, dot, minors
from square_three import G


def adjugate_square(edges,cycles,lengths,normal,strengths):
    count = len(strengths)
    if count != len(cycles)-1 or count < 2 or sum(strengths):
        raise ValueError("Expected at least two zero-sum line strengths and one additional cycle")
    columns,potentials,_ = compressed_network(edges,cycles,lengths,normal)
    mid = [[(x+y)/2 for x,y in zip(potentials[i],potentials[j])] for i,j in edges]
    diff = [[x-y for x,y in zip(potentials[i],potentials[j])] for i,j in edges]
    resistance = fmpq_mat([[str(sum(t*x*y for t,x,y in zip(lengths,c,d)))
                            for d in cycles] for c in cycles])
    rows = fmpq_mat([list(map(str,row)) for row in [*columns,*mid,*diff]])
    gram = rows*resistance.inv()*rows.transpose()
    size = len(edges)
    total = 0
    for i,t in enumerate(lengths):
        m,d = size+i,2*size+i
        average = gram[i,i]*gram[m,m]-gram[i,m]**2
        spread = gram[i,i]*gram[d,d]-gram[i,d]**2
        total += (count-1)*t*(average+spread/12)
    for i,j in combinations(range(size),2):
        m,n = size+i,size+j
        mix = ((gram[i,i]*gram[j,j]-gram[i,j]**2)*gram[m,n]
               -gram[i,i]*gram[m,j]*gram[n,j]-gram[j,j]*gram[m,i]*gram[n,i]
               +gram[i,j]*(gram[m,i]*gram[n,j]+gram[m,j]*gram[n,i]))
        total -= 2*lengths[i]*lengths[j]*mix
    scale = factorial(count-2)*prod(strengths)**2*sum(a*a for a in strengths)
    return F(str(resistance.det()*total))*F(scale,sum(lengths)**2)


def fixture(genus,marks):
    size = genus-1
    if genus == 3:
        core = list(combinations(range(4),2))
    else:
        core = [(v+i*size,(v+1)%size+i*size) for i in range(2) for v in range(size)]
        core += [(i,i+size) for i in range(size)]
    edges = []
    for h,(left,right) in enumerate(core):
        chain = [left,*[2*size+i for i,e in enumerate(marks) if e==h],right]
        edges.extend(zip(chain,chain[1:]))
    normal = [(0,0,0)]*(2*size)+list(NORMAL)
    cycles = G(len(normal),edges).C
    lengths = [1+(7*i*i+2*i)%17 for i in range(len(edges))]
    return edges,cycles,lengths,normal


def literal_subset(data,strengths):
    edges,cycles,lengths,normal = data
    columns,potentials,_ = compressed_network(*data)
    count = len(strengths)
    for subset in combinations(range(len(edges)),count):
        wedge = minors([columns[i] for i in subset])
        if not any(wedge):
            continue
        values = {dot(wedge,potentials[v]) for i in subset for v in edges[i]}
        if len(values)>1:
            break
    else:
        raise AssertionError('Expected an active subset of distinct segments')
    mid = [F(1,2)]*count
    total = F()
    for assigned in permutations(subset):
        center = restored_square(*data,assigned,mid,strengths)
        term = center
        for i in range(count):
            low,high = mid[:],mid[:]
            low[i],high[i] = F(1,4),F(3,4)
            term += F(2,3)*(restored_square(*data,assigned,low,strengths)
                             +restored_square(*data,assigned,high,strengths)-2*center)
        total += term
    assert total*prod(lengths[i] for i in subset) == collinear_square(*data,strengths,subset) > 0


if __name__ == '__main__':
    fusion_weights = {}
    for i,sign in ((0,1),(1,-1)):
        for j,wave in enumerate(NORMAL):
            if i != j:
                volume = abs(det([p for k,p in enumerate(NORMAL) if k not in (i,j)]))
                fusion_weights[volume] = fusion_weights.get(volume,0)+sign*wave[2]
    assert fusion_weights == {0:3,1:-4,2:1}
    for genus,marks,strengths in (
        (3,(0,1,2,3,4),(2,-2)),
        (4,(0,1,3,5,7),(1,1,-2)),
        (5,(0,2,5,7,9),(1,2,4,-7)),
    ):
        data = fixture(genus,marks)
        assert fusion(NORMAL,genus) == 4**genus-4
        value = adjugate_square(*data,strengths)
        assert value == collinear_square(*data,strengths) > 0
        literal_subset(data,strengths)
        print('Genus',genus,'complete fixed-length square',value,flush=True)
    for genus,marks,residue in (
        (6,(0,1,3,13,14),4708), (7,(0,4,7,8,10),1613),
        (8,(3,7,12,14,16),1202), (9,(3,12,13,19,20),6270), (10,(0,3,5,11,22),3167),
    ):
        strengths = (*range(1,genus-1),-(genus-2)*(genus-1)//2)
        value = adjugate_square(*fixture(genus,marks),strengths)
        assert value>0 and value.numerator*pow(value.denominator,-1,10007)%10007 == residue
        print('Genus',genus,'edge subsets',comb(3*genus+2,genus-1),'residue',residue)
