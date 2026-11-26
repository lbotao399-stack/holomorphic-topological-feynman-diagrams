"""Four-mark Li-Wang radial flux and its exact scalar-coefficient threefold ladder integral."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import sys

from flint import fmpq_mat
from theta_periods import blocks

IIB = Path(__file__).resolve().parents[1]/'iib-bcov'
sys.path.insert(0,str(IIB))
from point_graph import solve_graph
from scalar_fan import cofactor
from square_three import G, Graph, plan_generic, sp_monomial


def ladder_weights(genus):
    if not isinstance(genus,int) or genus<1:
        raise ValueError("Expected a positive integer genus")
    mass = F(6**(genus+1),factorial(2*genus)*factorial(2*genus+1)*factorial(2*genus+3))
    odd = sum(F(1,2*i+1) for i in range(genus))
    resistance = 4+F(2*genus,genus+2)+F(2*(2*genus+1),genus*(genus+2))*odd
    automorphisms = 8 if genus==1 else 4
    return mass,resistance,F(2,automorphisms)*mass*resistance


def shape(genus):
    return ([(2*i,2*i+1) for i in range(genus+1)] +
            [(2*i+s,2*i+s+2) for i in range(genus) for s in (0,1)])


def measured_resistances(edges, times, corners):
    nv = max(max(e) for e in edges)+1
    sources = [[int(i==a)-int(i==b) for a,b in combinations(corners,2)] for i in range(nv)]
    _,potential = solve_graph(edges,sources,times,0)
    return sum(potential[a][k]-potential[b][k] for k,(a,b) in enumerate(combinations(corners,2)))


def check_geometry(genus):
    edges = shape(genus)
    corners = (0,1,2*genus,2*genus+1)
    nv = 2*genus+2
    times = [1+(7*i+3*i*i)%29 for i in range(len(edges))]
    rung = times[:genus+1]
    rails = [sum(times[genus+1+2*i:genus+3+2*i]) for i in range(genus)]
    sources = [(0,0,0)]*nv
    simplex = [(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)]
    for k,p in zip(corners,simplex):
        sources[k] = p
    current,_ = solve_graph(edges,sources,times,0)
    cycles = G(nv,edges).C
    kirchhoff = fmpq_mat([[sum(a*b*t for a,b,t in zip(v,w,times))
                         for v in cycles] for w in cycles]).det()
    assert cofactor(current,cycles,times)*int(kirchhoff)**genus==prod(rung)**(genus-1)
    incidence = fmpq_mat([[int(j==i)-int(j==i+1) for i in range(genus)] for j in range(genus+1)])
    mat = fmpq_mat([[str(F(int(i==j),rung[i])) for j in range(genus+1)] for i in range(genus+1)])
    inv_rails = fmpq_mat([[str(F(int(i==j),rails[i])) for j in range(genus)] for i in range(genus)])
    green = (mat + incidence*inv_rails*incidence.transpose()).inv()
    assert F(str(green[0,genus]))==F(prod(rung),int(kirchhoff))
    difference = sum(green[i,i]+green[i+1,i+1]-2*green[i,i+1] for i in range(genus))
    mean = F(str(2*(green[0,0]+green[genus,genus])+difference/3))+F(2,3)*sum(rails)
    if genus<=3:
        # Integrate the independent rail fractions in each coordinate.
        # The resistance polynomial has degree two in each fraction. Use a
        # positive interior rule through degree two, avoiding zero-length edges.
        from itertools import product
        nodes = (F(1,6),F(1,2),F(5,6))
        weights = (F(3,8),F(1,4),F(3,8))
        value = 0
        for sample in product(range(3),repeat=genus):
            split = [x for i,s in enumerate(sample) for x in
                     (rails[i]*nodes[s],rails[i]*(1-nodes[s]))]
            value += prod(weights[s] for s in sample)*measured_resistances(edges,[*rung,*split],corners)
        assert value==F(str(mean))


def two_loop_descendant():
    result = []
    # The labeled Green matrices have already integrated the squared cofactor.
    for family,order in (((1,3,3),4),((2,2,3),4)):
        g = blocks(family)
        resistance = sum(g[i,i]+g[j,j]-2*g[i,j] for i,j in combinations(range(4),2))
        result.append(F(2,order)*F(str(resistance)))
    return result


def independent_integral(genus):
    edges = shape(genus)
    nv = 2*genus+2
    corners = (0,1,nv-2,nv-1)
    plan = plan_generic(Graph(nv,edges))
    result = 0
    for deleted in combinations(range(len(edges)),genus+1):
        parent = list(range(nv))
        def root(i):
            while parent[i]!=i:
                i = parent[i]
            return i
        cuts = set(deleted)
        for k,(a,b) in enumerate(edges):
            if k not in cuts:
                parent[root(a)] = root(b)
        roots = {root(i) for i in range(nv)}
        if len(roots)!=2:
            continue
        count = sum(root(v)==root(0) for v in corners)
        if count in (0,4):
            continue
        exponents = tuple(int(i in cuts)+(2*genus-2 if i<=genus else 0) for i in range(len(edges)))
        result += count*(4-count)*sp_monomial(exponents,2*genus+4,plan)
    return result*F(2,8 if genus==1 else 4)


def check():
    for genus in range(1,9):
        check_geometry(genus)
        mass,mean,value = ladder_weights(genus)
        normalized_lengths = sum(F(2*(2*genus+1)*(2*genus+3),
            (2*genus-2*i-1)*(2*genus-2*i+1)*(2*i+1)*(2*i+3)) for i in range(genus))
        assert mean == 4+F(2*(genus+1),2*genus+3)*normalized_lengths
        if genus<=4:
            assert independent_integral(genus)==value
            print("genus",genus,"ladder",value)
    sectors = two_loop_descendant()
    assert sectors[0]==ladder_weights(2)[2] and sum(sectors)==F(1,9720)
    print("complete genus-two descendant",*sectors,sum(sectors))


if __name__=='__main__':
    check()
