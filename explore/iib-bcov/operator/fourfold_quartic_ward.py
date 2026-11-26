"""Eight-primary fourfold Ward component on three zero-sum collinear waves."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path
import sys

from flint import fmpq_mat
from quartic_seven_plane_bare import compressed_network, dot, fixture, minors
from quartic_seven_plane_hexagon import det, normal_component as polygon, planar_factors, theta_planar_factor
from quartic_seven_plane_theta import normal_component as theta
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from point_graph import solve_graph
from scalar_fan import cofactor
from worldsheet_primary import fourfold_primary


NORMAL = ((1,0,0), (0,1,0), (0,0,1), (-1,0,1), (0,-1,-2))
STRENGTHS = (1,1,-2)


def fusion(normal,genus=4):
    out = 0
    for i, sign in ((0,1), (1,-1)):
        for j, p in enumerate(normal):
            if i != j:
                out += sign*p[2]*det([n for k,n in enumerate(normal) if k not in (i,j)])**(2*genus)
    return out


def lower_vertices():
    axes = [(0,*p) for p in NORMAL[:3]]
    planar = []
    for loops in (1,2):
        coeff = [F(), F()]
        for used in combinations(range(3), loops):
            rows = [(v,0,0,0) for i,v in enumerate(STRENGTHS) if i not in used] + axes
            vertex = fourfold_primary(4-loops, rows + [tuple(-sum(c) for c in zip(*rows))])
            for k in range(2):
                coeff[k] += vertex[k] * factorial(loops)*prod(STRENGTHS[i]**2 for i in used)
        planar.append(coeff)
    assert planar == [[F(22599229,252194342400), -F(143249,1945944000)], [F(1,1296), 0]]
    # Four extra identical fivefold waves multiply each sector by the Abel degree 24.
    coeff5 = ((F(-841479311,211843247616000),F(2104429,8826801984000),F(-70595227,31384184832000)),
              (F(53359,16345929600),F(-21479,108972864000),F(40357,21794572800)))
    lift = [(1,0)]*4 + [(0,x) for x in STRENGTHS]
    assert [-x/24 for x in planar_factors(lift, coeff5)] == planar[0]
    assert theta_planar_factor(lift)/24 == planar[1][0]
    pentagon = -polygon(NORMAL,0,1)
    branch, interior = (theta(NORMAL,kind) for kind in ('branch','interior'))
    assert (pentagon, branch, interior) == (F(3,280), -F(159293,91683900000), -F(1423,63504000))
    return tuple(pentagon*p+(branch+interior)*q for p,q in zip(*planar))


def collinear_square(edges,cycles,lengths,normal,strengths,eligible=None):
    """Integrate one fewer collinear mark than cycles, summing all segment assignments."""
    count = len(strengths)
    if count != len(cycles)-1 or count < 2 or sum(strengths):
        raise ValueError("Expected at least two zero-sum line strengths and one additional cycle")
    cols, potentials, _ = compressed_network(edges,cycles,lengths,normal)
    total = F()
    for subset in combinations(range(len(edges)) if eligible is None else eligible, count):
        wedge = minors([cols[i] for i in subset])
        if not any(wedge):
            continue
        mid, spread = [], F()
        for i in subset:
            x,y = [dot(wedge,potentials[v]) for v in edges[i]]
            mid.append((x+y)/2)
            spread += (count-1)*(x-y)**2/12
        total += prod(lengths[i] for i in subset)*(sum((x-y)**2 for x,y in combinations(mid,2))+spread)
    return total*factorial(count-2)*prod(strengths)**2*sum(a*a for a in strengths)/sum(lengths)**2


def restored_square(edges,cycles,lengths,normal,assigned,positions,strengths=STRENGTHS):
    full_edges, times, parents = [], [], []
    for i,(a,b) in enumerate(edges):
        cuts = sorted((s,len(normal)+k) for k,(e,s) in enumerate(zip(assigned,positions)) if e==i)
        nodes = [a,*[v for _,v in cuts],b]
        steps = [0,*[s for s,_ in cuts],1]
        for h,edge in enumerate(zip(nodes,nodes[1:])):
            full_edges.append(edge)
            times.append((steps[h+1]-steps[h])*lengths[i])
            parents.append(i)
    sources = [(0,*p) for p in normal]+[(s,0,0,0) for s in strengths]
    currents,_ = solve_graph(full_edges,[list(map(str,p)) for p in sources],times,0)
    return cofactor(currents,[[c[i] for i in parents] for c in cycles],times)**2


def check_bare():
    edges,cycles,lengths,normal = fixture()
    normal[-5:] = NORMAL
    # Exact degree-two cubature: the squared cofactor is the square of an affine
    # function of all three positions. Use its value and its three diagonal curvatures.
    for triple in combinations(range(len(edges)),3):
        reduced = collinear_square(edges,cycles,lengths,normal,STRENGTHS,triple)
        if reduced:
            break
    assert reduced
    observed = F()
    for assigned in permutations(triple):
        mid = [F(1,2)]*3
        center = restored_square(edges,cycles,lengths,normal,assigned,mid)
        averaged = center
        for j in range(3):
            a,b = mid[:],mid[:]
            a[j],b[j] = F(1,4),F(3,4)
            averaged += F(2,3)*(restored_square(edges,cycles,lengths,normal,assigned,a)
                                  + restored_square(edges,cycles,lengths,normal,assigned,b)-2*center)
        observed += prod(lengths[i] for i in assigned)*averaged
    assert observed == reduced
    full = collinear_square(edges,cycles,lengths,normal,STRENGTHS)
    assert full > reduced
    tree = fmpq_mat([[str(sum(t*u*v for t,u,v in zip(lengths,a,b))) for b in cycles] for a in cycles]).det()
    _, phi = solve_graph(edges,[list(map(str,p)) for p in normal],lengths,0)
    density = full*(phi[len(normal)-5][2]-phi[len(normal)-4][2])/F(str(tree))**4
    assert density
    return triple, density


def integrated_placement():
    """The 17-edge two-double graph, with the three collinear labels permuted."""
    from bulk_three import forest_data
    from quartic_seven_plane_double import CORE, tree_flow
    from quartic_seven_plane_double_three_two import integrate_polynomials
    from quartic_seven_plane_lift import fixture as subdivide
    from square_three import G
    from current_minor import numerator

    marks = (0,0,0,2,5,7,7,7)
    edges = subdivide(CORE,marks)
    graph = G(14,edges)
    paths, offset = [], 0
    for i in range(9):
        count = marks.count(i)+1
        paths.append(tuple(range(offset,offset+count)))
        offset += count
    groups = (sum(paths[:4],()), paths[4], sum([paths[i] for i in (7,8,5,6)],()))
    splits = ((len(paths[0]),len(paths[1])),(len(paths[7]),len(paths[8])))
    ctx,_,dist = forest_data(edges,14)
    waves = [(s,0,0,0) for s in STRENGTHS]+[(0,*q) for q in NORMAL]
    order = [4,0,6,5,2,1,3,7]
    plus, minus = (6+order.index(i) for i in (3,4))
    voltage = -sum(waves[i][3]*(dist[plus][6+j]-dist[minus][6+j]) for j,i in enumerate(order))/2
    assert all(c.denominator == 1 for _,c in voltage.terms())
    times = [2+(i*i+3*i)%17 for i in range(17)]
    tree = F(str(fmpq_mat([[sum(t*a*b for t,a,b in zip(times,c,d)) for d in graph.C] for c in graph.C]).det()))
    out = []
    for negative in (1,4,5):
        indices = order[:]
        indices[4],indices[negative] = indices[negative],indices[4]
        sources = [(0,)*4]*6+[waves[i] for i in indices]
        num = numerator(graph,sources,tree_flow(graph,sources))
        volt = num.context().from_dict({m:int(c) for m,c in voltage.terms()})
        currents,potentials = solve_graph(edges,sources,times,0)
        assert abs(int(num(*times))) == cofactor(currents,graph.C,times)*tree**4
        assert int(volt(*times)) == tree*(potentials[plus][3]-potentials[minus][3])
        # Integrate one squared polynomial against the voltage, with U^13 below.
        gram = [[int(i==j==0) for j in range(4)] for i in range(4)]
        zeros = [num.context().constant(0)]*3
        out.append(integrate_polynomials((graph,ctx,None,gram,dist,groups,splits),[num,*zeros],volt))
    assert out == [-F(1873392799,321949590346137600000000),
                   -F(1529408143,143088706820505600000000),
                   -F(376313593,257559672276910080000000)]
    return 2*sum(out)


if __name__ == '__main__':
    assert fusion(NORMAL) == 252
    lower = lower_vertices()
    assert lower == (F(1057669573,1123411161600000), -F(143249,181621440000))
    print("Tree coefficient 252; lower sectors in 1,zeta(3) order", *lower)
    print("Bare placements, literal replay triple and full fixed network density", *check_bare())
    labeled = integrated_placement()
    assert labeled == -F(80346571,2235761044070400000000)
    print("Six fully integrated labeled bare graphs", labeled)
