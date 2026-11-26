"""Local bubble and tetrahedron moments for six-mark single-double Ward chambers."""
from fractions import Fraction
from collections import defaultdict
from itertools import combinations,product
from math import gcd

from flint import fmpz_mpoly_ctx
from fourfold_mixed_soft_ward import V,integrated_soft_placements
from fourfold_quartic_single import project_squares
from mixed_bare_ranks import hard_waves
from point_graph import solve_graph
from quartic_seven_plane_double import G,fixture,forest_data,numerator,tree_flow
from quartic_seven_plane_single import CORE,CYCLES


def chamber_data(marks):
    edges = fixture(CORE,marks)
    graph = G(12,edges)
    paths = []
    start = 0
    for e in range(9):
        stop = start+1+marks.count(e)
        paths.append(tuple(range(start,stop)))
        start = stop
    graph.C = [[c[e] for e in range(9) for _ in paths[e]] for c in CYCLES]
    ctx = fmpz_mpoly_ctx.get(tuple('t'+str(i) for i in range(15)))
    _,u,dist = forest_data(edges,12)
    times = ctx.gens()
    lengths = [sum(times[i] for i in p) for p in paths]
    return graph,ctx,u,lengths,dist,paths,(len(paths[0]),len(paths[1]))


def chamber_polynomials(data,color):
    graph,ctx,_,lengths,dist,_,_ = data
    (u,t,z),hard = hard_waves(color)
    sources = [(0,)*4]*6+[(int(a),int(b+c),int(d),int(c)) for a,b,c,d in hard]
    flow = tree_flow(graph,sources)
    minors = {(i,j):numerator(graph,sources,flow,[(3,i),(3,j)]) for i,j in combinations(range(4),2)}
    squares = []
    for i,j in combinations(range(9),2):
        det = sum(((CYCLES[a][i]*CYCLES[b][j]-CYCLES[a][j]*CYCLES[b][i])*m
                  for (a,b),m in minors.items()),ctx.constant(0))
        if not det: continue
        coefficient = int(det.content())*(1 if det.leading_coefficient()>0 else -1)
        root = det//coefficient
        factor = lengths[i]*lengths[j]*coefficient**2
        for row in squares:
            if row[0]==root:
                row[1] += factor
                break
        else:
            squares.append([root,factor])
    def resistance(left,right):
        return sum(dist[6+i][6+j] for i,j in product(left,right))
    voltage = 4*dist[6+u[0]][6+u[1]]+resistance(t,z)-resistance(u,z)-resistance(t,u)
    assert all(c.denominator==1 for _,c in voltage.terms())
    return squares,ctx.from_dict({k:int(c) for k,c in voltage.terms()})


def check_polynomials(data,color,squares,voltage):
    graph,_,u,_,_,_,_ = data
    groups,waves = hard_waves(color)
    lengths = list(range(1,16))
    tree = int(u(*lengths))
    current,_ = solve_graph(graph.edges,[(0,)*4]*6+waves.tolist(),lengths,0)
    expected = integrated_soft_placements(current,graph.C,lengths,V)
    actual = 2*sum(int(f(*lengths))*int(root(*lengths))**2 for root,f in squares)
    assert actual==tree**8*expected and actual
    a,b = groups[0]
    t,z = groups[1:]
    from itertools import permutations
    from mixed_bare_ranks import WAVES
    direct_voltage = Fraction()
    for pair1,pair2,pair3 in product(permutations((a,b)),permutations(t),permutations(z)):
        ordinals = (pair1[0],*pair2,*pair3,pair1[1])
        src = [(0,)*4]*12
        for i,h in enumerate(ordinals):
            src[6+h] = WAVES[i]
        _,phi = solve_graph(graph.edges,src,lengths,0)
        direct_voltage += phi[6+pair1[0]][3]-phi[6+pair2[0]][3]
    assert 2*int(voltage(*lengths))==tree*direct_voltage
    # U is affine in each segment. An exact nonzero numerator on U=0 rules out
    # canceling a uniform Kirchhoff factor in the full signed chamber density.
    for k in range(15):
        point = lengths[:]
        point[k] = 0
        intercept = int(u(*point))
        point[k] = 1
        slope = int(u(*point))-intercept
        if not slope: continue
        common = gcd(slope,intercept)
        point = [i*(slope//common) for i in lengths]
        point[k] = -intercept//common
        assert not u(*point)
        sample = sum(int(f(*point))*int(root(*point))**2 for root,f in squares)*int(voltage(*point))
        if sample:
            return k,point,sample%1000003
    return None


def check_first(data,squares,voltage):
    from kite_native import coefficients
    from quartic_seven_plane_single import cuts,reduce_polynomial,tetrahedron_bank
    from square_three import checked_cut
    combined = sum(root**2*factor for root,factor in squares)*voltage
    assert len(combined)==151767 and combined.total_degree()==37
    weights = reduce_polynomial(data,combined)
    assert len(weights)==185
    for p in (min(weights),max(weights)):
        assert checked_cut(p)==tetrahedron_bank()[p]
        coefficients.cache_clear()
    return cuts(weights)


def shape_chambers():
    import numpy as np
    from mixed_bare_ranks import census,integral_orbits
    _,_,physical = census(CORE,np.indices((9,)*6,dtype=np.int8).reshape(6,-1))
    groups = defaultdict(list)
    for marks,color,_,weight in integral_orbits(CORE,physical):
        groups[marks].append((color,weight))
    assert len(groups)==102
    return groups


def integrate_shape(marks,rows):
    data = chamber_data(marks)
    result = []
    for color,weight in rows:
        squares,voltage = chamber_polynomials(data,color)
        first = marks==(0,0,0,8,8,8) and color==(0,0,2,1,1,2)
        if first:
            sample = check_polynomials(data,color,squares,voltage)
            assert sample[0]==0 and sample[2]==897822
            reference = check_first(data,squares,voltage)
        if not squares:
            value = (Fraction(),Fraction())
        else:
            value = project_squares(data,squares,voltage)
        if first:
            assert value==reference
        # mean4 is four times the voltage and each pair has two soft labelings.
        result.append((color,weight,*(weight*x/2 for x in value)))
    return result


if __name__=='__main__':
    from argparse import ArgumentParser
    from quartic_scalar_response import interval
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--paths',nargs=6,type=int,default=(0,0,0,8,8,8))
    args = parser.parse_args()
    marks = tuple(args.paths)
    groups = shape_chambers()
    if marks not in groups:
        parser.error('Expected six canonical hard paths of a retained single-double placement.')
    rows = integrate_shape(marks,groups[marks])
    for row in rows: print(*row,flush=True)
    subtotal = tuple(sum(row[j] for row in rows) for j in (2,3))
    if marks==(0,0,0,8,8,8):
        assert len(rows)==17
        assert [m for m in groups if len(set(m))==2]==[marks]
        assert subtotal==(-Fraction(65542158940903,19982114331379200000000),Fraction(53663,19666196550000))
        assert tuple(x*10**20 for x in interval(subtotal,20))==(223747,223748)
    print('Single-double subtotal',marks,*subtotal,flush=True)
