"""Recover the four-loop elementary probes and D3 contact from the five- and six-mark moments."""
from fractions import Fraction as F
import json
from pathlib import Path
from math import ceil, floor
import sys

from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
IIB = HERE.parent/'iib-bcov'
sys.path.insert(0,str(IIB))
from quartic_prism_physical_ibp import ordinary_zeta_interval
from quartic_scalar_response import BIPARTITE, total_channels
from worldsheet_primary import restricted_symbols


def elementary_period_weights(plane):
    """Coefficients of K4, J4 and W for seven two-component scalar momenta."""
    q = tuple(tuple(map(F,p)) for p in plane)
    if len(q)!=7 or any(len(p)!=2 for p in q):
        raise ValueError('Expected seven pairs of rational coordinates in the e1,e2 plane.')
    a,b,d = restricted_symbols(q)
    ambiguity = (22*a-5*b+18*d)/288
    return b/576+3*ambiguity/8,(8*a-b)/9+8*ambiguity,24*ambiguity


def elementary_four_loop(plane):
    """Return (1,zeta(3),zeta(5)) coefficients of the heat-stripped Q4 phi(0).

    The absorbed momentum is e4, the vector at e3 has polarization e_4, and
    the two additional scalar waves are e5 and the momentum balancing all slots.
    """
    weights = elementary_period_weights(plane)
    data = json.loads((HERE/'four_loop_arakelov_periods.json').read_text())
    if 'equal_strength_W' not in data:
        raise ValueError('The shared genus-four table has not yet been summed.')
    return tuple(weights[0]*k+weights[1]*j+weights[2]*F(w) for (k,j),w in
                 zip(total_channels(BIPARTITE),data['equal_strength_W']))


def direct_interval(coefficients, bounds, scale):
    lo = hi = coefficients[0]
    for c,(a,b) in zip(coefficients[1:],bounds):
        lo += min(c*a,c*b)
        hi += max(c*a,c*b)
    return floor(lo*scale),ceil(hi*scale)


def check_completed_graphs():
    data = json.loads((HERE/'four_loop_arakelov_periods.json').read_text())
    descendant = json.loads((HERE/'four_loop_descendant.json').read_text())
    expected = {'two-double':3,'one-double':8,'prism':28,'bipartite':11}
    assert data['checked_graph_count']==50 and data['expected_graphs']==expected
    rows = []
    for core,count in expected.items():
        assert set(data['graphs'][core])==set(map(str,range(count)))
        rows.extend(tuple(map(F,v)) for v in data['graphs'][core].values())
    auxiliary = tuple(map(sum,zip(*rows)))
    assert auxiliary == tuple(map(F,data['sum']))
    metric = total_channels(BIPARTITE)
    five,six = [],[]
    for a,d,(k,j) in zip(auxiliary,map(F,descendant['sum']),metric):
        curvature = a-4*d-2*j
        five.append(curvature/8-k/1152)
        six.append(11*curvature/12+j/3)
    w = tuple(18*i-3*h for i,h in zip(five,six))
    centered = tuple(12*h for h in six)
    axial = tuple(144*i for i in five)
    for value,name in ((w,'equal_strength_W'),(centered,'centered_V_1_1_minus_2'),
                       (axial,'stripped_axial_pole_72C4')):
        assert value == tuple(map(F,data[name]))
    u,v,t = (1,0),(0,1),(-4,-3)
    waves = (([u]*3+[v]*3+[t]),([u]*4+[v]*2+[(-4,-2)]),([u]*4+[v]*3))
    symbols = [restricted_symbols(p) for p in waves]
    assert symbols == [(3492,216,-1728),(288,0,0),(36,216,0)]
    assert [elementary_period_weights(p) for p in waves] == [
        (F(117,2),4320,3720),(F(33,4),432,528),(0,0,-24)]
    traces = [elementary_four_loop(p) for p in waves]
    assert traces[2]==tuple(-24*h for h in w)
    plane = (u,v,t,(-1,3),(2,-5),(F(2,3),1),(1,F(-1,7)))
    mapped = [(F(3,2)*x+F(5,7)*y,-F(2,9)*x+F(4,11)*y) for x,y in plane]
    factor = (F(3,2)*F(4,11)+F(5,7)*F(2,9))**8
    assert elementary_four_loop(mapped) == tuple(factor*x for x in elementary_four_loop(plane))
    assert elementary_four_loop([u]*7) == (0,0,0)
    assert all(3*p+t==6*q for p,t,q in zip(axial,traces[2],centered))
    matrix = fmpq_mat([[str(x) for x in row] for row in traces])
    prime = 999999937
    residue = lambda x: int(x.numerator)*pow(int(x.denominator),-1,prime)%prime
    determinant = residue(matrix.det())
    assert determinant==99665298
    print('Three elementary period rows modulo',prime,
          [[residue(x) for x in row] for row in traces],'determinant',determinant)
    inverse = matrix.inv()
    assert inverse*matrix==fmpq_mat([[int(i==j) for j in range(3)] for i in range(3)])
    bounds = [ordinary_zeta_interval(s) for s in (3,5)]
    signs = {}
    for core,entries in data['graphs'].items():
        intervals = [direct_interval(tuple(map(F,v)),bounds,10**24) for v in entries.values()]
        signs[core] = (sum(lo>0 for lo,hi in intervals),sum(hi<0 for lo,hi in intervals))
    assert signs == {core:(n,0) for core,n in expected.items()}
    print('Auxiliary graph signs (positive,negative)',signs)
    # The lower-genus polygon and theta give the other terms of this fourfold Ward component.
    tree_and_lower = (F(1057669573,1123411161600000),-F(143249,181621440000),F())
    bare = tuple(-252*h-x for h,x in zip(centered,tree_and_lower))
    for name,row in zip(('W','centered','72 C4','-24 W','bare centered Ward integral'),
                        (w,centered,axial,traces[2],bare)):
        enclosure=direct_interval(row,bounds,10**20)
        print(name,tuple(map(str,row)),'in',enclosure,'/10^20')
    assert direct_interval(traces[2],bounds,10**20)[1] < 0
    assert all(direct_interval(row,bounds,10**20)[0] > 0 for row in (five,six,w,axial))


if __name__=='__main__':
    check_completed_graphs()
