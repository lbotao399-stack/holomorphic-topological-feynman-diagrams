"""All hard color orders on the triple/triple two-star prism occupancy."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json

from flint import fmpz_mpoly_ctx
from fourfold_quartic_prism import project_squares
from mixed_bare_primitive_pole import PrimitiveChamber, WAVES, hard_resistances, image, numerator, star_period, tree_flow
from quartic_seven_plane_prism import cut_weights
from quartic_seven_plane_single import cuts


MARKS = (3,3,3,7,7,7)


def labeled_period(probe):
    """Accumulate the original signed voltage on all 720 distinct hard orders."""
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    times = ctx.gens()
    paths = []
    for e in range(9):
        start = sum(map(len,paths))
        paths.append(tuple(range(start,start+MARKS.count(e)+1)))
    lengths = [sum(times[j] for j in p) for p in paths]
    cycles = probe.cycles
    weight = sum((cycles[1][p]*cycles[3][q]-cycles[1][q]*cycles[3][p])**2*lengths[p]*lengths[q]
                 for p,q in combinations(range(9),2))
    data = probe.graph,ctx,None,lengths,None,paths,None
    grouped = {}
    count = 0
    for order in permutations(range(6)):
        source = [(0,)*4]*6+[WAVES[k] for k in order]
        sheared = [(a,b+c,d,c) for a,b,c,d in source]
        det = numerator(probe.graph,sheared,tree_flow(probe.graph,sheared),[(3,1),(3,3)])
        if not det:
            continue
        count += 1
        content = int(det.content())*(1 if det.leading_coefficient()>0 else -1)
        primitive = det//content
        key = tuple((tuple(map(int,m)),int(k)) for m,k in primitive.terms())
        if key not in grouped:
            grouped[key] = primitive,defaultdict(int)
        _,coeff = grouped[key]
        plus,minus = order.index(0),order.index(1)
        for i,charge in enumerate(wave[3] for wave in source[6:]):
            if i!=plus:
                coeff[tuple(sorted((i,plus)))] -= charge*content**2
            if i!=minus:
                coeff[tuple(sorted((i,minus)))] += charge*content**2
    assert count==432 and len(grouped)==36
    dist = hard_resistances(tuple(probe.edges))
    transformed = None
    for det,coeff in grouped.values():
        voltage = sum(k*dist[6+i][6+j] for (i,j),k in coeff.items())
        if not voltage:
            continue
        v = ctx.from_dict({m:int(k) for m,k in voltage.terms()})
        # Six hard path placements / twelve core symmetries, twice for the soft
        # labels, and one half in the resistance expression for the raw voltage.
        term = project_squares(data,[(det,weight)],v,2)
        transformed = term if transformed is None else transformed+term
    assert len(transformed)==140164
    bank = cut_weights(transformed)
    assert len(bank)==741
    return cuts(bank)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--integrate',action='store_true',help='Reintegrate each color representative')
    parser.add_argument('--labeled',action='store_true',help='Recompute the independent 720-order numerator contraction')
    args = parser.parse_args()
    stored = json.loads(Path(__file__).with_suffix('.json').read_text())
    probe = PrimitiveChamber(3,MARKS)
    stabilizer = [a for a in probe.actions if sorted(a[h][0] for h in MARKS)==list(MARKS)]
    assert len(stabilizer)==2 and len(probe.actions)==12
    orbits = {}
    for order in permutations(range(6)):
        color = tuple((0,1,1,2,2,0)[i] for i in order)
        config = tuple(tuple(color[i] for i in group) for group in probe.groups)
        words = {sum(image(config,g),()) for g in stabilizer}
        key = ''.join(map(str,min(words)))
        orbits[key] = len(words)
    rows = stored['chambers']
    assert orbits=={color:size for color,size,_,_ in rows} and len(rows)==48
    assert 6*8*sum(orbits.values())==4320
    partial = sum(F(1,k**3) for k in range(1,1025))
    lo,hi = partial+F(1,2*1025**2),partial+F(1,2*1024**2)
    total = (F(),F())
    signs = Counter()
    for word,_,r,s in rows:
        value = F(r),F(s)
        if args.integrate:
            assert star_period(tuple(map(int,word)))==value
        a,b = value
        bounds = sorted((a+b*lo,a+b*hi))
        if bounds==[0,0]:
            signs[0] += 1
        elif bounds[0]>0:
            signs[1] += 1
        else:
            assert bounds[1]<0
            signs[-1] += 1
        total = tuple(x+y for x,y in zip(total,value))
    assert total==tuple(map(F,stored['sum'])) and signs=={1:20,-1:10,0:18}
    a,b = total
    assert 72<(a+b*lo)*10**16 < (a+b*hi)*10**16 <73
    if args.labeled:
        assert labeled_period(probe)==total
    print('Rational and zeta(3) coefficients:',*total)


if __name__=='__main__':
    main()
