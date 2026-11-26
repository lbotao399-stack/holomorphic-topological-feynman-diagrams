"""Prism hard placements admitting two disjoint source-free three-edge stars."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from fractions import Fraction as F
from pathlib import Path
from math import prod
import argparse
import json

from flint import fmpq, fmpz_mpoly_ctx
from mixed_bare_primitive_pole import METRIC, PRIME, WAVES, PrimitiveChamber, hard_resistances, image, numerator, residue
from fourfold_quartic_prism import project_slices, project_squares
from quartic_seven_plane_prism import cut_weights
from quartic_seven_plane_prism_check import dualize, kirchhoff, pair_graph
from quartic_seven_plane_single import cuts
from mixed_bare_ranks import witness_counts


def placements():
    probe = PrimitiveChamber(3)
    core, actions = probe.core, probe.actions
    pairs = set(combinations(range(6),2))-set(core)
    empty_complements = [{e for e,ends in enumerate(core) if u not in ends and v not in ends} for u,v in pairs]
    assert len(pairs)==6 and {3,5,7} in empty_complements
    feasible = {m for m in combinations_with_replacement(range(9),6)
                if any(set(m)<=allowed for allowed in empty_complements)}
    assert len(feasible)==159
    while feasible:
        marks = min(m for m in feasible if set(m)<={3,5,7})
        layouts = {tuple(sorted(g[e][0] for e in marks)) for g in actions}
        feasible -= layouts
        stabilizer = [g for g in actions if sorted(g[e][0] for e in marks)==list(marks)]
        assert len(stabilizer)*len(layouts)==12
        colors = set(permutations((0,0,1,1,2,2)))
        rows = []
        while colors:
            color = min(colors)
            cfg = tuple(tuple(color[i] for i,h in enumerate(marks) if h==e) for e in range(9))
            words = {sum(image(cfg,g),()) for g in stabilizer}
            colors -= words
            rows.append((marks,color,len(words),F(8*len(words)*len(layouts),12)))
        assert sum(n for _,_,n,_ in rows)==90
        yield marks, len(layouts), rows


def classify(check_zeros):
    expected = {(3,2,1):(54,54), (3,1,2):(54,54), (3,0,3):(30,54),
                (2,3,1):(54,54), (2,2,2):(38,72)}
    total = Counter()
    output = {}
    for marks,layouts,rows in placements():
        occ = tuple(marks.count(i) for i in (3,5,7))
        unresolved = []
        probe = PrimitiveChamber(3,marks)
        count = sum(witness_counts(probe.core,rows,unresolved))
        missing = {c for _,c,_,_ in unresolved}
        words = sum(n for _,c,n,_ in rows if c not in missing)
        assert (count,words)==expected.get(occ,(0,0))
        if check_zeros:
            for _,color,_,_ in unresolved:
                p = PrimitiveChamber(3,marks,color)
                assert all(not numerator(p.graph,p.sheared,p.flow,[(3,i),(3,j)])
                           for i,j in combinations(range(4),2))
        total.update(dict(layouts=layouts,shapes=1,nonzero=count,zero=len(unresolved),
                          labeled=8*layouts*words))
        output[occ] = (layouts,len(rows),count,words)
    assert total['layouts']==159 and total['shapes']==16 and total['nonzero']==230
    return output,total


@cache
def prism_data(marks,core_index=3):
    probe = PrimitiveChamber(core_index,marks)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    times = ctx.gens()
    paths = []
    for e in range(9):
        start = sum(map(len,paths))
        paths.append(tuple(range(start,start+marks.count(e)+1)))
    lengths = [sum(times[j] for j in path) for path in paths]
    return (probe.graph,ctx,None,lengths,None,paths,None),soft_polynomials(probe.cycles,lengths,ctx)


def soft_polynomials(cycle,lengths,ctx):
    soft = {}
    for p,q in combinations(range(9),2):
        coeffs = [(i,j,cycle[i][p]*cycle[j][q]-cycle[i][q]*cycle[j][p]) for i,j in combinations(range(4),2)]
        coeffs = tuple(triple for triple in coeffs if triple[-1])
        if not coeffs:
            continue
        if coeffs[0][-1]<0:
            coeffs = tuple((i,j,-k) for i,j,k in coeffs)
        soft[coeffs] = soft.get(coeffs,ctx.constant(0)) + lengths[p]*lengths[q]
    return soft


def square_roots(probe, soft, determinant=None):
    mins = {(i,j):numerator(probe.graph,probe.sheared,probe.flow,[(3,i),(3,j)],determinant=determinant)
            for i,j in combinations(range(4),2)}
    yield from normalized_roots(mins,soft)


def normalized_roots(mins,soft):
    for coeffs in soft:
        d = sum(k*mins[i,j] for i,j,k in coeffs)
        if d:
            c = int(d.content())*(1 if d.leading_coefficient()>0 else -1)
            yield d//c,coeffs,c*c


def chamber_period(marks,color,sliced=False):
    data,soft = prism_data(marks)
    ctx = data[1]
    probe = PrimitiveChamber(3,marks,color)
    groups = []
    for polynomial,coeffs,scale in square_roots(probe,soft):
        for pair in groups:
            if pair[0]==polynomial:
                pair[1] += scale*soft[coeffs]
                break
        else:
            groups.append([polynomial,scale*soft[coeffs]])
    if not groups:
        return F(),F()
    times = ctx.gens()
    voltage = sum(c*prod(times[i] for i in ex) for ex,c in probe.voltage_terms)
    new,jac = dualize(METRIC)
    new_u = kirchhoff(pair_graph(),new)
    values = tuple(fmpq(str(t)) for t in new)
    slices = project_slices(data,groups,voltage,1) if sliced else (project_squares(data,groups,voltage,1),)
    bank = defaultdict(F)
    at_metric = F()
    for poly in slices:
        at_metric += F(str(poly(*values)))
        for ex,c in cut_weights(poly).items():
            bank[ex] += c
    check = residue(at_metric*probe.weight/2/new_u**13)
    assert check == probe.numerator(METRIC,check=True)*pow(probe.kirchhoff(METRIC),-13,PRIME)*residue(jac)%PRIME
    return tuple(v*probe.weight/2 for v in cuts(bank))


def exact_labeled_sum(marks):
    """Collect the raw signed resistances for all 720 hard orders before projecting."""
    data,soft = prism_data(marks)
    ctx = data[1]
    orders = defaultdict(list)
    for order in permutations(range(6)):
        orders[tuple((0,1,1,2,2,0)[i] for i in order)].append(order)
    grouped = defaultdict(list)
    for color,group in orders.items():
        probe = PrimitiveChamber(3,marks,color)
        voltages = Counter()
        for order in group:
            plus,minus = order.index(0),order.index(1)
            for i,lab in enumerate(order):
                charge = WAVES[lab][3]
                if plus!=i:
                    voltages[tuple(sorted((i,plus)))] -= charge
                if minus!=i:
                    voltages[tuple(sorted((i,minus)))] += charge
        for d,rays,scale in square_roots(probe,soft):
            # The two evaluation keys only route comparisons. Polynomial equality
            # fixes the bucket and the content was removed before taking the square.
            key = (len(d),int(d(*range(1,16))),int(d(*range(3,18))))
            bucket = grouped[key]
            for existing,coeff in bucket:
                if d==existing:
                    break
            else:
                coeff = Counter()
                bucket.append((d,coeff))
            for pair,k in voltages.items():
                coeff[rays,pair] += scale*k
    probe = PrimitiveChamber(3,marks)
    layouts = len({tuple(sorted(g[h][0] for h in marks)) for g in probe.actions})
    dist = hard_resistances(tuple(probe.edges))
    bank = defaultdict(F)
    for bucket in grouped.values():
        for d,coeff in bucket:
            weight = ctx.constant(0)
            for rays,times in soft.items():
                v = sum(k*dist[6+i][6+j] for (r,(i,j)),k in coeff.items() if r==rays)
                if v:
                    weight += ctx.from_dict({m:int(c) for m,c in v.terms()})*times
            if weight:
                transformed = project_squares(data,[(d,weight)],ctx.constant(1),1)
                for ex,c in cut_weights(transformed).items():
                    bank[ex] += F(c)*F(layouts,len(probe.actions))
    return cuts(bank),sum(map(len,grouped.values())),len(bank)


def signed_group_counts():
    """Three nonsingular metrics separate signed cofactors up to rational scale."""
    from random import Random
    import numpy as np
    from flint import nmod_mat
    from mixed_bare_ranks import hard_waves, modular_network
    prime = 1000003
    units = np.zeros((17,4),dtype=np.int64)
    units[-4:] = np.eye(4,dtype=np.int64)
    for marks,_,rows in placements():
        occ = tuple(marks.count(i) for i in (3,5,7))
        if occ not in ((3,1,2),(2,2,2)):
            continue
        probe = PrimitiveChamber(3,marks)
        networks = [modular_network(probe.core,marks,probe.cycles,Random(k).sample(range(2,300),15),prime)
                    for k in (103,187,311)]
        across = Counter()
        for _,color,_,_ in rows:
            _,waves = hard_waves(color)
            sheared = np.column_stack((waves[:,0],waves[:,1]+waves[:,2],waves[:,3],waves[:,2]))
            minors = []
            for flow,_,cycles in networks:
                rays = flow@sheared%prime
                hard = np.vstack((np.ones((1,15),dtype=np.int64),
                                  (rays.T[:,None,:]*cycles[None,:,:]).reshape(16,15)))
                minors.append({pair:int(nmod_mat((np.column_stack((hard,units[:,pair]))%prime).tolist(),prime).det())
                               for pair in combinations(range(4),2)})
            local = set()
            for coeffs in prism_data(marks)[1]:
                v = tuple(sum(k*d[a,b] for a,b,k in coeffs)%prime for d in minors)
                pivot = next((a for a in v if a),None)
                if pivot:
                    local.add(tuple(a*pow(pivot,-1,prime)%prime for a in v))
            across.update(local)
        assert set(across.values())=={1}
        assert len(across)=={(3,1,2):156,(2,2,2):207}[occ]
        print(occ,'distinct projective cofactor signatures:',len(across))


def integrated_shapes(options):
    stored = json.loads(Path(__file__).with_suffix('.json').read_text())
    all_shapes = {tuple(m.count(i) for i in (3,5,7)):(m,n,rows) for m,n,rows in placements()}
    for label,record in stored.items():
        if options.shape and options.shape!=label:
            continue
        marks,layouts,rows = all_shapes[tuple(map(int,label))]
        expected = [(w,n,F(a),F(b)) for w,n,a,b in record['chambers']]
        assert [(tuple(map(int,w)),n) for w,n,_,_ in expected] == [(c,n) for _,c,n,_ in rows]
        total = tuple(sum(r[i] for r in expected) for i in (2,3))
        assert total==tuple(map(F,record['sum']))
        assert 8*layouts*sum(r[1] for r in expected)==720*layouts
        if options.integrate:
            for w,_,a,b in expected:
                assert chamber_period(marks,tuple(map(int,w)))==(a,b)
        if options.labeled:
            answer,groups,orbits = exact_labeled_sum(marks)
            assert answer==total
            print(label,'primitive squared polynomials and tetrahedron moment orbits:',groups,orbits)
        print(label,'rational and zeta(3) coefficients:',*total)


def two_star_sum(reintegrate=False, selected=None):
    """Check each saved orbit against its stabilizer and add completed occupancies."""
    from quartic_prism_physical_ibp import zeta_interval

    # Check the signed resistance coefficients before choosing any graph metric.
    # The eight raw voltages for each word contain a factor 1/2 per resistance.
    voltages = defaultdict(Counter)
    for order in permutations(range(6)):
        word = tuple((0,1,1,2,2,0)[i] for i in order)
        plus,minus = order.index(0),order.index(1)
        for i,label in enumerate(order):
            if i!=plus:
                voltages[word][tuple(sorted((i,plus)))] -= WAVES[label][3]
            if i!=minus:
                voltages[word][tuple(sorted((i,minus)))] += WAVES[label][3]
    assert len(voltages)==90
    for word,coefficients in voltages.items():
        expected = {}
        for i,j in combinations(range(6),2):
            colors = tuple(sorted((word[i],word[j])))
            expected[i,j] = 4*{(0,0):4,(1,2):1,(0,1):-1,(0,2):-1}.get(colors,0)
        assert coefficients == Counter(expected)

    here = Path(__file__)
    records = json.loads(here.with_suffix('.json').read_text())
    records['303'] = json.loads(here.with_name('mixed_bare_prism_colors.json').read_text())
    for name in ('312','222'):
        saved = json.loads(here.with_name(f'{here.stem}_{name}.json').read_text())
        assert saved['occupancy'] == name and saved['completed'] == len(saved['chambers'])
        records[name] = {'chambers': [(w,*row) for w,row in saved['chambers'].items()],
                         'sum': saved['partial_sum'], 'total': saved['total']}
    if selected is not None and selected not in records:
        raise ValueError(f'Unknown two-star occupancy {selected}')

    zeta = zeta_interval(3)
    complete = []
    for marks,layouts,rows in placements():
        name = ''.join(str(marks.count(i)) for i in (3,5,7))
        if name not in records:
            continue
        record = records[name]
        expected = {''.join(map(str,w)):n for _,w,n,_ in rows}
        assert len(expected)==record.get('total',len(rows))
        seen,signs = set(),Counter()
        a = b = F()
        for word,n,x,y in record['chambers']:
            assert word not in seen and expected[word]==n
            seen.add(word)
            r,s = F(x),F(y)
            if reintegrate and selected in (None,name):
                assert chamber_period(marks,tuple(map(int,word)),sliced=True)==(r,s)
            lo,hi = sorted(r+s*v for v in zeta)
            if lo == hi == 0:
                signs[0] += 1
            elif lo > 0:
                signs[1] += 1
            else:
                assert hi < 0
                signs[-1] += 1
            a += r
            b += s
        assert (a,b)==tuple(map(F,record['sum']))
        print(name,f'{len(seen)}/{len(expected)}',dict(sorted(signs.items())),a,b)
        if seen == expected.keys():
            assert sum(expected.values())==90
            labeled = 8*layouts*sum(n for w,n,x,y in record['chambers'] if F(x) or F(y))
            complete.append((a,b,labeled,signs[1]+signs[-1]))
    if len(complete)==5:
        a,b,labeled,nonzero = map(sum,zip(*complete))
        assert labeled==21600 and nonzero==230
        lower,upper = sorted(a+b*v for v in zeta)
        assert -428900592 < 10**20*lower <= 10**20*upper < -428900591
        assert -80418861 < 10**20*F(3,16)*lower <= 10**20*F(3,16)*upper < -80418860
        print('All two-star chambers:',a,b,'with nonzero labeled orders:',labeled)
        return a,b


if __name__=='__main__':
    args = argparse.ArgumentParser()
    args.add_argument('--zeros',action='store_true',help='Prove each witness-free numerator vanishes symbolically')
    args.add_argument('--signed-groups',action='store_true',help='Test proportionality before integrating the two dense populations')
    args.add_argument('--sliced-check',action='store_true',help='Check the marked-degree slices on one dense (3,1,2) chamber')
    args.add_argument('--two-star-sum',action='store_true',help='Check all five signed color tables and sum complete two-star periods')
    args.add_argument('--shape',help='Select a recorded occupancy, such as 231')
    args.add_argument('--integrate',action='store_true',help='Reintegrate the individual color-orbit contributions')
    args.add_argument('--labeled',action='store_true',help='Recontract the 720 original signed voltages before integration')
    options = args.parse_args()
    if options.two_star_sum:
        two_star_sum(options.integrate,options.shape)
        raise SystemExit
    table,total = classify(options.zeros)
    print('Nonzero and vanishing colored chamber orbits:',total['nonzero'],total['zero'])
    integrated_shapes(options)
    if options.signed_groups:
        signed_group_counts()
    if options.sliced_check:
        small = (3,3,5,5,5,7)
        data,soft = prism_data(small)
        probe = PrimitiveChamber(3,small,(0,0,1,1,2,2))
        groups = [(d,v*soft[r]) for d,r,v in square_roots(probe,soft)]
        times = data[1].gens()
        voltage = sum(k*prod(times[i] for i in ex) for ex,k in probe.voltage_terms)
        columns = list(project_slices(data,groups,voltage,1))
        assert len(columns)==35 and sum(map(len,columns))==53354
        assert sum(columns)==project_squares(data,groups,voltage,1)
        marks = (3,3,3,5,7,7)
        value = chamber_period(marks,(0,2,1,0,1,2),sliced=True)
        assert value == (-F(674063179358209292219,37641028122134928691200000000000),
                         F(7747407785891,518642913940351200000000))
        print('Sliced (3,1,2) chamber:',*value)
