"""Zero-sum line placements and hard-label orbits in the eight-input fourfold Ward equation."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import lcm, prod
from pathlib import Path
from random import Random
import json
import sys

from flint import fmpq, fmpz_mpoly_ctx, nmod_mat
import numpy as np
from fourfold_quartic_ward import NORMAL, collinear_square
from quartic_seven_plane_bare import dot, minors
from quartic_seven_plane_double import CORE, graph_data, image, local_moment, tree_flow
from quartic_seven_plane_double_three_two import integrate_polynomials
from quartic_seven_plane_gale import gale
from quartic_seven_plane_lift import fixture
from quartic_seven_plane_primitive import admissible
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from current_minor import numerator
from quartic_planar_cores import cubic_cores
from square_three import G


def oriented_symmetries(core):
    destinations = defaultdict(list)
    for i,edge in enumerate(core):
        destinations[edge].append(i)
    for vertex in permutations(range(6)):
        pairs = [(vertex[a],vertex[b]) for a,b in core]
        ordered = [tuple(sorted(p)) for p in pairs]
        if Counter(ordered) != Counter(core):
            continue
        for options in product(*(permutations(destinations[k]) for k in destinations)):
            targets = {k:iter(p) for k,p in zip(destinations,options)}
            yield tuple((next(targets[k]),a!=k) for k,a in zip(ordered,pairs))


def stabilizer(core, marks):
    return [action for action in oriented_symmetries(core)
            if tuple(sorted(action[h][0] for h in marks)) == marks]


def label_images(order, marks, subgroup):
    cfg = tuple(tuple(order[i] for i,h in enumerate(marks) if h==j) for j in range(9))
    return {sum(image(cfg,a),()) for a in subgroup}


def hard_orbits(core, marks):
    group = stabilizer(core,marks)
    orders = set(permutations(range(5)))
    while orders:
        order = min(orders)
        images = label_images(order,marks,group)
        orders -= images
        yield order, F(len(images),len(group))


def voltage_groups(orders):
    """Collect equal normal volume data; retain the precise weighted hard-slot voltage."""
    coordinates = gale(NORMAL)+[0]
    groups = {}
    for order,weight in orders:
        key = tuple(coordinates[i] for i in order)
        if key not in groups:
            groups[key] = (order,[[F() for _ in range(5)] for _ in range(5)])
        weights = groups[key][1]
        for sign,k in ((1,0),(-1,1)):
            for j,h in enumerate(order):
                weights[order.index(k)][j] += weight*sign*NORMAL[h][2]
    return sorted(groups.items())


def grouped_orbits(core,marks):
    # Pick the same Gale ordering for geometric images before combining their voltages.
    z = gale(NORMAL)+[0]
    group = stabilizer(core,marks)
    ordered = []
    for order,weight in hard_orbits(core,marks):
        representative = min(label_images(order,marks,group),key=lambda p:(tuple(z[i] for i in p),p))
        ordered.append((representative,weight))
    return voltage_groups(ordered)


def endpoint_polynomials(graph, order):
    # The four-component auxiliary source has one longitudinal unit dipole.
    normal = [(0,0,0)]*6 + [NORMAL[i] for i in order]
    ctx = fmpz_mpoly_ctx.get(tuple('t'+str(i) for i in range(14)))
    points = [[ctx.constant(0)]*4]
    for v in range(1,11):
        source = [(int(i==v)-int(i==0),*q) for i,q in enumerate(normal)]
        flow = tree_flow(graph,source)
        points.append([(-1)**r*numerator(graph,source,flow,[(0,j) for j in range(4) if j!=r])
                       for r in range(4)])
    return points


def squared_cofactor_terms(graph, points):
    """Yield factors in the 14-edge numerator before division by U^8."""
    ctx = points[0][0].context()
    t = ctx.gens()
    columns = list(zip(*graph.C))
    for triple in combinations(range(14),3):
        k = minors([columns[i] for i in triple])
        if not any(k):
            continue
        means, differences = [], []
        for i in triple:
            a,b = graph.edges[i]
            left,right = dot(k,points[a]),dot(k,points[b])
            means.append(left+right)
            differences.append(left-right)
        yield 2*prod(t[i] for i in triple), [(3,x-y) for x,y in combinations(means,2)] + [(2,x) for x in differences]


def squared_cofactor_polynomial(graph,order):
    points = endpoint_polynomials(graph,order)
    total = points[0][0].context().constant(0)
    for lengths,terms in squared_cofactor_terms(graph,points):
        total += lengths*sum(w*p**2 for w,p in terms if p)
    return total


def separated_squares(graph,order):
    """Collect identical primitive factors before squaring any 14-variable polynomial."""
    points = endpoint_polynomials(graph,order)
    samples = ([2+(i*i+3*i)%17 for i in range(14)], [1+(5*i*i+3*i)%29 for i in range(14)])
    buckets = defaultdict(list)
    for lengths,terms in squared_cofactor_terms(graph,points):
        for weight,p in terms:
            if not p:
                continue
            content = int(p.content())*(1 if p.leading_coefficient()>0 else -1)
            primitive = p//content
            key = (len(primitive),*(int(primitive(*t)) for t in samples))
            coefficient = weight*content*content*lengths
            for entry in buckets[key]:
                if entry[0] == primitive:
                    entry[1] += coefficient
                    break
            else:
                buckets[key].append([primitive,coefficient])
    return [entry for bucket in buckets.values() for entry in bucket]


def integrated_group(data,order,weights,check=False,separate=False):
    graph,ctx,u,_,dist,_,_ = data
    voltage = -sum(fmpq(str(weights[i][j]))*dist[6+i][6+j] for i in range(5) for j in range(5)
                   if weights[i][j])/2
    if not voltage:
        return F()
    den = lcm(*(int(v.denominator) for _,v in voltage.terms()))
    voltage = ctx.from_dict({m:int(c*den) for m,c in voltage.terms()})
    zero = ctx.constant(0)
    if separate:
        gram = [[1,0,0,0],[0]*4,[0]*4,[0]*4]
        reduced = (*data[:3],gram,*data[4:])
        total = F()
        for p,weight in separated_squares(graph,order):
            total += integrate_polynomials(reduced,[p,zero,zero,zero],weight*voltage)
            if local_moment.cache_info().currsize > 25000:
                local_moment.cache_clear()
        local_moment.cache_clear()
        return total/den
    poly = squared_cofactor_polynomial(graph,order)
    if check:
        times = [2+(i*i+3*i)%17 for i in range(14)]
        direct = collinear_square(graph.edges,graph.C,times,[(0,0,0)]*6+[NORMAL[i] for i in order],(1,1,-2))
        assert F(int(poly(*times)),int(u(*times))**8) == direct
    if not poly:
        return F()
    # The existing separated-arm integrator takes a quadratic form in four polynomials.
    # Put B and 1 in the two off-diagonal entries to integrate 2 B times the voltage.
    gram = [[0,1,0,0],[1,0,0,0],[0]*4,[0]*4]
    one,zero = ctx.constant(1),ctx.constant(0)
    result = integrate_polynomials((*data[:3],gram,*data[4:]),[poly,one,zero,zero],voltage)/(2*den)
    local_moment.cache_clear()
    return result


def recorded_shapes():
    data = json.loads(Path(__file__).with_name('fourfold_quartic_periods.json').read_text())['two_double']
    assert [tuple(r['hard_paths']) for r in data['shapes']] == [m for m,_ in admissible(CORE)]
    for row in data['shapes']:
        groups = grouped_orbits(CORE,tuple(row['hard_paths']))
        assert [k for k,_ in groups] == [tuple(k) for k,_ in row['gale_groups']]
        assert sum(F(x) for _,x in row['gale_groups']) == F(row['total']) < 0
    assert sum(F(row['total']) for row in data['shapes']) == F(data['total']) == -F(950396813,399243043584000000)
    return data


def integrate_shape(index):
    marks,_ = list(admissible(CORE))[index]
    data = graph_data(marks)
    total = F()
    recorded = recorded_shapes()['shapes'][index]
    for (key,(order,weights)),(_,expected) in zip(grouped_orbits(CORE,marks),recorded['gale_groups'],strict=True):
        value = integrated_group(data,order,weights,separate=index>=2)
        assert value == F(expected)
        total += value
        print('group',index,key,value,flush=True)
    return total


def census():
    observed = []
    for core in cubic_cores():
        shapes = list(admissible(core))
        sizes = Counter()
        for marks,_ in shapes:
            sizes.update(weight.denominator for _,weight in hard_orbits(core,marks))
        observed.append((len(shapes),sum(s for _,s in shapes),dict(sizes)))
    assert observed == [(11,69,{1:510,2:15}),(0,0,{}),(36,171,{1:2550,2:30}),
                        (80,786,{1:7830,2:60}),(26,1134,{1:1890})]
    return observed


def metric_sum(core,marks,times,groups,check=False):
    """All normal labelings at one metric, modulo 10007, with three line insertions integrated."""
    prime = 10007

    def field(a):
        return nmod_mat(a.tolist(),prime)

    def array(m):
        return np.asarray(m.tolist(),dtype=np.int64)

    edges = fixture(core,marks)
    left,right = (np.asarray(x) for x in zip(*edges))
    cycles = np.asarray(G(11,edges).C,dtype=np.int64)
    incidence = np.asarray([[int(x==k)-int(y==k) for k in range(1,11)] for x,y in edges])
    t = np.asarray(times,dtype=np.int64)
    current = np.asarray([pow(int(x),-1,prime) for x in t])[:,None]*incidence
    green = array(field(incidence.T@current % prime).inv())
    flow = current@green % prime
    tree = int(field((cycles*t)@cycles.T % prime).det())
    triples,wedges = [],[]
    for triple in combinations(range(14),3):
        k = minors([tuple(v) for v in cycles.T[list(triple)].tolist()])
        if any(k):
            triples.append(triple)
            wedges.append(k)
    triples,wedges = np.asarray(triples),np.asarray(wedges)
    lengths = t[triples].prod(axis=1) % prime
    total = 0
    for _,(order,weights) in groups:
        voltage = int((np.asarray(weights,dtype=np.int64)*green[5:,5:]).sum()) % prime
        if not voltage:
            continue
        hard = flow[:,5:] @ np.asarray([NORMAL[i] for i in order]) % prime
        rows = np.vstack([np.ones(14,dtype=np.int64),
                          (hard.T[:,None,:]*cycles[None,:,:]).reshape(12,14)]) % prime
        kernel,dim = field(rows).nullspace()
        if dim != 1:
            continue
        k = array(kernel)[:,0]
        i = int(np.flatnonzero(k)[0])
        scale = int(field(np.delete(rows,i,axis=1)).det()) * pow(int(k[i]),-1,prime) % prime
        k = k*scale % prime
        potentials = np.vstack([np.zeros(4,dtype=np.int64),flow.T@(k[:,None]*cycles.T) % prime])
        start,end = potentials[left[triples]],potentials[right[triples]]
        means = np.einsum('ijk,ik->ij', start+end,wedges) % prime
        diff = np.einsum('ijk,ik->ij', start-end,wedges) % prime
        variance = (3*((means-np.roll(means,1,axis=1))**2).sum(axis=1)+2*(diff*diff).sum(axis=1)) % prime
        square = 2*int((lengths*variance).sum())*pow(sum(times),-2,prime) % prime
        if check:
            direct = collinear_square(edges,cycles.tolist(),times,[(0,0,0)]*6+[NORMAL[i] for i in order],
                                      (1,1,-2))
            assert square == direct.numerator*pow(direct.denominator,-1,prime) % prime
        total += square*voltage
    return total*pow(tree,-4,prime) % prime


def check_every_shape():
    rng = Random(7561)
    groups = voltage_groups((p,1) for p in permutations(range(5)))
    assert len(groups) == 30
    counts = []
    for j,core in enumerate(cubic_cores()):
        shapes = list(admissible(core))
        for i,(marks,_) in enumerate(shapes):
            for _ in range(3):
                times = [rng.randrange(1,39) for _ in range(14)]
                try:
                    if metric_sum(core,marks,times,groups,check=(i==5 and j==0)):
                        break
                except ZeroDivisionError:
                    continue
            else:
                raise ArithmeticError(('no nonzero modular witness',j,marks))
        counts.append(len(shapes))
    assert counts == [11,0,36,80,26]
    return counts


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--shape',type=int,choices=range(11))
    args = parser.parse_args()
    if args.shape is not None:
        value = integrate_shape(args.shape)
        print('shape',args.shape,value)
    else:
        print('Hard path shapes, path multisets, labeled orbits by stabilizer',census())
        print('Nonzero rational-integrand witnesses by core',check_every_shape())
        print('Two-double integrated total', recorded_shapes()['total'])
        data = graph_data((0,0,0,7,7))
        _,(order,weights) = voltage_groups([((0,4,1,2,3),1)])[0]
        value = integrated_group(data,order,weights,True)
        assert value == integrated_group(data,order,weights,separate=True) == -F(53,70413235200000)
        print('Complete three-line placement integral at one hard order',value)
