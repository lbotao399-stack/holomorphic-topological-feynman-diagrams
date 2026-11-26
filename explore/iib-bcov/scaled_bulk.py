"""Three-loop transverse bulk periods for two scaled pairs after a cotree reduction."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, combinations_with_replacement, product
from math import factorial, prod
from pathlib import Path

from fourfold_bulk import CORES, symmetries
from labeled_bulk import LabeledBulk


def source_and_tests(strength):
    # Cases 1 and 2 scale the second B slot, 3 lifts along the removed axis,
    # and 4 uses two independent projected B slots.
    b_scale = 2 if strength in (2,3) else 1
    axes = (0,0,0,1,1,2,2,3,4)
    p = [tuple(s*int(i==j) for j in range(5)) for i,s in
         zip(axes,(1,1,1,1,2,1,b_scale,1,1))]
    if strength==3:
        p[3:7] = [(i+1,*q[1:]) for i,q in enumerate(p[3:7])]
    if strength==4:
        p[6] = (0,1,1,0,0)
    p.append(tuple(-sum(c) for c in zip(*p)))
    def transfer(i,j,v):
        out = [[0]*5 for _ in p]
        out[i] = list(v)
        out[j] = [-x for x in v]
        assert all(sum(x*y for x,y in zip(vi,pi))==0 for vi,pi in zip(out,p))
        return out
    if strength==1:
        return p,[transfer(4,3,(0,0,1,0,0)),transfer(3,0,(0,0,1,0,0)),
                  transfer(3,9,(-2,0,3,0,0))]
    if strength==2:
        return p,[transfer(3,9,(-1,0,1,0,0))]
    if strength==3:
        v = [[0]*5 for _ in p]
        v[3],v[5],v[9] = [3,-3,0,0,0],[0,-10,0,0,0],[-3,13,0,0,0]
        assert all(sum(x*y for x,y in zip(q,s))==0 for q,s in zip(p,v)) and not any(map(sum,zip(*v)))
        return p,[v]
    if strength==4:
        return p,[transfer(5,0,(0,1,0,0,0)),transfer(3,9,(-2,0,3,0,0))]
    raise ValueError('Unknown momentum fixture.')


def placements(core,strength):
    ids = range(len(core))
    visited = set()
    perms = set(symmetries(core))
    for a in permutations(ids,2):
        b_choices = (product(ids,repeat=2) if strength==4 else
                     (combinations if strength==1 else permutations)(ids,2))
        for b in b_choices:
            for t in combinations_with_replacement(ids,3):
                if (a,b,t) in visited:
                    continue
                orbit = set()
                for perm in perms:
                    x,y,z = [tuple(perm[i] for i in group) for group in (a,b,t)]
                    orbit.add((x,tuple(sorted(y)) if strength==1 else y,tuple(sorted(z))))
                visited.update(orbit)
                labeled = (2 if strength==1 else 1)*6//prod(factorial(n) for n in Counter(t).values())
                yield a+b+t,F(len(orbit)*labeled*6,len(symmetries(core)))
    assert len(visited)=={1:25200,2:50400,3:50400,4:60480}[strength]


def geometry(core,placement,p):
    branch = max(max(e) for e in core)+1
    source = [(0,)*4]*branch+[q[1:] for q in p[3:]]
    edges = []
    for i,(a,b) in enumerate(core):
        path = [a,*[branch+j for j,k in enumerate(placement) if i==k],b]
        edges.extend(zip(path,path[1:]))
    return edges,source


def weights(p,polar):
    # A transverse transfer, modulo a common shift, is minus one half of the
    # resistance sum with coefficients p_i(v_j)+p_j(v_i).
    answer = Counter()
    for i,j in combinations(range(10),2):
        weight = -F(sum(x*y+s*t for x,y,s,t in zip(p[i],polar[j],p[j],polar[i])),2)
        if not weight:
            continue
        assert not (i<3 and j<3)
        choices = []
        for order in permutations((8,9,10)):
            # Four branch vertices, seven retained slots. The three deleted
            # vertices all have the marginal of a uniform cotree point.
            site = (-1,-1,-1,4,5,6,7,*order)
            choices.append(tuple(sorted((site[i],site[j]))))
        for pair in choices:
            answer[pair] += weight/6
    return {key:v for key,v in answer.items() if v}


def periods(core,strength,checkpoint=None):
    p,tests = source_and_tests(strength)
    rows = [weights(p,v) for v in tests]
    channels = [any(i==-1 for i,j in w) for w in rows]
    total = [F(0)]*(2*len(rows))
    saved = {}
    if checkpoint and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            label,*values = line.split()
            if label in saved or len(values)!=len(total):
                raise ValueError('Invalid checkpoint entry.')
            saved[label] = tuple(map(F,values))
    n = 0
    for placement,factor in placements(core,strength):
        label = ''.join(map(str,placement))
        n += 1
        if label in saved:
            value = saved.pop(label)
        else:
            edges,source = geometry(core,placement,p)
            value = [F(0)]*len(total)
            if len({tuple(sorted(e)) for e in edges})==len(edges):
                graph = LabeledBulk(edges,source,all_orders=True)
                if graph.moments.base:
                    value = [factor*x for row,ct in zip(rows,channels)
                             for x in graph.resistance_pairing(row,cotree=ct)]
            if checkpoint:
                with checkpoint.open('a') as f:
                    print(label,*value,sep='\t',file=f)
        for k,v in enumerate(value):
            total[k] += v
    if saved:
        raise ValueError('Checkpoint contains a placement outside this calculation.')
    return n,total


def trace_rows(strength,exe):
    import sys
    directory = Path(__file__).resolve().parent
    sys.path.insert(0,str(directory/'operator'))
    from three_loop_soft_jet import exact_boundary
    from bulk_three_hexagon import potentials
    p,tests = source_and_tests(strength)
    a = exact_boundary(exe,p)
    hexagon = potentials(p)
    dot = lambda x,y:sum(i*j for v,w in zip(x,y) for i,j in zip(v,w))
    return [[dot(symbol,v) for symbol in a] for v in tests],[dot(hexagon,v) for v in tests]


def restored_check(strength):
    from point_graph import solve_graph
    from scalar_fan import cofactor
    from square_three import EDGES,G
    from flint import fmpz_mat
    p,tests = source_and_tests(strength)
    placement = (0,2,4,3,1,5,4)
    core_edges,source = geometry(EDGES,placement,p)
    graph = G(len(source),core_edges)
    times = [2+(7*i*i+3*i)%11 for i in range(len(core_edges))]
    flow,_ = solve_graph(core_edges,source,times,0)
    base = cofactor(flow,graph.C,times)
    assert base
    cotrees = [s for s in combinations(range(len(times)),3)
               if fmpz_mat([[row[j] for j in s] for row in graph.C]).det()]
    for cut in (cotrees[0],cotrees[-1]):
        edges,split,ancestor = [],[],[]
        for j,(a,b) in enumerate(core_edges):
            if j in cut:
                v = len(source)+cut.index(j)
                edges += [(a,v),(v,b)]
                split += [F(times[j],2)]*2
                ancestor += [j,j]
            else:
                edges.append((a,b)); split.append(times[j]); ancestor.append(j)
        cycles = [[row[j] for j in ancestor] for row in graph.C]
        nv = len(source)+3
        # Exact grounded Green matrix on this restored five-dimensional graph.
        load = [[int(i==j) for j in range(1,nv)] for i in range(nv)]
        _,green = solve_graph(edges,load,split,0)
        def resistance(i,j):
            return green[i][i-1]+green[j][j-1]-2*green[i][j-1]
        expected = []
        for test in tests:
            total = F(0)
            for (i,j),w in weights(p,test).items():
                mean = (sum(resistance(l,j) for l in range(nv-3,nv))/3
                        if i==-1 else resistance(i,j))
                total += w*mean
            expected.append(total)
        actual = [F(0)]*len(tests)
        for t in permutations((8,9,10)):
            site = (11,12,13,4,5,6,7,*t)
            charge = [(0,)*5]*nv
            for v,q in zip(site,p):
                charge[v] = q
            flow,phi = solve_graph(edges,charge,split,0)
            assert cofactor(flow,cycles,split)==base
            for l in permutations((11,12,13)):
                site = (*l,4,5,6,7,*t)
                for k,test in enumerate(tests):
                    actual[k] += sum(x*y for v,q in zip(site,test) for x,y in zip(phi[v],q))/36
        assert actual==expected and any(actual)
    print('Restored five-dimensional cofactors and transverse voltages:',strength,flush=True)


if __name__=='__main__':
    cli = ArgumentParser()
    cli.add_argument('core',choices=('tetrahedron','double-square','boundary','restored'))
    cli.add_argument('strength',type=int,choices=(1,2,3,4))
    cli.add_argument('--checkpoint',type=Path)
    args = cli.parse_args()
    if args.core=='restored':
        restored_check(args.strength)
    elif args.core=='boundary':
        from subprocess import run
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as td:
            exe = Path(td)/'nine'
            cpp = Path(__file__).resolve().parent/'operator'/'nine_sextics.cpp'
            run(['c++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
            print(trace_rows(args.strength,exe),flush=True)
    else:
        print(args.core,args.strength,*periods(CORES[args.core],args.strength,args.checkpoint),flush=True)
