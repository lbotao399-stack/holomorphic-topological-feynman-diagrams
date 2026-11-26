"""One cotree reduction of transverse BCOV bulk potentials at two and three loops."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
from argparse import ArgumentParser

from labeled_bulk import LabeledBulk
from square_three import EDGES


THETA = ((0,1),)*3
DOUBLE_SQUARE = ((0,1),(0,1),(0,3),(1,2),(2,3),(2,3))
CORES = {'theta':THETA,'tetrahedron':EDGES,'double-square':DOUBLE_SQUARE}


@cache
def symmetries(core):
    nv = max(max(e) for e in core)+1
    actions = []
    for vertex in permutations(range(nv)):
        target = [tuple(sorted((vertex[a],vertex[b]))) for a,b in core]
        if Counter(target)!=Counter(core):
            continue
        opts = [tuple(permutations([i for i,e in enumerate(core) if e==pair]))
                for pair in sorted(set(core))]
        for blocks in product(*opts):
            mapping = {}
            for pair,block in zip(sorted(set(core)),blocks):
                mapping.update(zip([i for i,e in enumerate(target) if e==pair],block))
            actions.append(tuple(mapping[i] for i in range(len(core))))
    return tuple(actions)


def configurations(core):
    n = 2 if core==THETA else 3
    ids = range(len(core))
    visited = set()
    for a in combinations(ids,2):
        for b in combinations(ids,n-1):
            for t in combinations_with_replacement(ids,3):
                key = (a,b,t)
                if key in visited:
                    continue
                orbit = set()
                for perm in symmetries(core):
                    x,y,z = (tuple(sorted(perm[i] for i in q)) for q in key)
                    orbit.add((x,y,z))
                    if n==3:
                        orbit.add((y,x,z))
                visited.update(orbit)
                assignments = factorial(len(a))*factorial(len(b))*6//prod(
                    factorial(v) for v in Counter(t).values())
                yield key,F(len(orbit)*assignments*factorial(n),len(symmetries(core)))
    assert len(visited)==(12600 if n==3 else 90)


def subdivision(core,config):
    a,b,t = config
    nv = max(max(e) for e in core)+1
    source = [(0,)*4]*nv
    basis = [tuple(int(i==j) for i in range(4)) for j in range(4)]
    balance = (-len(a),-len(b),-1,-1)
    groups = [[],[],[]]
    edges = []
    for i,(left,right) in enumerate(core):
        path = [left]
        for group,(assigned,waves) in enumerate(zip(config,
                ([basis[0]]*len(a),[basis[1]]*len(b),[basis[2],basis[3],balance]))):
            for edge,wave in zip(assigned,waves):
                if edge==i:
                    v = len(source)
                    path.append(v)
                    groups[group].append(v)
                    source.append(wave)
        path.append(right)
        edges.extend(zip(path,path[1:]))
    return edges,source,groups


def pair_weights(groups):
    a,b,t = groups
    weights = {}
    for left,right,sign in ((a,b,-1),(b,t,1),(a,t,1)):
        for i,j in product(left,right):
            weights[tuple(sorted((i,j)))] = F(sign*len(b),2*len(left)*len(right))
    for i,j in combinations(t,2):
        weights[i,j] = -F(len(b),6)
    return weights


def cotree_weights(groups,exchange):
    # Phi(cotree) - Phi(A) from all B charges and their balance.
    a,b,t = groups
    weights = Counter()
    if not exchange:
        # In the two-loop control use Phi(L)-Phi(C) projected along A.
        for left,right,sign in (((-1,),t,1),(t,a,1),((-1,),a,-1)):
            for i,j in product(left,right):
                weights[tuple(sorted((i,j)))] += F(sign*len(a),2*len(left)*len(right))
        for i,j in combinations(t,2):
            weights[i,j] -= F(len(a),6)
        return weights
    pairs = ((a,b),(b,a))
    for target,load in pairs:
        for left,right,sign in (((-1,),t,1),(target,load,1),
                                ((-1,),load,-1),(target,t,-1)):
            weight = F(sign*len(load),2*len(left)*len(right)*len(pairs))
            for i,j in product(left,right):
                weights[tuple(sorted((i,j)))] += weight
    return weights


def averaged_voltage(bulk,groups,load_group=1):
    # The squared cofactor is fixed by permuting the retained basis and balance.
    a,b,t = groups
    charge = groups[load_group]
    # Reciprocity and the squared-cofactor symmetry allow each retained charge
    # assignment to be tested with this same path basis.
    original = bulk.source
    expected = [F(0),F(0)]
    for i in a:
        for j,_,ground in permutations(t):
            bulk.source = [(int(v in charge)-len(charge)*int(v==ground),)
                           for v in range(len(original))]
            for k,x in enumerate(bulk.voltage(i,j,(1,))):
                expected[k] += F(x,len(a)*6)
    bulk.source = original
    return tuple(expected)


def fixture(n):
    basis = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    p = [basis[0]]*n+[basis[1]]*2+[basis[2]]*(n-1)+basis[3:]
    return p+[tuple(-sum(x) for x in zip(*p))]


def collect(name,checkpoint=None,check=False,cotree=False):
    core = CORES[name]
    n = 2 if name=='theta' else 3
    saved = {}
    if checkpoint and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            idx,x,y = line.split()
            saved[int(idx)] = (F(x),F(y))
    total = [F(0),F(0)]
    tested = False
    longitudinal = F(0)
    for i,(config,factor) in enumerate(configurations(core)):
        if i in saved:
            values = saved[i]
        else:
            edges,source,groups = subdivision(core,config)
            if len({tuple(sorted(edge)) for edge in edges})<len(edges):
                values = (F(0),F(0))
            else:
                bulk = LabeledBulk(edges,source,all_orders=True)
                values = (bulk.resistance_pairing(cotree_weights(groups,n==3),cotree=True) if cotree
                          else bulk.resistance_pairing(pair_weights(groups)))
                if check and bulk.moments.base:
                    if not tested and not cotree:
                        assert values==averaged_voltage(bulk,groups)
                        tested = any(values)
                    if n==2:
                        longitudinal += factor*averaged_voltage(bulk,groups,0)[0]
                values = tuple(x*factor for x in values)
            if checkpoint:
                with checkpoint.open('a') as out:
                    print(i,*values,sep='\t',file=out)
        for k,v in enumerate(values):
            total[k] += v
    print(name,i+1,*total,flush=True)
    if name=='theta':
        from bulk_eight import potentials
        phi = potentials(fixture(n))
        control = phi[0][1]-phi[-3][1] if cotree else phi[n][2]-phi[-3][2]
        assert not cotree or control==-F(7,77760)
        assert total==[control,0]
        if check and not saved:
            assert longitudinal==phi[n][1]-phi[-3][1]==F(17,77760)
            print('Longitudinal theta normalization:',longitudinal)
    else:
        known = {
            ('tetrahedron',False):(-F(9261411221,56491532697600),F(459733,3405402000)),
            ('double-square',False):(-F(829523,1975228416000),F(0)),
            ('tetrahedron',True):(-F(85970041,38516954112000),F(23449,14859936000)),
            ('double-square',True):(-F(23503,219469824000),F(0))}
        assert tuple(total)==known[name,cotree]
    return tuple(total)


def boundary_check():
    import sys
    from tempfile import TemporaryDirectory
    from subprocess import run
    from flint import fmpz_mat
    from bulk_three_hexagon import potentials
    directory = Path(__file__).resolve().parent/'operator'
    sys.path.insert(0,str(directory))
    from three_loop_soft_jet import exact_boundary
    p = fixture(3)
    hexagon = potentials(p)
    assert (hexagon[3][2]-hexagon[7][2],hexagon[0][2]-hexagon[3][2])==(
        -F(1,1451520),-F(1,13063680))
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(directory/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        symbols = exact_boundary(exe,p)
    rows = []
    for coefficients in ((-4,-2,4,1),(-12,-16,12,5)):
        row = [0]*18
        for i,c in zip((2,7,8,17),coefficients):
            row[i] = c
        rows.append(row)
    rows += [[a[i][2]-a[j][2] for a in symbols] for i,j in ((3,7),(0,3))]
    assert fmpz_mat(rows).rank()==4
    print('Two planar and two nonplanar rows:',rows,flush=True)


if __name__=='__main__':
    cli = ArgumentParser()
    cli.add_argument('core',choices=[*CORES,'boundary'])
    cli.add_argument('--checkpoint',type=Path)
    cli.add_argument('--check',action='store_true')
    cli.add_argument('--cotree',action='store_true')
    args = cli.parse_args()
    if args.core=='boundary':
        boundary_check()
    else:
        collect(args.core,args.checkpoint,args.check,args.cotree)
