"""Symbolic two-terminal reductions of the last mixed four-loop metric ranks."""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import sys

import numpy as np
from flint import fmpq_mat,fmpz_mpoly_ctx
from mixed_bare_ranks import census,cubic_cores,hard_waves,integral_orbits,witness_counts
from square_three import G


def symbolic_rank(matrix):
    """Bareiss elimination with row and column pivoting over integer polynomials."""
    previous = None
    rank = 0
    while matrix and matrix[0]:
        choices = [(len(p),i,j) for i,row in enumerate(matrix) for j,p in enumerate(row) if p]
        if not choices:
            break
        _,i,j = min(choices)
        matrix[0],matrix[i] = matrix[i],matrix[0]
        for row in matrix:
            row[0],row[j] = row[j],row[0]
        head = matrix[0]
        pivot = head[0]
        next_matrix = []
        for row in matrix[1:]:
            reduced = [pivot*x-row[0]*y for x,y in zip(row[1:],head[1:])]
            if previous is not None and previous!=1:
                reduced = [p/previous for p in reduced]
            next_matrix.append(reduced)
        matrix = next_matrix
        previous = pivot
        rank += 1
    return rank


def macro(core_index,grouping):
    if grouping=='right':
        return [(0,1)]*3,[(0,),(1,),(2,3,4,5,6,7,8)],[(),(),(4,5,6,7,8)]
    if grouping=='left' and core_index==0:
        return [(0,1)]*3,[(7,),(8,),(0,1,2,3,4,5,6)],[(),(),(0,1,2,3,4)]
    if core_index==0:
        return [(0,1)]*3,[(0,1,2,3),(4,),(5,6,7,8)],[(0,1),(),(7,8)]
    if core_index==2:
        return ([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
                [(0,1,2,3),(4,),(5,),(6,),(7,),(8,)],[(0,1),(),(),(),(),()])
    raise ValueError('Expected one of the two reducible cores having undecided chambers')


def symbolic_certificate(core_index,paths,color,grouping='short'):
    edges,parts,blocks = macro(core_index,grouping)
    nodes = max(v for edge in edges for v in edge)+1
    cycle = G(nodes,edges).C
    pairs = tuple(combinations_with_replacement(range(1,nodes),2))
    names = tuple('s'+str(i) for i in range(6))+tuple('g'+str(i)+str(j) for i,j in pairs)
    ctx = fmpz_mpoly_ctx.get(names)
    split = ctx.gens()[:6]
    green = dict(zip(pairs,ctx.gens()[6:]))
    if nodes==2:
        green[(1,1)] = ctx.constant(1)
    zero = ctx.constant(0)
    source = [[zero]*4 for _ in range(nodes)]
    columns = []
    _,waves = hard_waves(color)
    for e,((a,b),part,block) in enumerate(zip(edges,parts,blocks)):
        inside = [i for i,path in enumerate(paths) if path in part]
        vectors = [[int(x) for x in waves[i]] for i in inside if paths[i] not in block]
        total = [sum(int(waves[i,mu]) for i in inside) for mu in range(4)]
        left = [sum((split[i]*int(waves[i,mu]) for i in inside),zero) for mu in range(4)]
        if block:
            vectors.append([sum(int(waves[i,mu]) for i in inside if paths[i] in block)
                            for mu in range(4)])
            vectors.append(left)
        for mu in range(4):
            source[a][mu] += left[mu]
            source[b][mu] += total[mu]-left[mu]
        for q in vectors:
            columns.append([zero+c[e]*x for c in cycle for x in q])
    potential = [[zero]*4]
    for i in range(1,nodes):
        potential.append([sum((green[tuple(sorted((i,j)))]*source[j][mu]
                               for j in range(1,nodes)),zero) for mu in range(4)])
    for e,(a,b) in enumerate(edges):
        columns.append([c[e]*(x-y) for c in cycle for x,y in zip(potential[a],potential[b])])
    # One independent scale for each block and one independent length for each other segment.
    parameters = 15-sum(len(block)-1+sum(path in block for path in paths) for block in blocks if block)
    assert len(columns)==parameters
    return parameters,symbolic_rank(columns)


def direct_control(core_index,paths,color,grouping,expected_rank):
    """Recover the six hard-source splits, terminal potentials, and variation rank directly."""
    from point_graph import solve_graph
    from quartic_seven_plane_lift import fixture
    core = cubic_cores()[core_index]
    edges,parts,blocks = macro(core_index,grouping)
    if grouping=='right':
        terminals = (0,1)
    elif grouping=='left':
        terminals = (4,5)
    else:
        terminals = (2,3) if core_index==0 else (2,3,4,5)
    full = fixture(core,paths)
    parents = [e for e in range(9) for _ in range(paths.count(e)+1)]
    lengths = [2+(i*i+3*i)%17 for i in range(15)]
    _,waves = hard_waves(color)
    sources = [(0,)*4]*6+waves.tolist()
    flow,phi = solve_graph(full,sources,lengths,terminals[0])
    deposits = [[Fraction() for _ in range(4)] for _ in terminals]
    resistances = []
    for (a,b),part in zip(edges,parts):
        ia,ib = terminals[a],terminals[b]
        ids = [i for i,e in enumerate(parents) if e in part]
        vs = sorted({v for i in ids for v in full[i]})
        subnet = [(vs.index(full[i][0]),vs.index(full[i][1])) for i in ids]
        dipole = [[int(v==ia)-int(v==ib)] for v in vs]
        _,unit = solve_graph(subnet,dipole,[lengths[i] for i in ids],vs.index(ib))
        resistance = unit[vs.index(ia)][0]
        resistances.append(resistance)
        left = [Fraction() for _ in range(4)]
        for i,p in enumerate(paths):
            if p not in part: continue
            weight = unit[vs.index(6+i)][0]/resistance
            for mu,x in enumerate(waves[i]):
                left[mu] += weight*int(x)
                deposits[b][mu] += (1-weight)*int(x)
        for mu in range(4):
            deposits[a][mu] += left[mu]
            initial = sum((int(full[i][0]==ia)-int(full[i][1]==ia))*flow[i][mu] for i in ids)
            assert initial==(phi[ia][mu]-phi[ib][mu])/resistance-left[mu]
    _,reduced_phi = solve_graph(edges,[[str(x) for x in row] for row in deposits],resistances,0)
    assert all(reduced_phi[i]==phi[v] for i,v in enumerate(terminals))
    cycle = G(6,core).C
    cols = [[lengths[i]*c[e]*x for c in cycle for x in flow[i]] for i,e in enumerate(parents)]
    blocked = {e for block in blocks for e in block}
    scales = [[sum(col[j] for col,e in zip(cols,parents) if e in block) for j in range(16)]
              for block in blocks if block]
    scaled = scales+[col for col,e in zip(cols,parents) if e not in blocked]
    assert fmpq_mat([[str(x) for x in col] for col in scaled]).rank()==expected_rank


def classified_orbits(core_index):
    """Return every nonzero chamber with its labeled multiplicity and symmetry factor."""
    core = cubic_cores()[core_index]
    assignments = np.indices((9,)*6,dtype=np.int8).reshape(6,-1)
    _,_,physical = census(core,assignments)
    if not np.any(physical):
        return [],Counter()
    candidates = list(integral_orbits(core,physical))
    missing = []
    witness_counts(core,candidates,missing)
    excluded = set()
    ranks = Counter()
    for item in missing:
        paths,color,_,_ = item
        for choice in (('right','left','short') if core_index==0 else ('right','short')):
            n,k = symbolic_certificate(core_index,paths,color,choice)
            if n-k>=2:
                if (choice,n,k) not in ranks:
                    direct_control(core_index,paths,color,choice,k)
                ranks[choice,n,k] += 1
                excluded.add(item)
                break
        else:
            raise ArithmeticError(f'No witness or all-length exclusion for {paths}, {color}')
    return [item for item in candidates if item not in excluded],ranks


if __name__=='__main__':
    for index in ((int(sys.argv[1]),) if len(sys.argv)>1 else range(5)):
        records,ranks = classified_orbits(index)
        assert len(records)==(797,0,3995,11732,2575)[index]
        assert sum(ranks.values())==(24,0,178,0,0)[index]
        print(index,'nonzero chambers',len(records),'polynomial exclusions',sorted(ranks.items()),flush=True)
