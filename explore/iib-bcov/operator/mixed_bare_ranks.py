"""Length-independent current-rank exclusions in the mixed four-loop Ward integral."""
from functools import cache
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from flint import fmpz_mat,nmod_mat

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from quartic_planar_cores import cubic_cores


PAIRS = tuple(combinations(range(9),2))
WAVES = ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(0,0,0,1),
         (-1,-1,-1,-2),(0,-1,1,0),(0,-1,1,0))


def supported_cycles(core):
    for k in range(2,10):
        for paths in combinations(range(9),k):
            graph = nx.MultiGraph(core[e] for e in paths)
            h = k-len(graph)+1
            need = 4*h+1-k
            if need<=0 or not nx.is_connected(graph) or nx.has_bridges(graph):
                continue
            incidence = fmpz_mat([[int(core[e][0]==v)-int(core[e][1]==v)
                                   for e in paths] for v in sorted(graph)])
            kernel,dimension = incidence.nullspace()
            assert dimension==h
            cols = [(0,)*h]*9
            for i,e in enumerate(paths):
                vec = tuple(int(kernel[i,j]) for j in range(h))
                if next(x for x in vec if x)<0:
                    vec = tuple(-x for x in vec)
                cols[e] = vec
            yield paths,cols,need


def rank_conditions(cols,need):
    """The first and last momentum coordinates separate from the coupled middle pair."""
    zero = (0,)*len(cols[0])
    @cache
    def rank(rows):
        return fmpz_mat(rows).rank()

    one = np.array([int(any(v)) for v in cols],dtype=np.int8)
    two = np.array([[rank((i,j)) for j in cols] for i in cols],dtype=np.int8)
    three = np.array([[[rank((i,j,k)) for k in cols] for j in cols] for i in cols],dtype=np.int8)
    masks = np.zeros((9,9,9,2,need+1),dtype=np.uint64)
    @cache
    def pair_ranks(b,c,balance,u,v):
        # Hard waves 1,2 carry the second and third coordinates separately.
        # Both soft columns carry their difference. The balance has equal components.
        mid = (b+zero,zero+c,tuple(-x for x in u)+u,tuple(-x for x in v)+v)
        base = rank(mid)
        return base,rank((*mid,balance+balance))
    for b,c,k in product(range(9),repeat=3):
        passes = [[0]*(need+1) for _ in range(2)]
        for bit,(i,j) in enumerate(PAIRS):
            base,extended = pair_ranks(cols[b],cols[c],cols[k],cols[i],cols[j])
            for flag,size in enumerate((extended,base+1)):
                for threshold in range(min(size,need)+1):
                    passes[flag][threshold] |= 1<<bit
        masks[b,c,k] = passes
    # Direct full-matrix checks of the middle-coordinate elimination on each subgraph.
    for n in range(40):
        ids = tuple((n//(2+j)+3*n+j*j)%9 for j in range(6))
        i,j = PAIRS[(11*n+7)%len(PAIRS)]
        rows = [[a*b for a in cols[e] for b in q] for e,q in zip((*ids,i,j),WAVES)]
        a,b,c,d,f,k = ids
        flag = int(two[a,k]>one[a] or three[d,f,k]>two[d,f])
        threshold = max(0,need-int(one[a]+two[d,f]))
        assert (fmpz_mat(rows).rank()>=need)==bool(int(masks[b,c,k,flag,threshold])>>(PAIRS.index((i,j)))&1)
    return one,two,three,masks


def census(core,assignments):
    a,b,c,d,f,k = assignments
    full = (1<<len(PAIRS))-1
    geometric = np.full(a.size,full,dtype=np.uint64)
    physical = geometric.copy()
    boundary = geometric.copy()
    mark_ranks,no_new_direction = [],[]
    for mask in range(64):
        rows = [list(p) for i,p in enumerate(WAVES[:6]) if mask>>i&1]
        r = fmpz_mat(rows).rank() if rows else 0
        mark_ranks.append(r)
        no_new_direction.append(fmpz_mat([*rows,list(WAVES[6])]).rank()==r)
    for paths,cols,need in supported_cycles(core):
        # First retain only the counts, with the two soft legs on distinct core paths.
        inside = np.array([int(e in paths) for e in range(9)],dtype=np.int8)
        counts = sum(inside[row] for row in assignments)
        geometric_masks = [sum(1<<n for n,(i,j) in enumerate(PAIRS)
                               if inside[i]+inside[j]>=t) for t in range(need+1)]
        geometric &= np.asarray(geometric_masks,dtype=np.uint64)[np.maximum(0,need-counts)]
        if not np.any(geometric):
            return geometric,geometric.copy(),geometric.copy()
        one,two,three,masks = rank_conditions(cols,need)
        flag = np.asarray((two[a,k]>one[a]) | (three[d,f,k]>two[d,f]),dtype=np.int8)
        threshold = np.maximum(0,need-one[a]-two[d,f])
        physical &= masks[b,c,k,flag,threshold]
        # Every current in H is a linear combination of its hard sources and its
        # boundary injections. Their sum is prescribed by the internal sources.
        vertices = {v for e in paths for v in core[e]}
        exposed = vertices & {v for e in range(9) if e not in paths for v in core[e]}
        local_sources = sum(inside[q].astype(np.uint8)<<i for i,q in enumerate(assignments))
        h = len(cols[0])
        defect = 4-np.asarray(mark_ranks)[local_sources]-max(0,len(exposed)-1)
        needs_soft = sum(1<<bit for bit,(i,j) in enumerate(PAIRS) if two[i,j]==h)
        allows_soft = (defect==1) & ~np.asarray(no_new_direction)[local_sources]
        boundary &= np.where(defect<=0,np.uint64(full),
                             np.where(allows_soft,np.uint64(needs_soft),np.uint64(0)))
    assert not np.any(physical & ~geometric)
    return geometric,physical,physical & boundary


def totals(mask,assignments):
    eligible = mask!=0
    placement = np.sum(7**assignments.astype(np.int32),axis=0)
    ordered = np.ones(len(mask),dtype=np.int32)
    fac = np.array([factorial(i) for i in range(7)])
    for e in range(9):
        ordered *= fac[np.sum(assignments==e,axis=0)]
    return (int(eligible.sum()),len(np.unique(placement[eligible])),int(ordered[eligible].sum()))


def chamber_orbits(core,physical):
    """Six ordered sites modulo oriented core symmetries and the three label exchanges."""
    from fourfold_quartic_orbits import oriented_symmetries
    actions = tuple(oriented_symmetries(core))
    colors = sorted(set(permutations((0,0,1,1,2,2))))
    visited = set()
    for paths in combinations_with_replacement(range(9),6):
        if paths in visited:
            continue
        shapes = {tuple(sorted(g[e][0] for e in paths)) for g in actions}
        visited.update(shapes)
        stabilizer = []
        for g in actions:
            order = sorted(range(6),key=lambda i:(g[paths[i]][0],-i if g[paths[i]][1] else i))
            if tuple(g[paths[i]][0] for i in order)==paths:
                stabilizer.append(order)
        assert len(shapes)*len(stabilizer)==len(actions)
        seen = set()
        for color in colors:
            if color in seen:
                continue
            orbit = {tuple(color[i] for i in order) for order in stabilizer}
            seen.update(orbit)
            u,t,z = [[i for i,c in enumerate(color) if c==j] for j in range(3)]
            slots = [u[0],*t,*z,u[1]]
            index = sum(paths[i]*9**(5-j) for j,i in enumerate(slots))
            if physical[index]:
                multiplicity = 8*len(orbit)*len(shapes)
                yield paths,color,multiplicity,Fraction(multiplicity,len(actions))


def check_current(core,paths):
    from fourfold_mixed_soft_ward import V,integrated_soft_placements
    from point_graph import solve_graph
    from square_three import G
    core_cycles = G(6,core).C
    edges,parents = [],[]
    for e,(a,b) in enumerate(core):
        chain = [a,*[6+i for i,p in enumerate(paths) if p==e],b]
        edges.extend(zip(chain,chain[1:]))
        parents.extend([e]*(len(chain)-1))
    cycles = [[c[e] for e in parents] for c in core_cycles]
    lengths = [2+(i*i+3*i)%17 for i in range(15)]
    source = [(0,)*4]*6+list(WAVES[:6])
    currents,_ = solve_graph(edges,source,lengths,0)
    square = integrated_soft_placements(currents,cycles,lengths,V)
    if not square:
        return square
    # Grounded Green matrix on the hard sites, with compensating unit sinks at the root.
    impulses = [(-1,)*6,*[(0,)*6]*5,*[tuple(int(i==j) for i in range(6)) for j in range(6)]]
    _,green = solve_graph(edges,impulses,lengths,0)
    def resistance(i,j):
        return green[6+i][i]+green[6+j][j]-2*green[6+i][j]
    u,t,z = (0,5),(1,2),(3,4)
    def cross(left,right):
        return sum(resistance(i,j) for i in left for j in right)
    expected = resistance(*u)+(cross(t,z)-cross(u,z)-cross(t,u))/4
    total = 0
    for (a,b),(c,d),(f,g) in product(permutations(u),permutations(t),permutations(z)):
        order = [a,c,d,f,g,b]
        waves = source[:]
        for i,j in enumerate(order):
            waves[6+j] = WAVES[i]
        flow,phi = solve_graph(edges,waves,lengths,0)
        assert integrated_soft_placements(flow,cycles,lengths,V)==square
        total += phi[6+a][3]-phi[6+c][3]
    assert total==8*expected
    return square


def check_census(core,geometric,physical,assignments):
    rejected = np.flatnonzero((geometric!=0) & (physical==0))
    for i in rejected[::max(1,len(rejected)//3)][:3]:
        assert not check_current(core,assignments[:,i])
    retained = np.flatnonzero(physical)
    if len(retained):
        sample = retained[::max(1,len(retained)//32)][:32]
        assert any(check_current(core,assignments[:,i]) for i in sample)


def modular_network(core,paths,core_cycles,lengths,prime):
    edges,parents = [],[]
    for e,(a,b) in enumerate(core):
        chain = [a,*[6+i for i,p in enumerate(paths) if p==e],b]
        edges.extend(zip(chain,chain[1:]))
        parents.extend([e]*(len(chain)-1))
    incidence = np.array([[int(a==v)-int(b==v) for v in range(1,12)] for a,b in edges])
    conductance = np.array([pow(int(length),-1,prime) for length in lengths])
    weighted = conductance[:,None]*incidence
    laplace = nmod_mat((incidence.T@weighted%prime).tolist(),prime)
    if not laplace.det():
        return None
    green = np.array(laplace.inv().tolist(),dtype=np.int64)[:,5:]
    transfer = weighted@green%prime
    green = green[5:]
    resistance = (green.diagonal()[:,None]+green.diagonal()[None,:]-2*green)%prime
    return transfer,resistance,np.array(core_cycles)[:,parents]


def separating_pairs(core):
    for i,j in PAIRS:
        graph = nx.MultiGraph(core[e] for e in range(9) if e not in (i,j))
        parts = list(nx.connected_components(graph))
        if len(parts)==2:
            sides = [tuple(e for e in range(9) if e not in (i,j) and core[e][0] in part)
                     for part in parts]
            yield (i,j),*sides


def extra_length_mode(paths,waves,cuts):
    for connectors,left,right in cuts:
        def sourced(edges):
            return [waves[i].tolist() for i,e in enumerate(paths) if e in edges]
        marks,a,b = sourced(connectors),sourced(left),sourced(right)
        net = [sum(q[i] for q in a) for i in range(4)]
        fixed = [*marks,net]
        rank = fmpz_mat(fixed).rank()
        additional = all(fmpz_mat([*fixed,*side]).rank()>rank for side in (a,b))
        if min(4,1+rank+int(additional))<=2+len(marks):
            return True
    return False


def hard_waves(color):
    u,t,z = [[i for i,c in enumerate(color) if c==j] for j in range(3)]
    waves = np.empty((6,4),dtype=np.int64)
    waves[[u[0],*t,*z,u[1]]] = WAVES[:6]
    return (u,t,z),waves


def integral_orbits(core,physical):
    cuts = tuple(separating_pairs(core))
    for item in chamber_orbits(core,physical):
        paths,color,_,_ = item
        if not cuts or not extra_length_mode(paths,hard_waves(color)[1],cuts):
            yield item


def witness_counts(core,orbits,unresolved=None):
    """Nonzero residues certify nonzero rational metric integrands over Q."""
    from random import Random
    from square_three import G
    prime = 1000003
    core_cycles = G(6,core).C
    lengths = [Random(seed).sample(range(2,257),15) for seed in (13,19,29)]
    count = [0]*len(lengths)
    previous = None
    for item in orbits:
        paths,color,_,_ = item
        if paths!=previous:
            networks = [modular_network(core,paths,core_cycles,t,prime) for t in lengths]
            previous = paths
        (u,t,z),waves = hard_waves(color)
        for fixture,network in enumerate(networks):
            if network is None:
                continue
            transfer,resistance,cycles = network
            voltage = 4*resistance[u[0],u[1]]
            for sign,left,right in ((1,t,z),(-1,u,z),(-1,t,u)):
                voltage += sign*sum(resistance[i,j] for i in left for j in right)
            if voltage%prime==0:
                continue
            currents = transfer@waves%prime
            rows = currents.T[:,None,:]*cycles[None,:,:]
            one = np.ones((1,15),dtype=np.int64)
            projected = np.vstack((one,rows[0],rows[1]+rows[2],rows[3]))
            if nmod_mat((projected%prime).tolist(),prime).rank()<13:
                continue
            full = np.vstack((one,rows.reshape(16,15)))
            if nmod_mat((full%prime).tolist(),prime).rank()==15:
                count[fixture] += 1
                break
        else:
            if unresolved is not None:
                unresolved.append(item)
    return tuple(count)


def check_bipartite_two_stars():
    """Count hard placements permitting two successive unmarked core Y-to-delta moves."""
    core = next(c for c in cubic_cores() if len(set(c))==9 and nx.is_bipartite(nx.Graph(c)))
    graph = nx.Graph(core)
    stars = {v:{i for i,e in enumerate(core) if v in e} for v in graph}
    cuts = {frozenset(paths):need for paths,cols,need in supported_cycles(core)}
    obstructions = set()
    for v in graph:
        reduced = graph.copy()
        neighbors = list(reduced[v])
        reduced.remove_node(v)
        reduced.add_edges_from(combinations(neighbors,2))
        for w in reduced:
            if reduced.degree(w)==3:
                assert w not in neighbors
                union = frozenset(stars[v] | stars[w])
                assert len(union)==6 and cuts[union]==3
                obstructions.add(union)
    assert len(obstructions)==6
    placements = [p for p in combinations_with_replacement(range(9),6)
                  if any(cut.isdisjoint(p) for cut in obstructions)]
    types = {tuple(sorted((p.count(e) for e in set(p)),reverse=True)):p for p in placements}
    assert len(placements)==159 and len(types)==7
    assert all(check_current(core,p)==0 for p in types.values())
    return len(placements),sorted(types)


if __name__=='__main__':
    print('Bipartite two-star exclusions:',check_bipartite_two_stars(),flush=True)
    assignments = np.indices((9,)*6,dtype=np.int8).reshape(6,-1)
    expected = ((881,821,797),(0,0,0),(4467,4173,3995),
                (11732,11732,11732),(2575,2575,2575))
    for core,reference in zip(cubic_cores(),expected,strict=True):
        geometric,jumps,physical = census(core,assignments)
        first,second = totals(geometric,assignments),totals(physical,assignments)
        orbits = list(chamber_orbits(core,physical)) if second[0] else []
        assert sum(q[2] for q in orbits)==second[2]
        check_census(core,geometric,physical,assignments)
        survivors = list(integral_orbits(core,physical)) if second[0] else []
        witnesses = witness_counts(core,survivors)
        assert (len(orbits),len(survivors),sum(witnesses))==reference
        print(core,first,totals(jumps,assignments),second,'orbits',len(orbits),
              'after cuts',len(survivors),'witnesses',witnesses,flush=True)
