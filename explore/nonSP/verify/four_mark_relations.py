"""Four-loop holomorphic quadratic relations in the signed standard S4 source sector."""
from collections import Counter
from itertools import combinations
from math import gcd
from functools import reduce
import networkx as nx
from flint import fmpz_mat


def edge_order(graph):
    return sorted(tuple(sorted(e)) for e in graph.edges())


def sign(sequence):
    return (-1)**sum(x>y for i,x in enumerate(sequence) for y in sequence[i+1:])


def edge_sign(edges, image):
    return sign([tuple(sorted((image[i],image[j]))) for i,j in edges])


def pairings(vertices):
    for i,j,k,l in combinations(sorted(vertices),4):
        yield ((i,j),(k,l))
        yield ((i,k),(j,l))
        yield ((i,l),(j,k))


def map_pair(pair, image):
    (i,j),(k,l) = pair
    a,b,c,d = [image[v] for v in (i,j,k,l)]
    if len({a,b,c,d})!=4:
        return None,0
    left,right = tuple(sorted((a,b))),tuple(sorted((c,d)))
    orient = (-1)**((a>b)+(c>d)+(left>right))
    return tuple(sorted((left,right))),orient


def admissible(graph):
    edges = edge_order(graph)
    for k in range(3,len(graph)):
        for group in combinations(graph,k):
            if sum(i in group and j in group for i,j in edges)>2*k-3:
                return False
    return True


def types(v):
    return [g for g in nx.graph_atlas_g() if len(g)==v and g.number_of_edges()==(9 if v==6 else 10)
            and min(dict(g.degree()).values())>=2 and nx.is_connected(g) and admissible(g)]


def variables(graphs):
    assigned, reps = {}, []
    for number,graph in enumerate(graphs):
        edges = edge_order(graph)
        aut = [(p,edge_sign(edges,p)) for p in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter()]
        for pair in pairings(graph):
            if (number,pair) in assigned:
                continue
            orbit = set()
            for p,sgn in aut:
                other,weight = map_pair(pair,p)
                orbit.add((other,weight*sgn))
            if len({other for other,_ in orbit})<len(orbit):
                for other,_ in orbit:
                    assigned[number,other] = (0,0)
            else:
                for other,weight in orbit:
                    assigned[number,other] = (len(reps),weight)
                reps.append((number,pair))
    return reps,assigned


def contractions(graph, targets):
    result = []
    edges = edge_order(graph)
    for pos,(i,j) in enumerate(edges):
        collapsed = {v:(i if v==j else v) for v in graph}
        old = edges[:pos]+edges[pos+1:]
        mapped = [tuple(sorted((collapsed[a],collapsed[b]))) for a,b in old]
        if len(set(mapped))!=9:
            continue
        quotient = nx.Graph(mapped)
        for number,other in enumerate(targets):
            if sorted(dict(other.degree()).values())!=sorted(dict(quotient.degree()).values()):
                continue
            matcher = nx.algorithms.isomorphism.GraphMatcher(quotient,other)
            if matcher.is_isomorphic():
                image = {v:matcher.mapping[collapsed[v]] for v in graph}
                result.append((number,(-1)**pos*edge_sign(old,image),image))
                break
    return result


def matrix(pure_only=False,metadata=False):
    graphs,sliding = types(6),types(7)
    assert len(graphs)==13 and len(sliding)==60
    reps,assigned = variables(graphs)
    rows = {}
    for index,graph in enumerate(sliding):
        lower = [subset for k in (3,4,5) for subset in combinations(graph,k)
                 if graph.subgraph(subset).number_of_edges()==2*k-3]
        if pure_only and lower:
            continue
        terms = contractions(graph,graphs)
        for pair in pairings(graph):
            row = Counter()
            for number,weight,image in terms:
                other,mu = map_pair(pair,image)
                if mu:
                    k,w = assigned[number,other]
                    row[k] += weight*mu*w
            row = {k:v for k,v in row.items() if v}
            if row:
                norm = reduce(gcd,map(abs,row.values()))
                if row[min(row)]<0:
                    norm=-norm
                key = tuple(sorted((k,v//norm) for k,v in row.items()))
                witness = (len(lower),index,pair,norm)
                if key not in rows or witness<rows[key]:
                    rows[key]=witness
    records = [rows[key] for key in sorted(rows)]
    rows = sorted(rows)
    mat = fmpz_mat(len(rows),len(reps))
    for i,row in enumerate(rows):
        for j,k in row:
            mat[i,j]=k
    return (graphs,reps,mat,records) if metadata else (graphs,reps,mat)


def main():
    graphs,reps,mat=matrix()
    name = lambda g:nx.to_graph6_bytes(g,header=False).decode().strip()
    counts = Counter(i for i,_ in reps)
    print('variables per graph:', [(name(g),counts[i]) for i,g in enumerate(graphs)],flush=True)
    n=mat.ncols()
    basis,nullity = mat.nullspace()
    assert mat*basis == fmpz_mat(mat.nrows(),n)
    assert (mat.nrows(),n,n-nullity) == (1050,200,197)
    print('matrix',mat.nrows(),n,'rank',n-nullity,'nullity',nullity,flush=True)
    for i,(number,pair) in enumerate(reps):
        if name(graphs[number])=='ElUg':
            print('K33 pairing',pair,'row in nullspace',list(basis[i,j] for j in range(nullity)),flush=True)
    # Nullspace basis is implementation-dependent. Test the rank restriction
    # and the source-pair identities instead of a prescribed basis.
    restriction = fmpz_mat([[basis[i,j] for j in range(nullity)] for i in (195,196,197,198,199)])
    assert restriction.rank()==3
    assert fmpz_mat([[1,0,-1,0,0],[2,4,0,1,-1]])*restriction == fmpz_mat(2,nullity)
    assert fmpz_mat([[basis[i,j] for j in range(nullity)] for i in (195,196,199)]).rank()==3
    return graphs,reps,mat,basis,nullity


if __name__=='__main__':
    main()
