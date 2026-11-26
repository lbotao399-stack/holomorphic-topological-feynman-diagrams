"""Exact span of graph-local source polynomials on the four remaining graph types."""
import networkx as nx
import numpy as np
from flint import fmpz_mat, nmod_mat

from contractions import edge_order, edge_sign, record, substitution
from homogeneous import signed_basis
from lower_faces import G, tree_flows
from sampling import sampling

PRIME = 1009


def determinants(matrices):
    """Batched finite-field determinants by elimination with row pivoting."""
    m=matrices.copy()%PRIME
    rows=np.arange(len(m))
    out=np.ones(len(m),dtype=np.int64)
    inverse=np.array([0]+[pow(x,-1,PRIME) for x in range(1,PRIME)],dtype=np.int64)
    for j in range(m.shape[1]):
        pivots=np.argmax(m[:,j:,j]!=0,axis=1)+j
        swap=rows[(pivots!=j)&(m[rows,pivots,j]!=0)]
        k=pivots[swap]
        original=m[swap,k,:].copy()
        m[swap,k,:]=m[swap,j,:]
        m[swap,j,:]=original
        out[swap]*=-1
        p=m[:,j,j]
        out=out*p%PRIME
        if j+1<m.shape[1]:
            ratio=m[:,j+1:,j]*inverse[p,None]%PRIME
            m[:,j+1:,j:]=(m[:,j+1:,j:]-ratio[:,:,None]*m[:,j,None,j:])%PRIME
    return out


def fiber_dimension(graph):
    dense=signed_basis(graph)
    bivalent=[v for v,d in graph.degree if d==2]
    if not bivalent:
        return dense.shape[1]
    w,=bivalent
    sources=np.concatenate((np.eye(5,dtype=np.int64),-np.ones((1,5),dtype=np.int64)))
    difference=256*dense
    for u in graph.neighbors(w):
        mat=2*sources
        mat[w]=sources[w]
        mat[u]+=sources[w]
        difference-=substitution(tuple(tuple(int(x) for x in row) for row in mat))@dense
    return dense.shape[1]-fmpz_mat(difference.toarray().tolist()).rank()


def check():
    # The seven-node rule integrates every polynomial in alpha of degree <= 6.
    nodes=[j*pow(8,-1,PRIME)%PRIME for j in range(1,8)]
    vandermonde=nmod_mat([[pow(t,d,PRIME) for t in nodes] for d in range(7)],PRIME)
    moments=nmod_mat([[pow(d+1,-1,PRIME)] for d in range(7)],PRIME)
    weights=vandermonde.inv()*moments
    rng=np.random.default_rng(3461)
    probes=np.array(sampling()[0][:88],dtype=np.int64)
    _,graphs=record()
    for i in (8,10,11,12):
        original=graphs[i]
        edges=edge_order(original)
        g=G(6,edges)
        identity=np.concatenate((np.eye(5,dtype=np.int64),-np.ones((1,5),dtype=np.int64)))
        trees=tree_flows(g,identity)
        cotrees=np.array(list(trees))
        tree_sources=np.array(list(trees.values())).reshape(len(trees),-1)
        automorphisms=[(tuple(p[k] for k in range(6)),edge_sign(edges,p))
                       for p in nx.isomorphism.GraphMatcher(original,original).isomorphisms_iter()]
        sources=np.concatenate([probes[:,permutation,:] for permutation,_ in automorphisms])
        signs=np.array([s for _,s in automorphisms])
        cycles=np.array(g.C)
        def row(times):
            # At these moduli the pair products below fit in signed 64-bit integers.
            time_products=(times[cotrees[:,:2]].prod(axis=1)%PRIME)*(times[cotrees[:,2:]].prod(axis=1)%PRIME)%PRIME
            y=((time_products@tree_sources).reshape(9,5)%PRIME@sources[:,:5,:]).transpose(0,2,1)
            matrices=(y[:,:,None,1:]*cycles[None,None,:,1:]).reshape(-1,8,8)
            cofactors=determinants(matrices).reshape(len(automorphisms),-1)
            return (signs@cofactors)*pow(int(times[0]),-1,PRIME)%PRIME
        bivalent=[v for v,d in original.degree if d==2]
        if bivalent:
            a,b=[e for e,ends in enumerate(edges) if bivalent[0] in ends]
        bound=fiber_dimension(original)
        table=[]
        for _ in range(bound+12):
            times=rng.integers(1,PRIME,size=9)
            if bivalent:
                total=times[a]
                values=np.zeros(len(probes),dtype=np.int64)
                for j,alpha in enumerate(nodes):
                    times[a]=total*alpha%PRIME
                    times[b]=(total-times[a])%PRIME
                    values=(values+int(weights[j,0])*row(times))%PRIME
            else:
                values=row(times)
            table.append(values)
        result=nmod_mat(np.array(table).tolist(),PRIME).rank()
        assert result==bound
        print(i,'Kirchhoff source span',result,flush=True)


if __name__=='__main__':
    check()
