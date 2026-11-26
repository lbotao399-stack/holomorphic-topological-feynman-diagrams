"""Positive current-insertion witness on cubic product forests at root sources."""
from fractions import Fraction as F
from forest_orders import forests, graph
from forest_jets import edge
from vk import G, currents_exact, density_exact, fdet
import networkx as nx


def check(n):
    count=0
    for forest in forests(n):
        edges,a,b,roots,_=graph(forest)
        network=G(n+2,edges)
        source=[(0,0)]*n+[(1,0),(-1,-len(roots))]
        for v in roots:
            source[v]=(0,1)
        mark=edges.index((a,b))
        # Deterministic positive rational choices with unequal paired fractions.
        for offset in range(1,5):
            times=[F(1+(i+offset)**3) for i in range(len(edges))]
            y=currents_exact(network,times,source[:-1])
            m=[[sum(network.C[i][e]*network.C[j][e]*times[e] for e in range(len(edges)))
                for j in range(n)] for i in range(n)]
            minor=(-1)**mark*fdet([[cyc[e]*y[e][j] for e in range(len(edges)) if e!=mark]
                                    for j in range(2) for cyc in network.C])
            d=density_exact(network,times,source[:-1])
            assert minor==times[mark]*fdet(m)**2*d
            if minor:
                assert minor*d>0
                break
        else:
            raise ArithmeticError(f'No nonzero sample for {forest}.')
        internal=nx.Graph()
        internal.add_nodes_from(range(n))
        internal.add_edges_from((i,j) for i,j in edges if i<n and j<n)
        positions={p:i for i,p in enumerate(edges)}
        for i,j in list(internal.edges()):
            internal.remove_edge(i,j)
            child=nx.node_connected_component(internal,i)
            if set(roots)&child:
                child=nx.node_connected_component(internal,j)
            assert not set(roots)&child
            for k in (0,1):
                flow=sum(y[positions[v,t]][k] for v in child for t in (a,b) if (v,t) in positions)
                assert flow==(1 if min(i,j) not in child else -1)*y[positions[edge(i,j)]][k]
            internal.add_edge(i,j)
        count+=1
    print(n,'loops, positive cofactor at a rational point on',count,'forests',flush=True)


if __name__=='__main__':
    for n in range(1,6):
        check(n)
