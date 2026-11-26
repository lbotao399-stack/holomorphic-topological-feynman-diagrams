"""Primary background, zero-momentum bilinear mark and the genus-two bulk bracket."""
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from bulk_eight import BULK
from marked_boundaries import quartics
from point_graph import chain_graph, solve_graph
from scalar_fan import cofactor


def bivalent_zero(n):
    edges,nv = chain_graph(n)
    source = [(0,)*5 for _ in range(nv)]
    # Vertex 0 is bivalent, with zero source; choose the charge at 1 to balance.
    for i in range(2,nv):
        if sum(i in e for e in edges)==2:
            source[i] = tuple((i+1)**j for j in range(5))
    source[1] = tuple(-sum(q) for q in zip(*source))
    times = [1+(3*i*i+5*i)%19 for i in range(len(edges))]
    from square_three import G
    graph = G(nv,edges)
    flow,_ = solve_graph(edges,source,times,1)
    at = [i for i,e in enumerate(edges) if 0 in e]
    assert len(at)==2
    for cycle in graph.C:
        assert [cycle[at[0]]*x for x in flow[at[0]]]==[cycle[at[1]]*x for x in flow[at[1]]]
    assert cofactor(flow,graph.C,times)==0


def bulk_tree():
    p = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    p.append((2,1,-1,-2,1))
    p.append(tuple(-sum(col) for col in zip(*p)))
    value = sum(a*b for a,b in zip(BULK,quartics(p)))
    assert value
    for i in range(7):
        tree = [sum(p[j][k]*value for j in range(7) if j!=i) for k in range(5)]
        assert tree==[-x*value for x in p[i]]
    return value


if __name__=='__main__':
    for n in range(1,6):
        bivalent_zero(n)
    print('Bivalent integrated mark through loop five: all current cofactors zero.')
    print('Seven-primary longitudinal tree coefficient:',bulk_tree())
