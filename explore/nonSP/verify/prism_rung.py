"""The prism source pairing (05)|(12) from K3,3, 16 Ward rows, and three-loop moments."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
import json
from math import lcm
from pathlib import Path
import sys

import networkx as nx

from four_mark_relations import edge_order, matrix, sign, types
from three_loop_current import EDGES, integrate_numerator, numerator

# Indices refer to the sorted primitive rows produced by four_mark_relations.matrix.
CERT = {920:F(1,2),921:-F(1,2),930:-F(1,2),931:-1,932:1,933:-1,936:1,940:-1,
        944:-F(1,2),984:-2,985:1,995:1,998:-1,1003:1,1009:-1,1049:2}
MOMENTA = ((1,0),(0,1),(2,0),(-3,-1))
ZERO = (0,0)


def bracket(p,q):
    return p[0]*q[1]-p[1]*q[0]


def projection(pair):
    # Character projection of the pair stabilizer. Inducing this character
    # gives the signed standard S4 representation.
    p,q,r,s = MOMENTA
    amplitude = bracket(p,q)*bracket(r,s)*(bracket(r,s)**2-bracket(p,q)**2)
    for flip_i,flip_j,swap in product((0,1),repeat=3):
        left = pair[0][::-1] if flip_i else pair[0]
        right = pair[1][::-1] if flip_j else pair[1]
        sequence = right+left if swap else left+right
        yield dict(zip(sequence,MOMENTA)),F((-1)**(flip_i+flip_j+swap),8*amplitude)


def known_prism(pair):
    path = Path(__file__).resolve().parents[1]/'results_nonSP_C2.json'
    reference = json.loads(path.read_text())['prism_4src']['table']
    answer = Counter()
    for source,weight in projection(pair):
        a,b,c = [bracket(source[i],source[j]) for i,j in ((0,1),(0,2),(1,2))]
        for name,values in reference.items():
            i,j,k = map(int,name[-3:])
            for n,value in zip((0,3),values):
                answer[n] += weight*F(value)*a**i*b**j*c**k
    return answer


def faces(graph):
    edges = edge_order(graph)
    indexed = {tuple(sorted(e)):i for i,e in enumerate(EDGES)}
    for k in (3,4,5):
        for cut in combinations(graph,k):
            inner = [e for e in edges if set(e).issubset(cut)]
            if len(inner)!=2*k-3:
                continue
            assert k==3
            rest = [e for e in edges if e not in inner]
            shrink = {v:(cut[0] if v in cut else v) for v in graph}
            quotient = nx.Graph((shrink[i],shrink[j]) for i,j in rest)
            matcher = nx.algorithms.isomorphism.GraphMatcher(quotient,nx.Graph(EDGES))
            assert matcher.is_isomorphic()
            mapping = {v:matcher.mapping[shrink[v]] for v in graph}
            oriented = [(mapping[a],mapping[b]) for a,b in rest]
            ids = [indexed[tuple(sorted(e))] for e in oriented]
            # Split the ordered sliding edges, then put the quotient in the
            # precise edge order of the integrated seven-edge graph.
            face_sign = sign([edges.index(e) for e in inner+rest])*sign(ids)
            incidence = [Counter(),Counter()]
            for e,t,idx in zip(rest,oriented,ids):
                flip = 1 if t==EDGES[idx] else -1
                for pos,v in enumerate(cut[:2]):
                    incidence[pos][idx] += flip*(int(e[1]==v)-int(e[0]==v))
            yield cut,mapping,incidence,face_sign


def lower_sum(center,certificate,target_entries):
    _,reps,M,records=matrix(metadata=True)
    target = Counter()
    for row,weight in certificate.items():
        for j in range(M.ncols()):
            target[j] += weight*int(M[row,j])
    assert {i:v for i,v in target.items() if v} == target_entries
    assert [reps[i] for i in (195,196,198,199)] == [
        (11,((0,1),(2,3))),(11,((0,2),(1,3))),(11,((0,5),(1,2))),(12,((0,2),(1,3)))]
    graphs = types(7)
    total = 0
    jobs = []
    rational_scale = 1
    for row,weight in certificate.items():
        count,idx,pair,norm = records[row]
        graph = graphs[idx]
        found = list(faces(graph))
        assert len(found)==count
        for cut,mapping,incidence,face_sign in found:
            for source,projector in projection(pair):
                outer = [[0,0] for _ in range(5)]
                for v,P in source.items():
                    for i in range(2):
                        outer[mapping[v]][i] += P[i]
                # In project K conventions the (1,3) term has relative
                # sign (-1)^(1*3) against a (0,4) contraction. Moving the
                # term to the RHS leaves +K_triangle, namely -[tau_i,tau_j]/2.
                prefactor = -F(1,2)*face_sign*weight*projector/norm
                rational_scale = lcm(rational_scale,prefactor.denominator)
                jobs.append((tuple(map(tuple,outer)),source.get(cut[0],ZERO),source.get(cut[1],ZERO),
                             incidence,prefactor))
    @lru_cache(maxsize=4)
    def graph_numerator(source):
        return numerator(center,source)
    for outer,first,second,incidence,prefactor in sorted(jobs,key=lambda job:job[0]):
        q,y,det = graph_numerator(outer)
        if not det:
            continue
        shifted = []
        for p,terms in zip((first,second),incidence):
            shifted.append(tuple(p[j]*q+sum(c*y[e][j] for e,c in terms.items() if c) for j in (0,1)))
        total += int(prefactor*rational_scale)*det*bracket(*shifted)
    # Both star charts have relative orientation -1; all jobs have the
    # common denominator Q^9 and the rational scale was cleared above.
    result,margin = integrate_numerator(total,9,-F(1,rational_scale))
    print('center',center,'lower faces at projected sources',len(jobs),'polynomial terms',len(total),
          'margin',margin,'Ward RHS:',result,flush=True)
    return result


def run(center):
    lower = lower_sum(center,CERT,{195:2,196:4,198:1,199:-1})
    assert lower == {0:-F(9,6400),3:F(1,1200)}
    old = [known_prism(i) for i in (((0,1),(2,3)),((0,2),(1,3)))]
    k33 = {0:-F(539,14400),3:-F(1414,14400),5:F(2160,14400)}
    rung = {m:k33.get(m,0)+lower.get(m,0)-2*old[0][m]-4*old[1][m] for m in (0,2,3,4,5)}
    print('Known prism channels:',list(map(dict,old)),'K33 -> prism rung:',
          {m:v for m,v in rung.items() if v},flush=True)
    assert {m:v for m,v in rung.items() if v} == {0:-F(237,6400),3:-F(71,720),5:F(3,20)}
    difference = lower_sum(center,{1049:-1},{195:-1,197:1})
    assert difference == {0:F(1,2400),3:-F(1,4800)}
    last = {m:old[0][m]+difference.get(m,0) for m in (0,2,3,4,5)}
    print('Other missing prism pairing (01)|(24):',{m:v for m,v in last.items() if v})
    assert {m:v for m,v in last.items() if v} == {0:-F(1,28800),3:F(1,9600)}
    return rung,last


if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv)>1 else 0)
