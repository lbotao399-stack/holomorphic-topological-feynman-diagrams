"""Combine cyclic color numerators on each Kirchhoff polynomial and its variable-permutation orbits."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import lcm
from pathlib import Path
import json
import networkx as nx

from chain_moments import Reduction,graphs

def symmetries(red):
    bases=frozenset(k for k,v in red.minors.items() if v)
    g=nx.Graph()
    for i in range(9):
        g.add_node(i,count=sum(i in s for s in bases))
    for i,j in combinations(range(9),2):
        g.add_edge(i,j,count=sum(i in s and j in s for s in bases))
    matcher=nx.algorithms.isomorphism.GraphMatcher(g,g,node_match=lambda a,b:a==b,edge_match=lambda a,b:a==b)
    result=[]
    for p in matcher.isomorphisms_iter():
        if all(tuple(sorted(p[i] for i in s)) in bases for s in bases):
            result.append(tuple(p[i] for i in range(9)))
    return sorted(result)


def combine(start,end,folder):
    red=Reduction(graphs[start])
    group=symmetries(red)
    @cache
    def canonical(t):
        return min(tuple(t[i] for i in p) for p in group)
    totals=defaultdict(lambda:defaultdict(F))
    ranks={}
    for idx in range(start,end+1):
        rows=json.loads((folder/f'graph{idx:02}.json').read_text())
        ranks[idx]=[]
        for row in rows:
            assert row['edges']==[list(e) for e in red.core]
            local=defaultdict(F)
            for t,coeff in row['terms']:
                full=tuple(a+b for a,b in zip(t,row['h']))
                key=canonical(full)
                local[key]+=F(coeff,row['denominator'])
            if any(local.values()):
                ranks[idx].append(row['N_power'])
                for key,v in local.items():
                    totals[row['N_power']][key]+=v
    out=[]
    for k,coeff in sorted(totals.items()):
        coeff={t:c for t,c in coeff.items() if c}
        if not coeff:continue
        den=lcm(*(v.denominator for v in coeff.values()))
        out.append({'graphs':[start,end],'edges':red.core,'h':[0]*9,'N_power':k,
                    'symmetries':len(group),'orbit_representatives':True,'denominator':den,
                    'terms':[[list(t),int(c*den)] for t,c in sorted(coeff.items())]})
    target=folder/f'core{start:02}-{end:02}.json'
    target.write_text(json.dumps(out,separators=(',',':'))+'\n')
    return out,ranks


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('first',type=int)
    parser.add_argument('last',type=int)
    args=parser.parse_args()
    rows,ranks=combine(args.first,args.last,args.folder)
    for row in rows:
        signs=Counter(c>0 for _,c in row['terms'])
        print(row['N_power'],row['symmetries'],len(row['terms']),row['denominator'],dict(signs))
    print('surviving ranks by graph',ranks)
