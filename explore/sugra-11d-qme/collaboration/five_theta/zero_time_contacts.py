"""Exact common null vector of the sixty one-forms after collapsing a heat edge."""
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'verify'))
from vcommon import PQR,theta_graph
from replay import CASES,CUSTOM_VALUES,DEFAULT_AXES,SECTOR_VALUES


def coordinates(u,v,k,root):
    return [(sign,(w,k)) for w,sign in ((u,1),(v,-1)) if w!=root]


def product(x,y):
    terms=Counter()
    for a,i in x:
        for b,j in y:
            terms[tuple(sorted((i,j)))]+=a*b
    return terms


def check(shape):
    n,edges,_=theta_graph(*shape)
    cuts=0
    for cut,(first,second) in enumerate(edges):
        def merged(v):
            return first if v==second else v
        root=merged(n-1)
        vertices={merged(v) for v in range(n)}
        dimension=6*(len(vertices)-1)
        rows=0
        for e,(u,v) in enumerate(edges):
            if e==cut:
                continue
            u,v=merged(u),merged(v)
            assert u!=v  # No self-loop in these five contracted graphs.
            x=coordinates(u,v,0,root)
            for k in range(1,6):
                y=coordinates(u,v,k,root)
                # The five forms x_0 d q_k - x_k d q_0 span Ann(x)
                # on x_0 != 0, with x the edge's weighted Euler field.
                evaluation=product(x,y)
                evaluation.subtract(product(y,x))
                assert all(v==0 for v in evaluation.values())
                rows+=1
        assert rows==dimension==60 and len(vertices)==11
        cuts+=1
    assert cuts==13
    return cuts


if __name__=='__main__':
    counts={shape:check(shape) for shape in PQR}
    assert sum(counts.values())==65
    print('zero-time contacts in source edge order',counts)
    fixtures={(mask,DEFAULT_AXES):row for mask,row in SECTOR_VALUES.items()}
    fixtures.update(CUSTOM_VALUES)
    for (mask,axes),row in sorted(fixtures.items()):
        coefficient=-sum(x/case[1] for x,case in zip(row,CASES.values(),strict=True))/4096
        print(f'g^{(mask.bit_count()+1)//2} mask={mask:#x} axes={axes} '
              f'bare={coefficient}/pi^10 rooted={12*coefficient}/pi^10')
