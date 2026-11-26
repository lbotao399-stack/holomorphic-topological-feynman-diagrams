"""Two-loop open SO(N) theta graphs, in the repository's reduced-kernel and Wick conventions."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations,permutations,product
from math import factorial,lcm,prod
from pathlib import Path
import json
import sys

from hexagon import parity
from htk.linalg import det
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'twoloop-universal'))
import twoloop
SHAPES=((1,5,5),(2,4,5),(3,3,5),(3,4,4))
HERE=Path(__file__).resolve().parent


def theta(shape):
    return twoloop.Theta(*twoloop.theta_graph(*shape))


def sign_of_wick(shape):
    graph=theta(shape)
    half_edges=[]
    for e,(u,v) in enumerate(graph.edges):
        half_edges.extend([(e,u),(e,v)])
    external=[('x',v) for v in range(1,9)]
    slots=external+half_edges
    local=[]
    for v in range(10):
        if 0<v<9:
            local.append(('x',v))
        local.extend(h for h in half_edges if h[1]==v)
    return parity([slots.index(h) for h in local])


def cyclic(word):
    for i in range(len(word)):
        yield word[i:]+word[:i]


def canonical_trace(word):
    if len(word)==1:
        return (),0
    if not word:
        return (),1
    normal=min(cyclic(word))
    reverse=min(cyclic(word[::-1]))
    if reverse<normal:
        return reverse,(-1)**len(word)
    return normal,1


def rotate_to(words,h):
    idx=next(i for i,w in enumerate(words) if h in w)
    w=words[idx]
    j=w.index(h)
    return idx,w[j+1:]+w[:j]


def colour(shape):
    graph=theta(shape)
    # Negative letters are the two copies of each of the three adjoint basis elements.
    first=(-1,-3,-5)
    counts=Counter()
    for sides in product((0,1),repeat=8):
        second=[]
        for path,key in zip('ABC',(-2,-4,-6)):
            word=(key,)
            for u in graph.interior[path]:
                i=u-1
                word=(i,)+word if not sides[i] else word+(i,)
            second.extend(word)
        pending=[((first,tuple(second)),(-1)**sum(sides))]
        for a,b in zip(first,(-2,-4,-6)):
            result=[]
            for words,c in pending:
                i,u=rotate_to(words,a)
                j,v=rotate_to(words,b)
                if i!=j:
                    rest=tuple(w for k,w in enumerate(words) if k not in (i,j))
                    result.extend([(rest+(u+v,),c),(rest+(u[::-1]+v,),-(-1)**len(u)*c)])
                else:
                    cut=u.index(b)
                    u,v=u[:cut],u[cut+1:]
                    rest=tuple(w for k,w in enumerate(words) if k!=i)
                    result.extend([(rest+(u,v),c),(rest+(u[::-1]+v,),-(-1)**len(u)*c)])
            pending=result
        for words,c in pending:
            canonical=[canonical_trace(w) for w in words]
            traces=tuple(sorted((w for w,_ in canonical if w),key=lambda w:(len(w),w)))
            counts[len([w for w in words if not w]),traces]+=c*prod(x for _,x in canonical)
    # Both trivalent endpoints contribute 2 tr(ABC); each of the three Fierz contractions costs 1/2.
    return {k:F(sign_of_wick(shape)*v,2) for k,v in counts.items() if v}


def geometric(shape):
    g=theta(shape)
    a,b,c=([i-1 for i in g.interior[k]] for k in 'ABC')
    weights=[twoloop.mu_weights(g)[v] for v in range(1,9)]
    factor=twoloop.c0(shape,5)*twoloop.sgn_detBC(10,g.edges,9,g.C3())**5
    ans=Counter()
    for upper_b in combinations(b,4-len(a)):
        lower_b=[k for k in b if k not in upper_b]
        upper=a+list(upper_b)
        lower=lower_b+c
        selection=[0]+[2+i for i in upper]+[1]+[2+i for i in lower]
        for x,y in product(range(8),repeat=2):
            if x in upper or y in lower:
                continue
            row1=(x,*upper)
            row2=(y,*lower)
            key=tuple(sorted((tuple(sorted(row1)),tuple(sorted(row2)))))
            ans[key]+=factor*weights[x]*weights[y]*parity(selection)*parity(row1)*parity(row2)
    return {k:v for k,v in ans.items() if v}


def verify_geometry():
    from random import Random
    r=Random(25532)
    for shape in SHAPES:
        geom=geometric(shape)
        g=theta(shape)
        for _ in range(3):
            p=[[r.randrange(-8,9) for _ in range(5)] for _ in range(7)]
            p.append([-sum(c) for c in zip(*p)])
            minors={s:det([p[i] for i in s]) for s in combinations(range(8),5)}
            value=sum(v*minors[a]*minors[b] for (a,b),v in geom.items())
            expected=twoloop.kernel_htk_closed(g,[[0]*5,*p,[0]*5],5,0)
            assert value==expected,(shape,value,expected)


def relabelings(traces):
    lengths=list(map(len,traces))
    target=[]
    for i,n in enumerate(lengths):
        first=sum(lengths[:i])
        target.append(tuple(range(first,first+n)))
    for assignment in permutations(range(len(target))):
        if any(len(target[j])!=n for n,j in zip(lengths,assignment)):
            continue
        possibilities=[]
        for j in assignment:
            w=target[j]
            orients={x:1 for x in cyclic(w)}
            orients.update({x:(-1)**len(w) for x in cyclic(w[::-1])})
            possibilities.append(list(orients.items()))
        for choices in product(*possibilities):
            perm=[None]*8
            for src,(dst,_) in zip(traces,choices):
                for i,j in zip(src,dst):
                    perm[i]=j
            yield tuple(perm),parity(perm)*prod(sign for _,sign in choices)


@cache
def reduce_minor(seq):
    """Replace p7 by momentum conservation, then pass to complementary two-column labels."""
    reduced=Counter()
    for k in range(7) if 7 in seq else [None]:
        row=tuple(k if i==7 else i for i in seq)
        if len(set(row))!=5:
            continue
        complement=tuple(i for i in range(7) if i not in row)
        reduced[complement]+=(-1 if k is not None else 1)*parity(row+complement)
    return tuple(reduced.items())


@cache
def straight(a,b):
    if a>b:
        return straight(b,a)
    i,j=a
    k,l=b
    if j>l:
        return ((((i,l),(k,j)),1),(((i,k),(l,j)),-1))
    return (((a,b),1),)


def graph_coefficients(shape):
    graph=geometric(shape)
    denominator=lcm(*(v.denominator for v in graph.values()))
    graph=[(a,b,int(denominator*c)) for (a,b),c in graph.items()]
    color=colour(shape)
    symmetry=2*prod(factorial(shape.count(n)) for n in set(shape))
    output=defaultdict(Counter)
    # Combine repeated color assignments before acting on the two five-minors.
    perms=defaultdict(Counter)
    for (rank,traces),c in color.items():
        for perm,sign in relabelings(traces):
            perms[(rank,tuple(map(len,traces)))][perm]+=int(2*c)*sign
    for target,maps in perms.items():
        raw=Counter()
        for perm,count in maps.items():
            if not count:continue
            for a,b,c in graph:
                x=tuple(perm[i] for i in a)
                y=tuple(perm[i] for i in b)
                for u,i in reduce_minor(x):
                    for v,j in reduce_minor(y):
                        raw[tuple(sorted((u,v)))]+=count*c*i*j
        for pair,c in raw.items():
            for key,sign in straight(*pair):
                output[target][key]+=c*sign
    return {k:{m:F(v,2*denominator*symmetry) for m,v in poly.items() if v}
            for k,poly in output.items()}


if __name__=='__main__':
    verify_geometry()
    total=defaultdict(Counter)
    all_graphs={}
    for shape in SHAPES:
        values=graph_coefficients(shape)
        print(shape,'color words',len(colour(shape)),'surviving coefficients',[(k,len(v)) for k,v in values.items() if v],flush=True)
        for k,row in values.items():
            total[k].update(row)
        all_graphs[''.join(map(str,shape))]=values
    all_graphs['sum']={k:{m:v for m,v in row.items() if v} for k,row in total.items()}
    def encode(data):
        return [dict(rank=rank,lengths=list(l),terms=[[list(a),list(b),str(v)] for (a,b),v in sorted(poly.items())])
                for (rank,l),poly in sorted(data.items()) if poly]
    result={shape:encode(data) for shape,data in all_graphs.items()}
    path=HERE/'open_theta.json'
    if '--save' in sys.argv:
        path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    assert result==json.loads(path.read_text())
    print('sum',[(k,len(v)) for k,v in all_graphs['sum'].items() if v])
