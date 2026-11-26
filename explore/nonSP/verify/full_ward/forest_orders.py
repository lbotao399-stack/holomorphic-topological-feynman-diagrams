"""Enumerate marked-product forests and classify the zero-displacement kernels."""
from functools import cache
from flint import fmpq_mpoly_ctx
from forest_jets import chart, edge, evaluate, integrate_factors


@cache
def trees(n):
    # A,B denote the two ordered insertions. Parentheses enclose unordered child ports.
    if n==0:
        return ('A','B')
    result=[]
    for size in range((n-1)//2+1):
        for l in trees(size):
            for r in trees(n-1-size):
                if size==n-1-size and (l>r or (size==0 and l==r)):
                    continue
                result.append('('+l+r+')')
    return tuple(result)


def forests(n,smallest=(1,'')):
    if n==0:
        yield ()
    for size in range(smallest[0],n+1):
        for tree in trees(size):
            if (size,tree)<smallest:
                continue
            for tail in forests(n-size,(size,tree)):
                yield (tree,)+tail


def graph(forest):
    n=sum(t.count('(') for t in forest)
    a,b=n,n+1
    edges=[(a,b)]
    roots=[]
    next_label=0
    odd=False
    def parse(text,pos):
        nonlocal next_label,odd
        if text[pos] in 'AB':
            return (a if text[pos]=='A' else b),pos+1
        v=next_label
        next_label+=1
        l,middle=parse(text,pos+1)
        r,end=parse(text,middle)
        if text[pos+1:middle]==text[middle:end]:
            odd=True
        assert text[end]==')'
        edges.extend((edge(v,l),edge(v,r)))
        return v,end+1
    for t in forest:
        v,end=parse(t,0)
        assert end==len(t)
        roots.append(v)
    return tuple(sorted(edges)),a,b,roots,odd


def symbolic_chart(edges,a,b,roots):
    # On the open set det(lambda_A,lambda_root0) != 0, a change of basis
    # sends these two momenta to (1,0),(0,1). All other roots remain arbitrary.
    n=(len(edges)-1)//2
    names=tuple(f'r{i}' for i in range(2*n))+tuple(
        f'p{v}_{j}' for v in roots[1:] for j in range(2))
    ctx=fmpq_mpoly_ctx.get(names)
    source=[(0,0)]*(n+2)
    source[a]=(1,0)
    source[roots[0]]=(0,1)
    variables=ctx.gens()[2*n:]
    for i,v in enumerate(roots[1:]):
        source[v]=variables[2*i:2*i+2]
    source[b]=tuple(-sum(row[j] for row in source) for j in (0,1))
    return chart(edges,a,b,source,True,ctx)


def source_polynomial(edges,a,b,roots):
    factors,_,simplices=symbolic_chart(edges,a,b,roots)
    return integrate_factors(factors,simplices)


def run(n):
    total,symmetric,nonzero,other_zero=0,0,0,[]
    for f in forests(n):
        edges,a,b,roots,odd=graph(f)
        total+=1
        if odd:
            # Exchanging isomorphic nonempty branches permutes 2s+1 edge pairs
            # for s internal vertices and fixes all allowed vertex momenta.
            symmetric+=1
            continue
        source=[(0,0)]*(n+2)
        source[a]=(1,0)
        for v in roots:
            source[v]=(0,1)
        source[b]=(-1,-len(roots))
        value=evaluate(edges,a,b,source)
        if value or (len(roots)>1 and source_polynomial(edges,a,b,roots)):
            nonzero+=1
        else:
            other_zero.append(f)
    print(n,total,symmetric,nonzero,'additional identically zero forests',other_zero,flush=True)
    return total,symmetric,nonzero,other_zero


if __name__=='__main__':
    import sys
    for i in range(1,int(sys.argv[1]) + 1 if len(sys.argv)>1 else 8):
        run(i)
