"""Least total current-insertion degree on spokes of each cubic product forest."""
from collections import Counter
from itertools import combinations_with_replacement
from forest_jets import chart, integrate_factors
from forest_orders import forests,graph,symbolic_chart


def odd_swaps(forest):
    # Vertex labels agree with graph(): preorder of the trees, then A and B.
    next_label=0
    swaps=[]
    def parse(t,pos):
        nonlocal next_label
        if t[pos] in 'AB':
            return pos+1
        next_label+=1
        begin=next_label
        middle=parse(t,pos+1)
        half=next_label
        end=parse(t,middle)
        if t[pos+1:middle]==t[middle:end]:
            assert half>begin
            swaps.append(set(range(begin,next_label)))
        return end+1
    for t in forest:
        parse(t,0)
    return swaps


def least_degree(forest,limit):
    edges,a,b,roots,_=graph(forest)
    swaps=odd_swaps(forest)
    spokes=[e for e in edges if e[0]<a<=e[1]]
    source=[(0,0)]*a+[(1,0),(-1,-len(roots))]
    for v in roots:
        source[v]=(0,1)
    numerical=chart(edges,a,b,source,True)
    symbolic=None
    def witness(words,form):
        factors,current,simplices=form
        for word in words:
            value=integrate_factors(factors+[current[spokes[idx//2]][idx%2] for idx in word],simplices)
            if value:
                return True
        return False
    for degree in range(limit+1):
        # A leaf transposition untouched by every insertion still negates the integral.
        words=[w for w in combinations_with_replacement(range(2*len(spokes)),degree)
               if all(any(spokes[j//2][0] in s for j in w) for s in swaps)]
        if not words:
            continue
        if witness(words,numerical):
            return degree
        if len(roots)>1:
            if symbolic is None:
                symbolic=symbolic_chart(edges,a,b,roots)
            if witness(words,symbolic):
                return degree
    return None


def run(n):
    counts=Counter()
    for f in forests(n):
        d=least_degree(f,2*n)
        if d is None:
            raise ArithmeticError(f'Missing the degree-{2*n} cofactor witness for {f}')
        counts[d]+=1
    print('loops',n,'forests by first nonzero insertion degree',dict(counts),flush=True)
    return counts


if __name__=='__main__':
    import sys
    for n in range(1,int(sys.argv[1])+1 if len(sys.argv)>1 else 8):
        run(n)
