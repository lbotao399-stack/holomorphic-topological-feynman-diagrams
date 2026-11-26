"""Exact evaluation and soft tests of polarized primary scalar counterterms."""
from functools import cache
from itertools import combinations, permutations
from subprocess import run

from higher_genus_two import add, det


@cache
def exterior_shapes(n,free):
    # Rows sharing an odd number of determinants cancel under their exchange.
    masks=[[sum(1<<perm[i] for i in range(4) if h>>i&1) for h in range(16)]
           for perm in permutations(range(4))]
    def fills(h, left, degrees):
        if h==16:
            if left==0 and not any(degrees):
                yield ()
            return
        axes=[i for i in range(4) if h>>i&1]
        upper=min(left,*(degrees[i] for i in axes))
        if h.bit_count()%2:
            upper=min(1,upper)
        for m in range(upper+1):
            rest=tuple(x-m*(i in axes) for i,x in enumerate(degrees))
            if min(rest)<0 or max(rest)>left-m:
                continue
            for tail in fills(h+1,left-m,rest):
                yield (h,)*m+tail
    out=[]
    for h0 in ((0,1,3,7,15) if free else (0,)):
        stabilizer=[p for p in masks if p[h0]==h0]
        for tail in fills(1,n-int(free),tuple(5-(h0>>i&1) for i in range(4))):
            if any(tuple(sorted(p[h] for h in tail))<tail for p in stabilizer):
                continue
            out.append(((16+h0,) if free else ())+tail)
    return tuple(out)


def flags_of(word):
    subsets=list(combinations(range(9),5))
    return tuple(sum(1<<j for j,s in enumerate(word) if i in subsets[s]) for i in range(9))


def orbit_weight(flags):
    # Ordered incidence assignments per distinct signed monomial.
    from collections import Counter
    from math import prod,factorial
    groups=Counter(flags)
    stabilizer=0
    for perm in permutations(range(4)):
        changed=[(h&16)+sum(1<<perm[i] for i in range(4) if h>>i&1) for h in flags]
        if Counter(changed)==groups:
            stabilizer+=1
    edges=[tuple(i for i,h in enumerate(flags) if h>>j&1) for j in range(4)]
    return stabilizer//prod(factorial(m) for m in Counter(edges).values())


def native(exe,shapes,tests,modular=False):
    ns=len(tests[0])
    subsets=list(combinations(range(ns),5))
    data=[(ns,len(tests),len(shapes),int(modular)),*((len(s),*s) for s in shapes)]
    for p in tests:
        data.extend(p)
        data.append([det([p[i] for i in q]) for q in subsets])
    stdin='\n'.join(' '.join(map(str,row)) for row in data)
    result=run([exe],text=True,input=stdin,capture_output=True,check=True)
    return [tuple(map(int,row.split())) for row in result.stdout.splitlines()]


def fused(exe,shapes,p):
    pairs=list(combinations(range(len(p)),2))
    tests=[[p[h] for h in range(len(p)) if h not in (i,j)]+[add(p[i],p[j])] for i,j in pairs]
    canonical=[tuple(sorted(q)) for q in tests]
    unique=list(dict.fromkeys(canonical))
    back={q:i for i,q in enumerate(unique)}
    short=native(exe,shapes,unique)
    coeffs=[[row[back[q]] for q in canonical] for row in short]
    out=[[[0]*5 for row in p] for s in shapes]
    for at,values in enumerate(coeffs):
        for (i,j),v in zip(pairs,values):
            for k in range(5):
                out[at][i][k]+=v*p[j][k]
                out[at][j][k]+=v*p[i][k]
    return out


def soft_test(exe,shapes,p):
    tests=[[s for j,s in enumerate(p) if j!=i] for i in range(len(p))]
    cols=native(exe,shapes,tests,modular=True)
    rows=[]
    for i,k in enumerate(p):
        pivot=next(a for a,x in enumerate(k) if x)
        for a in range(5):
            if a!=pivot:
                rows.append([(k[pivot]*v[5*i+a]-k[a]*v[5*i+pivot])%1009 for v in cols])
    return rows
