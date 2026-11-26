"""Odd-coordinate constraints on the N=4 fan and the interaction trees in a marked product."""
from itertools import combinations, product


# Bit i records the presence of theta_{i+1} in the functional derivative.
letters = {'c':7,'gamma1':6,'gamma2':5,'gamma3':3,'beta1':1,'beta2':2,'beta3':4,'b':0}


def wedge(left,right):
    result = {}
    for a,c in left.items():
        for b,d in right.items():
            if a & b:
                continue
            parity = sum((a>>(i+1)).bit_count() for i in range(b.bit_length()) if b & (1<<i))
            key = a | b
            result[key] = result.get(key,0)+(-1)**parity*c*d
    return {a:c for a,c in result.items() if c}


def edge_delta(i,j):
    ans = {0:1}
    for k in range(3):
        ans = wedge(ans,{1<<(3*i+k):1,1<<(3*j+k):-1})
    return ans


def disjoint_count(present):
    """Maximum number of non-b inputs of a single fan allowed by its odd integral."""
    betas = [i for i in range(1,4) if 'beta'+str(i) in present]
    if len(betas)==3:
        return 3
    if len(betas)>=2 or any('gamma'+str(i) in present for i in betas):
        return 2
    return int(any(s!='b' for s in present))


def check():
    names = list(letters)
    allowed = []
    for mask in range(1<<len(names)):
        group = {x for i,x in enumerate(names) if mask & (1<<i)}
        if all(letters[x] & letters[y] for x,y in product(group,repeat=2)):
            allowed.append(group)
        nonempty = [letters[x] for x in group if x!='b']
        rank = max(k for k in range(len(nonempty)+1) if any(
            all(x&y==0 for x,y in combinations(choice,2)) for choice in combinations(nonempty,k)))
        assert disjoint_count(group)==rank
    maximal = [g for g in allowed if not any(g<h for h in allowed)]
    assert len(maximal)==4 and all(len(g)==4 for g in maximal)
    # Work on the unreduced odd coordinates of every vertex, including the
    # four-vertex branching tree. An allowed pair must really survive the deltas.
    for tree in ([],[(0,1)],[(0,1),(1,2)],[(0,1),(0,2),(0,3)]):
        kernel = {0:1}
        for edge in tree:
            kernel = wedge(kernel,edge_delta(*edge))
        vertices = range(len(tree)+1)
        for i,j in product(vertices,repeat=2):
            for x,y in product(letters.values(),repeat=2):
                value = wedge(wedge(kernel,{x<<(3*i):1}),{y<<(3*j):1})
                assert bool(value)==(x&y==0)
    print("Four maximal coordinate sectors and all 256 fan supports verified.")


if __name__ == '__main__':
    check()
