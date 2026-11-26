"""Direct Koszul ranks on small transverse jet spaces, without Weyl characters."""
from itertools import combinations

import sympy as sp


def ordered(items):
    if len(set(items)) != len(items):
        return None, 0
    sign = (-1)**sum(x>y for i,x in enumerate(items) for y in items[i+1:])
    return tuple(sorted(items)),sign


def increment(jet,i):
    axis, degree = jet
    exps = list(degree)
    exps[i] += 1
    if axis==0 and exps[0]:
        # In two dimensions the differentiated divergence has two terms.
        exps[0] -= 1
        exps[1] += 1
        return (1,tuple(exps)),-1
    return (axis,tuple(exps)),1


def basis(q,w,j):
    degree = 2*w-j
    jets = []
    for h in range(1,degree-q+2):
        for a in range(h+2):
            for axis in range(2):
                if a and axis==0:
                    continue
                jets.append((axis,(a,h+1-a)))
    result = []
    for mon in combinations(sorted(jets),q):
        if sum(sum(degree)-1 for axis,degree in mon) != degree:
            continue
        for axes in combinations(range(2),j):
            if all(sum(int(axis==i)-d[i] for axis,d in mon)-int(i in axes)==-w for i in range(2)):
                result.append((mon,axes))
    return result


def differential(domain,codomain):
    rows = {entry:i for i,entry in enumerate(codomain)}
    matrix = sp.zeros(len(codomain),len(domain))
    for col,(mon,axes) in enumerate(domain):
        for k,i in enumerate(axes):
            remainder = axes[:k]+axes[k+1:]
            for l,jet in enumerate(mon):
                new,sgn = increment(jet,i)
                term,sgn2 = ordered(mon[:l]+(new,)+mon[l+1:])
                if sgn2:
                    matrix[rows[term,remainder],col] += (-1)**k*sgn*sgn2
    return matrix


def check(q,w):
    spaces = [basis(q,w,j) for j in range(3)]
    first = differential(spaces[1],spaces[0])
    second = differential(spaces[2],spaces[1])
    assert first*second == sp.zeros(len(spaces[0]),len(spaces[2]))
    assert second.rank()==len(spaces[2])
    assert first.rank()==len(spaces[1])-len(spaces[2])
    print(q,w,tuple(map(len,spaces)),(first.rank(),second.rank()))


if __name__=='__main__':
    for q,w in ((2,2),(2,3),(3,2),(3,3)):
        check(q,w)
