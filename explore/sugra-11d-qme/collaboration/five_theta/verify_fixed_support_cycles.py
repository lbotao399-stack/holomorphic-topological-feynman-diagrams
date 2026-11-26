"""Exhibit an integer unit minor in the five-vector part of the CE boundary."""
from collections import Counter, defaultdict
from itertools import combinations, product
from math import prod

axes = range(5)
polarizations = list(product(*(tuple(j for j in axes if i != j) for i in axes)))


def shortest_cycle(f):
    cycles = set()
    for i in axes:
        trail = []
        while i not in trail:
            trail.append(i)
            i = f[i]
        loop = trail[trail.index(i):]
        cut = loop.index(min(loop))
        cycles.add(tuple(loop[cut:]+loop[:cut]))
    return min(cycles,key=lambda cycle:(len(cycle),cycle))


cycles = {f:shortest_cycle(f) for f in polarizations}
columns = sorted(polarizations,key=lambda f:(len(cycles[f]),f))
indices = {f:i for i,f in enumerate(columns)}


def row_key(f, i, j, component):
    return min(i,j),max(i,j),tuple(f[k] for k in axes if k not in (i,j)),component


def verify(signs):
    # For i<j the vector-field bracket is
    # signs[j]*delta[f[i],j]*e_(f[j]) - signs[i]*delta[f[j],i]*e_(f[i]).
    # A constant CE permutation sign per pair can be absorbed in its target row.
    rows = defaultdict(dict)
    for f in columns:
        for i,j in combinations(axes,2):
            if f[i]==j:
                rows[row_key(f,i,j,f[j])][f] = signs[j]
            if f[j]==i:
                rows[row_key(f,i,j,f[i])][f] = -signs[i]
    selected = set()
    diagonal = []
    for f in columns:
        cycle = cycles[f]
        i,j=cycle[:2]
        component=i if len(cycle)==2 else cycle[2]
        key=row_key(f,i,j,component)
        assert key not in selected
        selected.add(key)
        row=rows[key]
        assert abs(row[f])==1
        assert all(indices[g]<indices[f] for g in row if g!=f)
        diagonal.append(row[f])
    assert len(selected)==4**5 and abs(prod(diagonal))==1


if __name__ == "__main__":
    for signs in product((-1,1),repeat=5):
        verify(signs)
    print("32 signed axis choices, each with a triangular 1024 x 1024 minor of determinant ±1.")
    print("Pivot columns by shortest cycle length:",dict(sorted(Counter(map(len,cycles.values())).items())))
