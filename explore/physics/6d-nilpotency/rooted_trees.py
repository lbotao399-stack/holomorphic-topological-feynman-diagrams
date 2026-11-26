"""All-loop counts for the cubic rooted C^3|2 graph selection rules."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product as cartesian_product
from math import comb

import networkx as nx


def product(left, right, degree):
    return sum(left[i]*right[degree-i] for i in range(degree+1))


def unordered_pair(branch, degree):
    fixed = branch[degree//2] if degree % 2 == 0 else 0
    return (product(branch, branch, degree) + fixed)//2


def counts(order):
    # incoming: the parent edge points into this vertex. outgoing: it points out.
    incoming = [1] + [0]*order
    outgoing = [0]*(order+1)
    result = [1]
    for n in range(1, order+1):
        outgoing[n] = incoming[n-1] + product(incoming, outgoing, n-1)
        incoming[n] = outgoing[n] + unordered_pair(incoming, n-1)
        # Mark a vertex or an undirected edge, then subtract the incidences.
        # Each tree edge already has a canonical direction.
        sinks = outgoing[n] + unordered_pair(outgoing, n)
        splitters = unordered_pair(incoming, n-1) + sum(
            unordered_pair(incoming, n-1-i)*outgoing[i] for i in range(1,n))
        result.append(sinks + splitters - product(incoming, outgoing, n))
    return result


def symmetry_weight(n):
    if n == 0:
        return Fraction(1)
    return sum(Fraction(comb(n+2+i, i) * comb(2*n-2-i, n-1-i), 2**(n-1-i)*n*(n+1))
               for i in range(n))


def paired_branches(make, n):
    for i in range(n//2+1):
        if 2*i == n:
            yield from combinations_with_replacement(make(i), 2)
        else:
            yield from cartesian_product(make(i), make(n-i))


@cache
def incoming_branches(n):
    if n == 0:
        return ((True, ()),)
    marked = [(True, (b,)) for b in outgoing_branches(n)]
    unmarked = [(False, tuple(sorted(ab))) for ab in paired_branches(incoming_branches, n-1)]
    return tuple(sorted(marked + unmarked))


@cache
def outgoing_branches(n):
    if n == 0:
        return ()
    trees = [(False, (a,)) for a in incoming_branches(n-1)]
    for i in range(n-1):
        trees.extend((False, tuple(sorted((a,b)))) for a,b in
                     cartesian_product(incoming_branches(i), outgoing_branches(n-1-i)))
    return tuple(sorted(trees))


def graphs(n):
    """Yield one ordered edge list per class, with the operator at vertex 2*n+1."""
    roots = ((True, (b,)) for b in outgoing_branches(n))
    pairs = ((True, tuple(sorted(ab))) for ab in paired_branches(outgoing_branches,n))
    from itertools import chain
    for tree in chain(roots,pairs, ((True, ()),) if n == 0 else ()):
        neighbours = []
        marked = []
        edges = []
        def visit(node, parent):
            vertex = len(neighbours)
            neighbours.append([])
            if parent is not None:
                edges.append((parent,vertex))
                neighbours[parent].append(vertex)
                neighbours[vertex].append(parent)
            if node[0]:
                marked.append(vertex)
            for child in node[1]:
                visit(child,vertex)
        visit(tree,None)
        marks = set(marked)
        @cache
        def encoding(v,parent):
            return (v in marks, tuple(sorted(encoding(u,v) for u in neighbours[v] if u != parent)))
        if tree == min(encoding(v,None) for v in marked):
            yield edges + [(v,len(neighbours)) for v in marked]
        encoding.cache_clear()


def decorations(n):
    # Independent unlabelled-tree enumeration. Test every connected subtree,
    # deduplicate by the full automorphism group of the unmarked tree.
    for tree in nx.nonisomorphic_trees(2*n+1):
        if max(dict(tree.degree()).values()) > 3:
            continue
        leaves = {v for v, d in tree.degree if d == 1}
        internal = [v for v, d in tree.degree if d == 2]
        extra = n+1-len(leaves)
        if not 0 <= extra <= len(internal):
            continue
        auto = list(nx.isomorphism.GraphMatcher(tree,tree).isomorphisms_iter())
        remaining = {frozenset(leaves.union(c)) for c in combinations(internal, extra)}
        subtrees = [set(c) for k in range(2,len(tree)) for c in combinations(tree,k)
                    if nx.is_connected(tree.subgraph(c))]
        while remaining:
            marked = next(iter(remaining))
            orbit = {frozenset(g[v] for v in marked) for g in auto}
            remaining.difference_update(orbit)
            if all(2*len(marked & sub) <= len(sub)+1 for sub in subtrees):
                yield len(auto)//len(orbit)


if __name__ == '__main__':
    numbers = counts(15)
    for n in range(1,5):
        symmetries = Counter(decorations(n))
        assert sum(symmetries.values()) == numbers[n]
        assert sum(Fraction(k,g) for g,k in symmetries.items()) == symmetry_weight(n)
        print(n, 'automorphism order and multiplicity', sorted(symmetries.items()))
    for n in range(8):
        assert sum(1 for _ in graphs(n)) == numbers[n]
    print('unlabelled counts, n=0..15:', numbers)
    print('symmetry weights, n=0..6:', *(symmetry_weight(n) for n in range(7)))
