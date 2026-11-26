"""Ordered odd fibre contractions and automorphism orders of the rooted cubic graphs."""
from collections import Counter
from fractions import Fraction
from itertools import chain
from math import factorial, prod
import networkx as nx

from rooted_trees import graphs, symmetry_weight


def wick_sign(edges, root):
    """Remove the tail, then the head of each oriented edge with left derivatives."""
    word = []
    for v in range(root):
        incoming = [(i,v) for i,(_,head) in enumerate(edges) if head == v]
        outgoing = [(i,v) for i,(tail,_) in enumerate(edges) if tail == v]
        for slot in range(3-len(incoming)-len(outgoing)):
            incoming.append(('output',v,slot))
        word.extend(incoming + outgoing)
    spokes = [(i,root) for i,edge in enumerate(edges) if root in edge]
    word.extend(spokes)
    # A printed product of m left derivatives acts on the input word right to left.
    crossings = len(spokes)*(len(spokes)-1)//2
    for i, edge in enumerate(edges):
        for v in edge:
            position = word.index((i,v))
            crossings += position
            word.pop(position)
    return (-1)**crossings


def automorphisms(edges, root):
    """Automorphisms fixing the insertion are the mark-preserving tree automorphisms."""
    tree = nx.Graph()
    tree.add_nodes_from(range(root))
    tree.add_edges_from(e for e in edges if root not in e)
    marked = {v for e in edges if root in e for v in e if v != root}

    def branch(v, parent):
        children = sorted(branch(w,v) for w in tree[v] if w != parent)
        types = tuple(kind for kind, _ in children)
        multiplicity = prod(factorial(n) for n in Counter(types).values())
        return (v in marked, types), multiplicity*prod(count for _,count in children)

    center = nx.center(tree)
    if len(center) == 1:
        return branch(center[0],None)[1]
    left, right = center
    kind_left, a = branch(left,right)
    kind_right, b = branch(right,left)
    return a*b*(2 if kind_left == kind_right else 1)


def verify():
    from pathlib import Path
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parent/'verify'))
    from random import Random
    from vops import GRAPH_SPECS
    expected = (Fraction(1,2), Fraction(1,2), Fraction(1,2),
                Fraction(-1,2), Fraction(1,2), -1, Fraction(1,8), Fraction(-1,2))
    names = ('box','A','B','G1','G2','G3','G4','G5')
    for name, known in zip(names,expected):
        labels = GRAPH_SPECS[name]['edges']
        root = len(set(chain.from_iterable(labels)))-1
        def vertex(v):
            return root if v == 'O' else v-1
        edges = [(vertex(a),vertex(b)) for a,b in labels]
        assert Fraction(wick_sign(edges,root),automorphisms(edges,root)) == known
    assert wick_sign([(0,1)],1) == 1
    # Rebuild the ordered color slots and the two lists of odd external factors.
    def permutation_sign(order):
        return (-1)**sum(a > b for j,a in enumerate(order) for b in order[j+1:])

    def local_orders(edges, root):
        inputs, outputs, vertices = [], [], []
        for v in range(root):
            ins = [k for k,(a,b) in enumerate(edges) if b == v]
            outs = [k for k,(a,b) in enumerate(edges) if a == v]
            local = ins + outs
            if len(local) == 2:
                outputs.append(v)
                local = ins + [-1] + outs
            vertices.append(local)
        inputs = [k for k,e in enumerate(edges) if root in e]
        return inputs, outputs, vertices

    def relabelling(edges, root, rng):
        perm = list(range(root))
        order = list(range(len(edges)))
        rng.shuffle(perm)
        rng.shuffle(order)
        reverse = [rng.randrange(2) for _ in order]
        new = []
        for j,k in enumerate(order):
            edge = [root if v == root else perm[v] for v in edges[k]]
            new.append(tuple(edge[::-1] if reverse[j] else edge))
        old_inputs, old_outputs, old_legs = local_orders(edges,root)
        inputs, outputs, legs = local_orders(new,root)
        tensor = permutation_sign([old_inputs.index(order[j]) for j in inputs])
        tensor *= permutation_sign([old_outputs.index(perm.index(v)) for v in outputs])
        for v in range(root):
            updated = [order[j] if j != -1 else -1 for j in legs[perm[v]]]
            tensor *= permutation_sign([old_legs[v].index(j) for j in updated])
        kernel = permutation_sign(perm)*(-1)**sum(reverse)
        assert wick_sign(new,root)*tensor*kernel == wick_sign(edges,root)

    rng = Random(115)
    for n in range(1,7):
        for edges in graphs(n):
            relabelling(edges,2*n+1,rng)
    for n in range(1,8):
        root = 2*n+1
        orders = [automorphisms(edges,root) for edges in graphs(n)]
        assert sum(Fraction(1,g) for g in orders) == symmetry_weight(n)
        print(n, 'automorphism orders and counts', sorted(Counter(orders).items()))


if __name__ == '__main__':
    verify()
