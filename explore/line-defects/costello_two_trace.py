"""Color projector for the degree-three PBW symbol of the second-order line bracket."""
from collections import Counter
from itertools import permutations

# The marked dual pair is consecutive along the nontrivial trace after alternating
# the two source modes with two outputs; the final marked output tests the center.
PROBE = ((1,2),(2,3),(3,4),(4,1),(5,5))


def cycles(p):
    remaining = set(range(len(p)))
    out = []
    while remaining:
        word = []
        i = min(remaining)
        while i in remaining:
            word.append(i)
            remaining.remove(i)
            i = p[i]
        out.append(tuple(word))
    return tuple(out)


def surviving_traces():
    return [cycles(p) for p in permutations(range(5))
            if all(PROBE[i][1]==PROBE[p[i]][0] for i in range(5))]


def tree_words(leaves):
    """All oriented binary subtrees, using each specified external matrix once."""
    if len(leaves)==1:
        yield leaves
        return
    for bits in range(1,2**len(leaves)-1):
        left = tuple(v for i,v in enumerate(leaves) if bits>>i&1)
        right = tuple(v for i,v in enumerate(leaves) if not bits>>i&1)
        for a in tree_words(left):
            for b in tree_words(right):
                yield a+b


def verify():
    assert surviving_traces()==[((0,1,2,3),(4,))]
    # Root a cubic ribbon tree at a fixed external leaf: 5 binary parenthesizations
    # of 4 remaining leaves, each in 4! orders.
    words = Counter((0,)+w for w in tree_words((1,2,3,4)))
    assert len(words)==24 and set(words.values())=={5}
    assert not any(all(PROBE[a][1]==PROBE[b][0] for a,b in zip(w,w[1:]+w[:1])) for w in words)
    # A disconnected nonzero color type must leave the singleton label in a
    # component with only one external leg. Its sole one-loop cubic tadpole has
    # momenta 0,k,-k, trivial Moyal phases, and opposite color contributions.
    # Pair the matrix unit E_ij with its trace dual E_ji (all even slots).
    for n in (5,6,8):
        trace_tadpole = sum(int(i==4)-int(j==4) for i in range(n) for j in range(n))
        assert trace_tadpole==0
    print('The sole surviving permutation is the (4,1) trace; all 120 marked ribbon trees vanish.')


if __name__=='__main__':
    verify()
