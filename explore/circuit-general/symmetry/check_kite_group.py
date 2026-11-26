"""Orbit of a generic weight vector of K4 (= completed K5 with one vertex at infinity) under
vertex relabelling, completion/decompletion and Fourier (planar dual of K4: edge ij -> complementary edge),
together with the Gamma-factor bookkeeping: |orbit| of labelled weight vectors (cf. the Z2 x S6
symmetry, order 1440, of the two-loop master, Broadhurst 1986; Barfoot-Broadhurst 1988)."""
from itertools import permutations

K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
IDX = {e: i for i, e in enumerate(K4)}
q = 1000
a0 = (503, 487, 512, 491, 509, 498)        # generic, sum = 3000 = 3q


def relabels(a):
    for p in permutations(range(4)):
        b = [0] * 6
        for (u, v), x in zip(K4, a):
            e = tuple(sorted((p[u], p[v])))
            b[IDX[e]] = x
        yield tuple(b)


def fourier(a):
    comp = {(0, 1): (2, 3), (0, 2): (1, 3), (0, 3): (1, 2), (1, 2): (0, 3), (1, 3): (0, 2), (2, 3): (0, 1)}
    b = [0] * 6
    for e, x in zip(K4, a):
        b[IDX[comp[e]]] = q - x
    return tuple(b)


def decompletions(a):
    beta = {e: q - x for e, x in zip(K4, a)}
    for v in range(4):
        beta[(v, 4)] = sum(q - x for e, x in zip(K4, a) if v in e)
        beta[(v, 4)] = 2 * q - beta[(v, 4)]
    for c in range(4):
        verts = [x for x in range(5) if x != c]
        lab = {x: i for i, x in enumerate(verts)}
        b = [0] * 6
        for (u, v), w in beta.items():
            if c in (u, v):
                continue
            e = tuple(sorted((lab[u], lab[v])))
            b[IDX[e]] = q - w
        yield tuple(b)


seen = {a0}
st = [a0]
while st:
    a = st.pop()
    nxt = list(relabels(a)) + [fourier(a)] + list(decompletions(a))
    for b in nxt:
        assert sum(b) == 3 * q and min(b) > 0
        if b not in seen:
            seen.add(b)
            st.append(b)
print('labelled orbit size', len(seen))
seen2 = {a0}
st = [a0]
while st:
    a = st.pop()
    for b in list(relabels(a)) + list(decompletions(a)):
        if b not in seen2:
            seen2.add(b)
            st.append(b)
print('without Fourier', len(seen2))
seen3 = {a0}
st = [a0]
while st:
    a = st.pop()
    for b in list(relabels(a)) + [fourier(a)]:
        if b not in seen3:
            seen3.add(b)
            st.append(b)
print('without completion', len(seen3))
