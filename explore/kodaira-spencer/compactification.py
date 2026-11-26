"""Ring factor of a closed connected bulk graph for the K3 and four-torus reductions."""
from itertools import product


def torus_product(x,y):
    if x & y:
        return 0, 0
    inversions = sum((y & ((1 << j)-1)).bit_count() for j in range(4) if x & (1 << j))
    return (-1)**inversions, x|y


def k3_product(x,y):
    # Basis ordering: unit, eta, eta-bar, 20 orthogonal middle classes, top.
    if x == 0 or y == 0:
        return 1, x+y
    if {x,y} == {1,2} or (x == y and 3 <= x <= 22):
        return 1, 23
    return 0, 0


def cubic(multiply, top, x,y,z):
    first, xy = multiply(x,y)
    second, xyz = multiply(xy,z)
    return first*second if xyz == top else 0


def audit(multiply, top, parity):
    labels = range(top+1)
    dual = []
    handle = [0 for _ in labels]
    for i in labels:
        (j,) = [j for j in labels if cubic(multiply,top,i,j,0)]
        coefficient = cubic(multiply,top,i,j,0)
        dual.append(j)
        sign, ij = multiply(i,j)
        handle[ij] += (-1)**parity(i)*coefficient*sign
    assert not any(handle[:-1])
    # Direct matter labels on the two trivalent ends of a vacuum theta.
    # Insertions on a path multiply its dual endpoint label by the external labels.
    for i,j,k in product(labels,repeat=3):
        if cubic(multiply,top,i,j,k):
            assert cubic(multiply,top,dual[i],dual[j],dual[k]) == 0
    for label in labels:
        supertrace = sum((-1)**parity(j)*multiply(label,j)[0] for j in labels
                         if multiply(label,j)[1] == j)
        assert supertrace == handle[top]*cubic(multiply,top,top,label,0)
    return handle[top]


if __name__ == '__main__':
    assert audit(k3_product,23,lambda x:0) == 24
    assert audit(torus_product,15,lambda x:x.bit_count()%2) == 0
    print('ring-valued bulk: K3 one-loop factor 24 on units; T4 zero; both theta contractions vanish')
