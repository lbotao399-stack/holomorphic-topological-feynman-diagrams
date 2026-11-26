"""Dolbeault top form of a (3,3,4) theta with one puncture contraction."""
from math import prod

from sympy import Matrix


# The first two coordinates are tangential, the last three are normal.
POSITIONS = [
    (-2, 1, -2, 1, -1), (-2, -1, -2, 2, -1), (-1, -1, -2, -2, -1),
    (-1, 2, 1, 2, 1), (-1, 2, 2, 0, 0), (-1, 2, 0, 1, 0),
    (1, 0, 2, -1, -1), (-1, 1, 2, 2, -2), (-2, -1, 0, -1, -1),
]
# Parent of each vertex in the eight-edge tree rooted at 0.
PARENT = [None, 0, 1, 2, 0, 3, 0, 8, 3]


def path_to_root(vertex):
    path = set()
    while PARENT[vertex] is not None:
        path.add(vertex)
        vertex = PARENT[vertex]
    return path


def check():
    points = list(map(Matrix, POSITIONS))
    tree = [(j, PARENT[j]) for j in range(1, 9)]
    paths = list(map(path_to_root, range(9)))
    tree_displacements = [points[i] - points[j] for i, j in tree]

    def pullback(target, source, covector):
        # Express a form on the point difference in the 40 tree coordinates.
        return [((j in paths[target]) - (j in paths[source])) * value
                for j in range(1, 9) for value in covector]

    # Special edge (4,5): tangent diagonal, target angular form, source cutoff.
    target, source = 4, 5
    assert points[target][:2] == points[source][:2]
    assert all(any(p[2:]) for p in points)
    assert points[target][2:] == [2, 0, 0] and points[source][2:] == [0, 1, 0]
    unit = Matrix.eye(5)
    special = [pullback(target, source, unit[:, j]) for j in range(2)]
    # y cross x is the covector sending v to d^3 bar w(y,x,v).
    angular = Matrix(points[source][2:]).cross(Matrix(points[target][2:]))
    special += [pullback(target, 0, [0, 0, *angular]),
                pullback(source, 0, [0, 0, *points[source][2:]])]

    # The other non-tree edge carries the ordinary flat angular form.
    other = (6, 7)
    other_disp = points[other[0]] - points[other[1]]
    extra = [pullback(*other, unit[:, j]) for j in range(5)]

    # First evaluate just on the one missing direction at each tree edge.
    reduced = []
    for idx, covector in enumerate(special + extra):
        reduced.append([Matrix(covector[5*k:5*k+5]).dot(d)
                        for k, d in enumerate(tree_displacements)]
                       + [0 if idx < 4 else other_disp[idx-4]])
    assert abs(Matrix(reduced).det()) == 1116

    # Independent 40-dimensional exterior product. For nonzero pivot u_k,
    # the wedge of (u_k e^j-u_j e^k) equals (-1)^k u_k^3 i_u(volume).
    full = special.copy()
    scales = []
    for i, j in [*tree, other]:
        displacement = points[i] - points[j]
        k = next(k for k, value in enumerate(displacement) if value)
        scales.append((-1)**k * displacement[k]**3)
        for direction in range(5):
            if direction != k:
                form = displacement[k]*unit[:, direction] - displacement[direction]*unit[:, k]
                full.append(pullback(i, j, form))
    coefficient = Matrix(full).det(method='domain-ge') / prod(scales)
    assert abs(coefficient) == 1116
    print('One puncture edge, nine flat edges: |angular coefficient| = 1116 (9x9 and 40x40).')

    # Apply the tangent holomorphic divergence to the special diagonal kernel.
    # The coefficients below are those of the leading derivatives of its delta.
    def missing(edge):
        i,j = edge
        return Matrix.hstack(*[((k in paths[i])-(k in paths[j]))*tree_displacements[k-1]
                               for k in range(1,9)])

    blocks = missing((target,source)).col_join(missing(other))
    tail = Matrix.zeros(5,1)
    def odd_coefficient(direction):
        return blocks.row_join(direction.col_join(tail)).row_join(tail.col_join(other_disp)).det()

    assert [odd_coefficient(unit[:,k]) for k in range(5)] == [1584,-1026,558,1116,-1674]
    # Replacing the puncture direction by the flat difference restores the null line.
    assert odd_coefficient(points[target]-points[source]) == 0

    # Include all five common antiholomorphic directions. Put the flux at 1
    # and a test with anti-form dzbar1 dzbar2 dwbar1 and normal trivector at 2.
    def absolute_row(vertex, covector):
        tree_part = [int(j in paths[vertex])*covector.dot(tree_displacements[j-1])
                     for j in range(1,9)]
        return Matrix([*covector,*tree_part,0]).T

    def relative_row(i,j,covector):
        return absolute_row(i,covector)-absolute_row(j,covector)

    rows = [relative_row(target,source,unit[:,k]) for k in range(2)]
    rows += [absolute_row(target,Matrix([0,0,*angular])),
             absolute_row(source,Matrix([0,0,*points[source][2:]]))]
    for k in range(5):
        row = relative_row(*other,unit[:,k])
        row[-1] = other_disp[k]
        rows.append(row)
    normal = Matrix(points[1][2:])
    rows += [absolute_row(1,normal[0]*unit[:,k+2]-normal[k]*unit[:,2]) for k in (1,2)]
    rows += [absolute_row(2,unit[:,k]) for k in (0,1,2)]
    total_anti = Matrix.vstack(*rows).det(method='domain-ge') / normal[0]
    assert total_anti == -2232 and abs(total_anti*odd_coefficient(unit[:,0])) == 3535488
    print('One flux and normal trivector test: |first tangent delta-derivative coefficient| = 3535488.')


if __name__ == '__main__':
    check()
