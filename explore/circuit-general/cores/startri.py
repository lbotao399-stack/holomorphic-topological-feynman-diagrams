"""Exact star-triangle identity for projective periods at a unique vertex.

For a cubic vertex v with edges to x1,x2,x3 of powers a1,a2,a3 and a1+a2+a3 = q:
   P_G(a) = P_G'(a'),   G' = G - v + triangle (x2x3: q-a1, x1x3: q-a2, x1x2: q-a3).
(Triangle-to-star substitution x = S/r1, ...: U_tri = (x+y+z) U_star, unit log-Jacobian,
integrand factor S^(a_x+a_y+a_z-2q) = 1.)
"""


def star_to_triangle(edges, powers, v, q):
    inc = [i for i, e in enumerate(edges) if v in e]
    if len(inc) != 3:
        raise ValueError('vertex not cubic')
    if sum(powers[i] for i in inc) != q:
        raise ValueError('vertex not unique')
    nb = [edges[i][0] if edges[i][1] == v else edges[i][1] for i in inc]
    if len(set(nb)) != 3 or v in nb:
        raise ValueError('degenerate star')
    a = [powers[i] for i in inc]
    ne = [e for i, e in enumerate(edges) if i not in inc]
    npw = [p for i, p in enumerate(powers) if i not in inc]
    x1, x2, x3 = nb
    ne += [(x2, x3), (x1, x3), (x1, x2)]
    npw += [q - a[0], q - a[1], q - a[2]]
    return ne, npw


def unique_vertices(edges, powers, q):
    verts = sorted({x for e in edges for x in e})
    out = []
    for v in verts:
        inc = [i for i, e in enumerate(edges) if v in e]
        if len(inc) == 3 and sum(powers[i] for i in inc) == q:
            out.append(v)
    return out
