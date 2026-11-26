"""One six-hard-mark chamber of the genus-four mixed Ward projection."""

from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod

from mixed_bare_ranks import WAVES
from quartic_seven_plane_double import (
    Arm, CORE, CYCLES, G, Q, SCALE, decompose, fixture, forest_data,
    fmpz_mpoly_ctx, image, local_moment, numerator, symmetries, tree_flow,
)
from fourfold_mixed_soft_ward import V, bare_minor, integrated_soft_placements
from point_graph import solve_graph
from square_three import Graph, plan_generic, sp_monomial


MARKS = (0, 0, 4, 4, 7, 7)


def integrate(square_root, position_voltage, groups, splits, plan):
    left, right = (Arm(group, split) for group, split in zip((groups[0], groups[2]), splits))

    def factors(poly):
        return [(left.register(l), right.register(r), m, c)
                for l, r, m, c in decompose(poly, groups)]

    minors = factors(square_root)
    weights = factors(position_voltage)
    total = F()
    tested = False
    for i, j in combinations_with_replacement(range(len(minors)), 2):
        x, y = minors[i], minors[j]
        for z in weights:
            a, b, c = (h[0] for h in (x, y, z))
            d, e, f = (h[1] for h in (x, y, z))
            middle = tuple(sum(col) for col in zip(x[2], y[2], z[2]))
            sl = sum(left.degrees[h] for h in (a, b, c)) + len(left.group) - Q
            sr = sum(right.degrees[h] for h in (d, e, f)) + len(right.group) - Q
            sm = sum(middle) + len(middle)
            assert sl + sr + sm == 2 * Q and 0 < min(sl, sr, sm) <= max(sl, sr, sm) < Q
            outer = F(prod(factorial(Q-1-s) for s in (sl, sm, sr)), factorial(Q-1))
            mid = F(prod(map(factorial, middle)), factorial(sm-1))
            if not tested:
                powers = [0] * 15
                check = outer * mid
                for arm, ids in ((left, (a, b, c)), (right, (d, e, f))):
                    ex = tuple(map(sum, zip(*(next(iter(arm.polys[h].monoms())) for h in ids))))
                    check *= F(local_moment(arm.split, ex), SCALE)
                    for k, power in zip(arm.group, ex):
                        powers[k] = power
                for k, power in zip(groups[1], middle):
                    powers[k] = power
                assert check == sp_monomial(powers, Q, plan)
                tested = True
            total += ((1 + (i != j)) * x[3] * y[3] * z[3] * outer * mid
                      * F(left.integrate(a, b, c) * right.integrate(d, e, f), SCALE**2))
    return total


def chamber():
    paths = []
    for e in range(9):
        start = sum(map(len, paths))
        paths.append(tuple(range(start, start + MARKS.count(e) + 1)))
    edges = fixture(CORE, MARKS)
    graph = G(12, edges)
    graph.C = [[c[e] for e in range(9) for _ in paths[e]] for c in CYCLES]
    source = [(0,)*4]*6 + list(WAVES[:6])
    # This unimodular map sends the common soft direction to the last unit vector.
    sheared = [(a, b+c, d, c) for a,b,c,d in source]
    flow = tree_flow(graph, sheared)
    minors = {(i,j): numerator(graph, sheared, flow, [(3,i), (3,j)])
              for i,j in product((0,1), (2,3))}
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    times = ctx.gens()
    path_time = [sum(times[i] for i in path) for path in paths]
    _, u, resistances = forest_data(edges, 12)
    voltage = -sum(source[j][3]*(resistances[6][j]-resistances[7][j]) for j in range(12))/2

    def cross(left, right):
        return sum(resistances[6+i][6+j] for i,j in product(left, right))

    U,T,Z = (0,5),(1,2),(3,4)
    mean4 = 4*resistances[6][11] + cross(T,Z)-cross(U,Z)-cross(T,U)

    def integral_polynomial(poly):
        assert all(c.denominator == 1 for _,c in poly.terms())
        return ctx.from_dict({m:int(c) for m,c in poly.terms()})

    voltage, mean4 = map(integral_polynomial, (voltage, mean4))
    groups = sum(paths[:4], ()), paths[4], sum((paths[i] for i in (7,8,5,6)), ())
    splits = ((len(paths[0]), len(paths[1])), (len(paths[7]), len(paths[8])))
    plan = plan_generic(Graph(12, tuple(edges)))
    lengths = list(range(1,16))
    tree = int(u(*lengths))
    rays, phi = solve_graph(edges, source, lengths, 0)
    assert int(voltage(*lengths)) == tree*(phi[6][3]-phi[7][3])
    check_voltage = 0
    for (a,b),(c,d),(f,g) in product(permutations(U), permutations(T), permutations(Z)):
        assigned = source.copy()
        for j,wave in zip((a,c,d,f,g,b), WAVES[:6]):
            assigned[6+j] = wave
        _, potential = solve_graph(edges, assigned, lengths, 0)
        check_voltage += potential[6+a][3]-potential[6+c][3]
    assert 2*int(mean4(*lengths)) == tree*check_voltage
    pairs = tuple(product((0,1), (7,8)))
    raw, averaged = [], []
    squares = 0
    for p,q in pairs:
        det = sum(CYCLES[i][p]*CYCLES[j][q]*poly for (i,j),poly in minors.items())
        direct = bare_minor(rays, graph.C, lengths, [c[p] for c in CYCLES], [c[q] for c in CYCLES])
        assert abs(int(det(*lengths))) == abs(direct)*tree**4
        position = path_time[p]*path_time[q]
        squares += 2*int(position(*lengths))*direct**2
        raw.append(integrate(det, position*voltage, groups, splits, plan))
        averaged.append(integrate(det, position*mean4, groups, splits, plan)/4)
    assert squares == integrated_soft_placements(rays, graph.C, lengths, V)
    # The color word at the sites is U,T,T,Z,Z,U. Count the paths and
    # colored orders in its orbit, then restore the 2! choices in each pair.
    colors = (0,1,1,2,2,0)
    config = tuple(tuple(colors[j] for j,h in enumerate(MARKS) if h == e) for e in range(9))
    orbit = {image(config,g) for g in symmetries()}
    weight = F(8*len(orbit),len(symmetries()))
    return raw, 2*sum(raw), averaged, 2*weight*sum(averaged), weight


if __name__ == '__main__':
    raw, hard, averaged, weighted, weight = chamber()
    assert raw == [F(2741*k, 24638999261184000000) for k in (4,2,2,1)]
    assert hard == F(2741, 1368833292288000000)
    assert averaged == [F(56601521,301827740949504000000000), F(35073791,603655481899008000000000),
                        F(18537619,301827740949504000000000), F(3886877,301827740949504000000000)]
    assert weighted == F(858337, 167682078305280000000) and weight == 8
    print('Hard order with both soft assignments:', hard)
    print('Complete colored chamber orbit:', weighted)
