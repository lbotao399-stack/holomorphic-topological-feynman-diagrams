"""Two hard marks per bubble; two more together on one of the five other paths."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations, product
import argparse

from mixed_bare_double_chamber import (
    CORE, CYCLES, G, V, WAVES, bare_minor, fixture, forest_data, fmpz_mpoly_ctx,
    image, integrate, integrated_soft_placements, local_moment, numerator,
    plan_generic, Graph, solve_graph, symmetries, tree_flow,
)


def chambers():
    actions = symmetries()
    colors = sorted(set(permutations((0,0,1,1,2,2))))
    candidates = {l + (p,p) + r for l,r,p in product(
        combinations_with_replacement((0,1),2), combinations_with_replacement((7,8),2), range(2,7))}
    consumed = 0
    while candidates:
        marks = min(candidates)
        images = {tuple(sorted(g[e][0] for e in marks)) for g in actions}
        candidates -= images
        stabilizer = [g for g in actions if tuple(sorted(g[e][0] for e in marks)) == marks]
        assert len(stabilizer)*len(images) == len(actions)
        seen = set()
        rows = []
        for color in colors:
            if color in seen:
                continue
            config = tuple(tuple(color[j] for j,p in enumerate(marks) if p==e) for e in range(9))
            orbit = {sum(image(config,g),()) for g in stabilizer}
            seen.update(orbit)
            multiplicity = 8*len(images)*len(orbit)
            rows.append((color,multiplicity))
            consumed += multiplicity
        assert sum(h for _,h in rows) == 720*len(images)
        yield marks,rows
    assert consumed == 32400


def check_labeled_weights(all_shapes):
    """Recover all multiplicities by visiting every labeled chamber."""
    from fourfold_quartic_orbits import oriented_symmetries
    actions = tuple(oriented_symmetries(CORE))
    index = {}
    expected = {}
    for shape,(marks,rows) in enumerate(all_shapes):
        for color,count in rows:
            key = (shape,color)
            expected[key] = count
            config = tuple(tuple(color[i] for i,h in enumerate(marks) if h==e) for e in range(9))
            for g in actions:
                target = image(config,g)
                assert index.setdefault(target,key) == key
    actual = Counter()
    for left,right,path in product(combinations_with_replacement((0,1),2),
                                   combinations_with_replacement((7,8),2),range(2,7)):
        marks = left + (path,path) + right
        sites = [tuple(i for i,h in enumerate(marks) if h==e) for e in range(9)]
        for order in permutations(range(6)):
            color = tuple((0,1,1,2,2,0)[i] for i in order)
            config = tuple(tuple(color[i] for i in group) for group in sites)
            actual[index[config]] += 1
    assert actual == expected


def check_voltage(U,T,Z):
    """Compare the eight actual voltage rows as formal coefficients of each R_ij."""
    direct = Counter()
    for (a,b),(c,d),(f,g) in product(permutations(U),permutations(T),permutations(Z)):
        for mark,charge in ((b,-2),(f,1),(g,1)):
            direct[tuple(sorted((a,mark)))] -= charge
            direct[tuple(sorted((c,mark)))] += charge
    reference = Counter({tuple(sorted(U)):16})
    for left,right,sign in ((T,Z,1),(U,Z,-1),(T,U,-1)):
        for i,j in product(left,right):
            reference[tuple(sorted((i,j)))] += 4*sign
    # direct was multiplied by 2 to remove the 1/2 in the resistance formula.
    assert all(direct[k]==reference[k] for k in direct.keys()|reference.keys() if k[0]!=k[1])


def chamber_values(marks, colors, check_connectors=False, all_pairs=False):
    # The opposite bubble together with the central path and its two connectors
    # is a five-edge, two-cycle subgraph. At least four of the eight marks must
    # lie there. Only three do so when the two extra hard marks share a connector.
    if not all_pairs and marks[2] != 4 and not check_connectors:
        for color, multiplicity in colors:
            yield color, multiplicity, F()
        return
    paths = []
    for e in range(9):
        start = sum(map(len,paths))
        paths.append(tuple(range(start,start+marks.count(e)+1)))
    edges = fixture(CORE,marks)
    graph = G(12,edges)
    graph.C = [[c[e] for e in range(9) for _ in paths[e]] for c in CYCLES]
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    _,u,dist = forest_data(edges,12)
    resist = [[ctx.from_dict({m:int(c) for m,c in dist[6+i][6+j].terms()}) for j in range(6)] for i in range(6)]
    path_time = [sum(ctx.gens()[i] for i in path) for path in paths]
    groups = sum(paths[:4],()), paths[4], sum((paths[i] for i in (7,8,5,6)),())
    splits = ((len(paths[0]),len(paths[1])),(len(paths[7]),len(paths[8])))
    plan = plan_generic(Graph(12,tuple(edges)))
    lengths = list(range(1,16))
    tree = int(u(*lengths))
    for color,multiplicity in colors:
        U,T,Z = [tuple(i for i,c in enumerate(color) if c==k) for k in range(3)]
        check_voltage(U,T,Z)
        order = (U[0],*T,*Z,U[1])
        source = [(0,)*4]*12
        for j,wave in zip(order,WAVES[:6]):
            source[6+j] = wave
        sheared = [(a,b+c,d,c) for a,b,c,d in source]
        flow = tree_flow(graph,sheared)
        minors = {(i,j):numerator(graph,sheared,flow,[(3,i),(3,j)])
                  for i,j in (combinations(range(4),2) if all_pairs else product((0,1),(2,3)))}
        if not all_pairs and marks[2] != 4:
            assert not any(minors.values())
        def cross(x,y):
            return sum(resist[i][j] for i,j in product(x,y))
        mean4 = 4*resist[U[0]][U[1]] + cross(T,Z)-cross(U,Z)-cross(T,U)
        rays,_ = solve_graph(edges,source,lengths,0)
        squares = 0
        value = F()
        groups_of_squares = []
        for p,q in (combinations(range(9),2) if all_pairs else product((0,1),(7,8))):
            det = sum((CYCLES[i][p]*CYCLES[j][q]-CYCLES[j][p]*CYCLES[i][q])*poly
                      for (i,j),poly in minors.items())
            direct = bare_minor(rays,graph.C,lengths,[c[p] for c in CYCLES],[c[q] for c in CYCLES])
            assert abs(int(det(*lengths))) == abs(direct)*tree**4
            position = path_time[p]*path_time[q]
            squares += 2*int(position(*lengths))*direct**2
            if not det:
                continue
            content = int(det.content())
            if det.leading_coefficient() < 0:
                content = -content
            primitive = det//content
            for bucket in groups_of_squares:
                if bucket[0] == primitive:
                    bucket[1] += position*content**2
                    break
            else:
                groups_of_squares.append([primitive,position*content**2])
        for det,position in groups_of_squares:
            value += integrate(det,position*mean4,groups,splits,plan)
        assert squares == integrated_soft_placements(rays,graph.C,lengths,V)
        # Two soft label assignments; mean4 is four times the averaged voltage.
        yield color,multiplicity,value*F(multiplicity,2*len(symmetries()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list',action='store_true')
    parser.add_argument('--check-connectors',action='store_true')
    parser.add_argument('shapes',nargs='*',type=int)
    args = parser.parse_args()
    all_shapes = list(chambers())
    check_labeled_weights(all_shapes)
    if args.list:
        for j,(marks,rows) in enumerate(all_shapes):
            print(j,marks,len(rows),sum(h for _,h in rows))
        return
    subtotals = (F(), F(), -F(12917,226842638400000000), -F(258047,4990538044800000000),
                 F(), F(), -F(29399,3208203028800000000))
    counts = ((0,0,90),(0,0,54),(9,9,9),(8,9,13),(0,0,54),(0,0,33),(3,3,6))
    for j in args.shapes or range(len(all_shapes)):
        marks,colors = all_shapes[j]
        values = [v for _,_,v in chamber_values(marks,colors,args.check_connectors)]
        signs = sum(v>0 for v in values),sum(v<0 for v in values),sum(v==0 for v in values)
        assert sum(values) == subtotals[j] and signs == counts[j]
        local_moment.cache_clear()
        print(j,marks,signs,subtotals[j],flush=True)
    if not args.shapes:
        assert sum(subtotals) == -F(211663,1796593696128000000)
        print('Family subtotal:',sum(subtotals))


if __name__ == '__main__':
    main()
