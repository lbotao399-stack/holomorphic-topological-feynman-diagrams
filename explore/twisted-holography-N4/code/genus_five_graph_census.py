"""Classify rational electrical charts on genus-five three-mark stable graphs.

The ordered charges are (1,0), (0,1), (-1,-1). Coincident marks have one
independent current and admit scalar rational integration. A two-terminal
series-parallel minor admits polynomial charts at arbitrary input currents.
"""
from collections import Counter
from functools import cache
from itertools import combinations

from betti_series_parallel import face_charts, vertices
from betti_star_charts import plan
from planar_four_loop import graphs, image, labelings, marked_shapes, sp

LEVELS = ((5,0,0),(4,1,0),(3,2,0),(3,1,1),(2,2,1))


def eligible(degrees,marks,size):
    """Minimum stable genus supporting the marked cotangent degrees at each vertex."""
    slots = [degree+int(j in marks) for j,degree in enumerate(degrees)]
    for channel,level in enumerate(LEVELS):
        required = sum(max(0,(4-b)//2,(5+(level[marks.index(i)] if i in marks else 0)-b)//3)
                       for i,b in enumerate(slots))
        if required<=5-(size-len(degrees)+1):
            yield channel


@cache
def safe_roots(edges):
    """Every nonempty inactive piece omits at least one endpoint of the active root."""
    good = {v for v in {v for e in edges for v in e} if sp(tuple(e for e in edges if v not in e))}
    return tuple(r for r,(v,w) in enumerate(edges) if v in good and w in good)


def nested_faces(graph,possible,root):
    """With polynomial residual loads, 'possible' bounds their nonzero vertex support."""
    if sp(tuple((v,w) for i,v,w in graph)): return ()
    bad = []
    for active,groups,pieces,collapsed,source,tree,others,sign in face_charts(graph,root):
        merged = {i for i,g in enumerate(groups) if possible.intersection(g)}
        if len(merged)<2 or (len(merged)==2 and len(groups)>2): continue
        empty = frozenset(range(len(groups)))-merged
        if tree is None and plan(collapsed,empty) is None:
            bad.append(len(active))
            continue
        # Residual loads are original loads minus incident active currents.
        residual = possible|vertices(active)
        for piece in pieces:
            if piece: bad.extend(nested_faces(piece,residual & vertices(piece),piece[0][0]))
    return tuple(sorted(set(bad)))


@cache
def tetrahedron_chart(edges,marked):
    """Allow one terminal charged tetrahedron in the positive electrical chart."""
    verts = {v for e in edges for v in e}
    if len(verts)<=2 or (len(verts)==4 and len(edges)==6): return True
    for v in sorted(verts):
        near = {w for e in edges if v in e for w in e if w!=v}
        if len(near)!=2 and (len(near)!=3 or v in marked): continue
        remaining = frozenset(e for e in edges if v not in e)|frozenset(combinations(sorted(near),2))
        loaded = marked-{v}
        if v in marked: loaded=loaded|near
        if tetrahedron_chart(remaining,loaded): return True
    return False


@cache
def marked_quotient(collapsed,groups,marks):
    """Canonical ordered marking on a contracted nonpolynomial chart."""
    edges=tuple(sorted(tuple(sorted((u,v))) for i,u,v in collapsed))
    marked=tuple(next(i for i,g in enumerate(groups) if v in g) for v in marks)
    return min((image(edges,p),tuple(p[k] for k in marked)) for p in labelings(edges))


def unresolved_faces(edges,marks,allow_tetrahedron=False,quotients=False):
    graph = tuple((i,*e) for i,e in enumerate(edges))
    roots = safe_roots(edges)
    if not roots:
        results = [nested_faces(graph,set(marks),r) for r in range(len(graph))]
        if len(edges)==10 and marks==(4,5,0):
            assert results==[()]*10
        return min(results,key=lambda x:(len(x),x))
    best = None
    for root in roots:
        bad = []
        for active,groups,pieces,collapsed,source,tree,others,sign in face_charts(graph,root):
            if tree is not None: continue
            occupancy = [sum(v in group for v in marks) for group in groups]
            # With two coincident marks all merged currents are collinear. For V>=3
            # the energy differential then has rank <= E-V+1 < E-1. V=2 is parallel.
            if occupancy.count(1)!=3 or len(active)-1>2*(len(active)-len(groups)+1): continue
            empty = frozenset(i for i,n in enumerate(occupancy) if not n)
            if plan(collapsed,empty) is None:
                simple=frozenset(tuple(sorted((u,v))) for i,u,v in collapsed)
                marked=frozenset(range(len(groups)))-empty
                if not allow_tetrahedron or not tetrahedron_chart(simple,marked):
                    bad.append(marked_quotient(collapsed,groups,marks) if quotients else len(active))
        if not bad: return ()
        bad = tuple(sorted(bad) if quotients else sorted(set(bad)))
        if best is None or (len(bad),bad)<(len(best),best):
            best = bad
        if best==(len(edges),): break
    return best


def census(size):
    counts = Counter()
    failures = []
    for edges in graphs(size):
        if sp(edges): continue
        counts['nonSP'] += 1
        n = max(max(e) for e in edges)+1
        degrees = tuple(sum(i in e for e in edges) for i in range(n))
        for nv,marks,aut in marked_shapes(edges):
            if len(set(marks))!=3 or any(d==1 and i not in marks for i,d in enumerate(degrees)):
                continue
            channels = tuple(eligible(degrees,marks,size))
            if not channels: continue
            counts['stable'] += 1
            bad = unresolved_faces(edges,marks)
            if bad:
                assert not unresolved_faces(edges,marks,allow_tetrahedron=True),(edges,marks,bad)
                failures.append((edges,marks,aut,channels,bad))
    counts['unresolved'] = len(failures)
    print(size,counts,tuple(sum(k in row[3] for row in failures) for k in range(len(LEVELS))),
          Counter(row[4] for row in failures),flush=True)
    return counts,failures


def compare_nine_edge_quotients(nine,ten):
    """Match every surviving ten-edge boundary to a nine-edge full-face marking."""
    expected=Counter()
    actual=Counter()
    for rows,histogram in ((nine,expected),(ten,actual)):
        for graph,marks,aut,channels,bad in rows:
            for q in unresolved_faces(graph,marks,quotients=True):
                histogram[q,channels]+=1
    assert len(expected)==6 and set(actual)==set(expected)
    assert sorted(expected.values())==[1]*6 and sorted(actual.values())==[1,1,2,2,3,3]
    assert len({q[0] for q,channels in expected})==2


def lower_edge_sums(last=5):
    """Exact rational graph sums through the selected edge count, all five level sets."""
    from fractions import Fraction
    from planar_four_loop import moment, stencil
    total = [Fraction() for _ in LEVELS]
    for size in range(last+1):
        rows = {}
        for (g,q,p,i),v in stencil(5,size,LEVELS).items():
            rows.setdefault((g,q,p),[Fraction() for _ in LEVELS])[i] += v
        for (g,q,p),weights in rows.items():
            value = moment(g,q,p)
            for i,w in enumerate(weights):
                total[i] += w*value
        print('Cumulative through',size,'inputs',len(rows),*total,flush=True)
    if last==5:
        assert total==list(map(Fraction,(
            '-48720516388020370063/4463304708120536678400000',
            '-18539964005616146181511/785541628629214455398400000',
            '-1024727409432589428751/26184720954307148513280000',
            '-931963601816152738099/39277081431460722769920000',
            '-2491511767896700529023/78554162862921445539840000')))
    return tuple(total)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edges',nargs='+',type=int,default=(6,7,8,9,10))
    parser.add_argument('--list',action='store_true',help='print unreduced marked graphs')
    parser.add_argument('--lower',action='store_true',help='also sum all moments through five edges')
    args = parser.parse_args()
    if args.lower: lower_edge_sums()
    residues={}
    for size in args.edges:
        counts,rows = census(size)
        residues[size]=rows
        if 6<=size<=10:
            assert (counts["nonSP"],counts["stable"])==((1,1),(3,33),(20,304),(118,743),(714,714))[size-6]
        if args.list:
            for row in rows: print(row,flush=True)

        if size in (9,10):
            expected=(0,3,0,3,0) if size==9 else (0,6,0,6,0)
            assert tuple(sum(i in row[3] for row in rows) for i in range(5))==expected
            assert all(row[-1]==(9,) for row in rows)
        if size in (6,7):
            assert not rows
        if size==8:
            graph = tuple(e for e in combinations(range(5),2) if e not in ((0,1),(0,2)))
            assert rows==[(graph,(0,1,2),2,(3,),(8,)),(graph,(1,2,0),2,(1,),(8,))]
    if 9 in residues and 10 in residues:
        compare_nine_edge_quotients(residues[9],residues[10])
