"""Independent graph fixtures and all-six irrationality witness."""
from fractions import Fraction as F
from kernel import K4SUB,evaluate

SOURCE = [(1,0,0),(0,1,0),(0,0,1),(1,1,-1),(1,-1,1),(-1,1,1),(-1,-1,0),(-1,-1,-2)]
EXPECTED = [(-20383,1209600,83,4200),(1333,57600,-23,1200),(-7495,217728,5,168),
            (394519,10886400,-307,8400),(8663,259200,-29,1200),(-409,226800,11,8400)]


if __name__=='__main__':
    for (name,edges),(a,b,c,d) in zip(K4SUB,EXPECTED):
        assert evaluate(edges,SOURCE)==(F(a,b),F(c,d))
    edges = [(0,1),(0,2),(0,3),(1,4),(4,2),(1,3),(2,3)]
    source = [(1,2),(2,-1),(-1,3),(-3,-4),(1,0)]
    expected = [(F(-71,12),0),(F(222317,25200),-F(3021,700)),(-F(4003781,302400),F(3772,525))]
    for m,value in enumerate(expected):
        assert evaluate(edges,source,[(0,0)]*m)==value
    print('Six eight-source kernels and three independent seven-edge moments agree.')
    from htk.catalog import is_series_parallel,tight_graphs
    from htk.spred import sp_reduce
    catalog = tight_graphs(3,3)
    checked = 0
    for g in catalog:
        if is_series_parallel(g):
            result = evaluate(g.edges,SOURCE)
            exact = sp_reduce(g,3).integrate(SOURCE[:-1])
            assert result==(exact,0)
            checked += 1
    assert (checked,len(catalog))==(10,16)
    print('Ten series-parallel kernels agree with the independent ratio-coordinate engine.')
