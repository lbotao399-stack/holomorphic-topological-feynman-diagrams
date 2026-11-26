"""Mixed Newton moments by boundary integration, series-parallel division and empty-star elimination."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial,prod

from flint import fmpq,fmpq_mpoly_ctx,fmpz_mat


def components(vertices,edges):
    left=set(vertices)
    groups=[]
    while left:
        reached={min(left)}
        for _ in range(len(left)):
            for _,u,v in edges:
                if u in reached and v in left: reached.add(v)
                if v in reached and u in left: reached.add(u)
        left-=reached
        groups.append(tuple(sorted(reached)))
    return tuple(groups)


def vertices(edges):
    return {v for _,u,w in edges for v in (u,w)}


@lru_cache(maxsize=4096)
def network(edges,source,sink):
    """Two-terminal decomposition. Each leaf stores its original orientation."""
    if len(edges)==1:
        i,u,v=edges[0]
        assert {u,v}=={source,sink}
        return ('edge',i,1 if u==source else -1)
    all_vertices=vertices(edges)
    for middle in sorted(all_vertices-{source,sink}):
        parts=components(all_vertices-{middle},edges)
        if len(parts)!=2 or any(source in c and sink in c for c in parts): continue
        left=next(c for c in parts if source in c)
        first=tuple(e for e in edges if any(v in left for v in e[1:]))
        second=tuple(e for e in edges if e not in first)
        inside=tuple(sorted(vertices(first)-{source,middle}))
        return ('series',(*inside,middle),network(first,source,middle),network(second,middle,sink))
    groups=components(all_vertices-{source,sink},edges)
    branches=[tuple(e for e in edges if any(v in part for v in e[1:])) for part in groups]
    branches.extend((e,) for e in edges if {e[1],e[2]}=={source,sink})
    if len(branches)<2:
        raise ValueError('The contracted block has no two-terminal series-parallel decomposition.')
    return ('parallel',(),network(branches[0],source,sink),
            network(tuple(sorted(e for branch in branches[1:] for e in branch)),source,sink))


def determinant(rows,one):
    if not rows: return one
    if len(rows)==1: return rows[0][0]
    return sum(((-1)**i*x*determinant([r[:i]+r[i+1:] for r in rows[1:]],one)
                for i,x in enumerate(rows[0]) if x),one*0)


@lru_cache(maxsize=4096)
def cycle_minors(edges,root):
    """Signed Laplace terms of the energy gradient in two sets of chord currents."""
    tree={v:[] for v in vertices(edges)}
    part={v:v for v in tree}
    chords=[]
    for e,u,v in edges:
        if part[u]==part[v]:
            chords.append((e,u,v))
            continue
        left,right=part[u],part[v]
        part={k:(left if value==right else value) for k,value in part.items()}
        tree[u].append((v,e,1))
        tree[v].append((u,e,-1))
    basis=[]
    for chord,u,v in chords:
        # Close u->v by the unique tree path from v back to u.
        previous={v:None}
        pending=[v]
        while u not in previous:
            point=pending.pop()
            for neighbor,e,sign in tree[point]:
                if neighbor not in previous:
                    previous[neighbor]=(point,e,sign)
                    pending.append(neighbor)
        row={chord:1}
        point=u
        while point!=v:
            point,e,sign=previous[point]
            row[e]=sign
        basis.append(row)
    other=tuple(e for e,u,v in edges if e!=root)
    h=len(chords)
    n=len(other)
    assert 0<=2*h-n<=1
    terms={}
    omitted=(None,) if n==2*h else range(2*h)
    for slot in omitted:
        xrows=[j for j in range(h) if j!=slot]
        yrows=[j for j in range(h) if j+h!=slot]
        count=len(xrows)
        for ids in combinations(range(n),count):
            chosen=tuple(other[i] for i in ids)
            rest=tuple(e for e in other if e not in chosen)
            def minor(rows,columns):
                return int(fmpz_mat([[basis[j].get(e,0) for e in columns] for j in rows]).det())
            signed=(-1)**(sum(ids)-count*(count-1)//2)*minor(xrows,chosen)*minor(yrows,rest)
            if signed:terms.setdefault((chosen,rest),[]).append((0 if slot is None else slot,signed))
    return tuple(e for e,u,v in chords),terms


def cycle_jacobian(edges,root,current,axes,one):
    """Factor the energy differential through its chord currents by Cauchy–Binet."""
    chords,terms=cycle_minors(edges,root)
    dim=2*len(chords)
    saturated=len(axes)==dim
    projected=[one*0 for _ in range(1 if saturated else dim)]
    for (left,right),weights in terms.items():
        value=prod(current[e][0] for e in left)*prod(current[e][1] for e in right)
        for i,weight in weights:projected[i]+=weight*value
    if not any(projected):return one*0
    deriv=[[current[e][d].derivative(j) for j in axes] for d in (0,1) for e in chords]
    if saturated:return projected[0]*determinant(deriv,one)
    return sum((v*determinant(deriv[:i]+deriv[i+1:],one)
                for i,v in enumerate(projected) if v),one*0)


def barycenter_direction(tree,dimension):
    """Log-resistance derivatives at the center of the cube in preorder coordinates."""
    slopes={}
    def visit(node,vector,index):
        if node[0]=='edge':
            slopes[node[1]]=tuple(vector)
            return index
        sign=1 if node[0]=='series' else -1
        l=vector.copy()
        r=vector.copy()
        l[index]=sign
        r[index]=-sign
        return visit(node[3],r,visit(node[2],l,index+1))
    assert visit(tree,[0]*dimension,0)==dimension
    return slopes


def currents(tree,charges,parameter,injection,zero):
    rank=len(injection)
    def upward(node,index):
        if node[0]=='edge':
            return [zero]*rank,(node[1],node[2]),index
        var=parameter[index]
        a,first,nxt=upward(node[2],index+1)
        b,second,end=upward(node[3],nxt)
        if node[0]=='parallel':
            effective=[x+y for x,y in zip(a,b)]
            return effective,('parallel',var,a,effective,first,second),end
        load=[sum((charges[v][k] for v in node[1]),zero) for k in range(rank)]
        effective=[var*x+(1-var)*(y-z) for x,y,z in zip(a,b,load)]
        return effective,('series',load,first,second),end
    _,schedule,end=upward(tree,0)
    assert end==len(parameter)
    flows={}
    def downward(step,I):
        if isinstance(step[0],int):
            flows[step[0]]=[step[1]*x for x in I]
        elif step[0]=='series':
            _,load,first,second=step
            downward(first,I)
            downward(second,[x+y for x,y in zip(I,load)])
        else:
            _,t,a,b,first,second=step
            left=[x+t*(y-z) for x,y,z in zip(a,I,b)]
            downward(first,left)
            downward(second,[x-y for x,y in zip(I,left)])
    downward(schedule,injection)
    return flows


@lru_cache(maxsize=8192)
def face_charts(edges,root,proper=False):
    """Regular Newton faces meeting the chosen projected area coordinate."""
    all_vertices=vertices(edges)
    remaining=[i for i,u,v in edges if i!=root]
    charts=[]
    for count in range(len(remaining)+int(not proper)):
        for chosen in combinations(remaining,count):
            keep={root,*chosen}
            active=tuple(e for e in edges if e[0] in keep)
            inactive=tuple(e for e in edges if e[0] not in keep)
            groups=components(all_vertices,inactive)
            where={v:j for j,group in enumerate(groups) for v in group}
            collapsed=tuple((i,where[u],where[v]) for i,u,v in active)
            if any(u==v for _,u,v in collapsed): continue
            # Independent length scales on distinct blocks lower the rank of the area map.
            if any(len(components(set(range(len(groups)))-{j},collapsed))>1 for j in range(len(groups))):
                continue
            _,source,sink=next(e for e in collapsed if e[0]==root)
            nonroot=tuple(i for i,_,_ in active if i!=root)
            try:
                tree=network(collapsed,source,sink)
            except ValueError:
                tree=None
            sign=None
            if tree is not None:
                dimension=len(nonroot)
                slopes=barycenter_direction(tree,dimension)
                jac=fmpz_mat([[s-r for s,r in zip(slopes[i],slopes[root])] for i in nonroot]).det()
                assert jac
                sign=(-1)**(dimension+1)*(1 if jac>0 else -1)
            pieces=tuple(tuple(e for e in inactive if e[1] in group) for group in groups)
            charts.append((active,groups,pieces,collapsed,source,tree,nonroot,sign))
    return charts


def integrate_cube(poly,offset,dimension,ctx):
    # Each term of the zero-constant antiderivative vanishes at the lower
    # endpoint. Integrating one variable natively merges surviving monomials
    # before advancing to the next coordinate.
    for axis in range(offset,offset+dimension):
        poly=poly.integral(axis).subs({axis:1})
    return poly


def coefficient(edges,charges,powers,root=0,proper=False,*,omit_active=()):
    """Exact divided branch moment on series-parallel or source-free star reduction charts."""
    if len(edges)!=len(powers) or any(k<0 for k in powers) or not 0<=root<len(edges):
        raise ValueError('Supply a nonempty graph, a root edge and one nonnegative power per edge.')
    if any(map(sum,zip(*charges))): raise ValueError('Conserve every target coordinate.')
    if any(u==v for u,v in edges): return Fraction()
    graph=tuple((i,u,v) for i,(u,v) in enumerate(edges))
    if set(range(len(charges)))!=vertices(graph) or len(components(vertices(graph),graph))!=1:
        raise ValueError('Supply a connected graph with consecutive vertex labels.')
    ctx=fmpq_mpoly_ctx.get(tuple('x'+str(i) for i in range(max(len(edges)-1,1))))
    one=ctx.constant(1)
    zero=ctx.constant(0)
    variables=ctx.gens()
    def evaluate(subgraph,loads,offset,selected,proper_faces=False):
        total=zero
        for active,groups,pieces,collapsed,source,tree,other,sign in face_charts(subgraph,selected,proper_faces):
            if omit_active and subgraph==graph and tuple(e for e,u,v in active) in omit_active: continue
            nc=len(other)
            merged=[tuple(map(sum,zip(*(loads[v] for v in group)))) for group in groups]
            # With two opposite injections all currents have one common direction
            # at fixed outer parameters. Their E-V+1 cycle coordinates cannot
            # support an (E-1)-dimensional projection when V exceeds two.
            support=sum(any(load) for load in merged)
            if support<2 or (support==2 and len(groups)>2): continue
            if nc>len(merged[0])*(len(active)-len(groups)+1): continue
            parameter=variables[offset:offset+nc]
            if tree is None:
                from betti_star_charts import plan,solve,face_sign
                chart=plan(collapsed,frozenset(i for i,p in enumerate(merged) if not any(p)))
                if chart is None:
                    raise ValueError('No polynomial electrical reduction for this face and its charges.')
                current=solve(chart,dict(enumerate(merged)),parameter,zero)
            else:
                current=currents(tree,merged,parameter,merged[source],zero)
            energies={e:sum((x*x for x in t),zero)/2 for e,t in current.items()}
            axes=range(offset,offset+nc)
            if len(merged[0])==2 and nc>=6 and 0<=2*(len(active)-len(groups)+1)-nc<=1:
                integrand=cycle_jacobian(collapsed,selected,current,axes,one)
            else:
                integrand=determinant([[energies[e].derivative(j) for j in axes] for e in other],one)
            if not integrand or not energies[selected]: continue
            if sign is None: sign=face_sign(integrand)
            integrand*=sign*energies[selected]**(powers[selected]+1)/(powers[selected]+1)
            for e in other: integrand*=energies[e]**powers[e]
            if not integrand: continue
            residual={v:list(t) for v,t in loads.items()}
            for i,u,v in active:
                for j,x in enumerate(current[i]):
                    residual[u][j]-=x
                    residual[v][j]+=x
            for piece in pieces:
                if piece:
                    integrand*=evaluate(piece,{v:residual[v] for v in vertices(piece)},offset+nc,piece[0][0])
                    if not integrand: break
            total+=integrate_cube(integrand,offset,nc,ctx)
        return total
    inputs={v:[ctx.constant(fmpq(str(x))) for x in t] for v,t in enumerate(charges)}
    p=evaluate(graph,inputs,0,root,proper)
    assert all(not any(k) for k,c in p.terms())
    return Fraction(str(p))/prod(map(factorial,powers))


def check():
    from mixed_diamond_moment import graph_value
    from betti_power_graphs import shapes,vector_coefficient
    diamond=((0,3),(0,1),(1,3),(0,2),(2,3))
    charges=((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))
    for k in ((0,)*5,(1,0,0,0,0),(0,1,1,0,1)):
        wanted=graph_value(diamond,charges,k)
        assert all(coefficient(diamond,charges,k,i)==wanted for i in range(5))
    for e in range(1,5):
        for n,edges,*_ in shapes(e):
            p=[(i+1,i*i-2) for i in range(n-1)]
            p.append(tuple(-x for x in map(sum,zip(*p))))
            for k in ((0,)*e,(1,)+(0,)*(e-1)):
                assert coefficient(edges,p,k)==vector_coefficient(edges,tuple(p),k),(edges,p,k)
    print('All diamond root projections and degree-four independent formulas passed.',flush=True)
    theta=((0,1),(1,4),(0,2),(2,4),(0,3),(3,4))
    sources=((1,1,0),(0,-1,1),(2,-1,0),(0,1,1),(-3,0,-2))
    values=[coefficient(theta,sources,(0,)*6,i) for i in range(6)]
    assert len(set(values))==1
    print('Three two-edge branches:',values[0],flush=True)
    tetrahedron=tuple(combinations(range(4),2))
    triangular_loads=((1,0),(0,1),(-1,-1),(0,0))
    for powers,expected in (
        ((0,)*6,Fraction(-4657,217728000)),
        ((1,0,0,0,0,0),Fraction(-89987,58378320000)),
        ((0,0,2,0,0,0),Fraction(-3231731,22190837760000)),
    ):
        assert all(coefficient(tetrahedron,triangular_loads,powers,root)==expected
                   for root in (0,5))
    print('Empty-star tetrahedron and two decorated moments passed.',flush=True)


if __name__=='__main__':
    check()
