"""Group all remaining genus-five nine-active-edge integrals before tetrahedron reduction."""
from math import prod

import sympy as sp

from genus_five_graph_census import (LEVELS, census, face_charts, image, labelings,
                                    marked_quotient, plan, safe_roots)
from planar_four_loop import stencil


def loads(marks,n):
    return tuple((int(v==marks[0])-int(v==marks[2]),int(v==marks[1])-int(v==marks[2]))
                 for v in range(n))


def cycle_currents(edges,marks):
    """Use the non-tree edge currents themselves as the eight coordinates."""
    n=max(max(e) for e in edges)+1
    b=sp.Matrix([[int(v==i)-int(v==j) for i,j in edges] for v in range(n-1)])
    tree=tuple(b.rref()[1])
    chords=tuple(j for j in range(len(edges)) if j not in tree)
    assert len(chords)==4
    variables=sp.symbols('x0:4 y0:4')
    free=sp.Matrix([[variables[j],variables[4+j]] for j in range(4)])
    fixed=b[:,tree].inv()*(sp.Matrix(loads(marks,n)[:-1])-b[:,chords]*free)
    currents=sp.zeros(len(edges),2)
    for ids,block in ((tree,fixed),(chords,free)):
        for i,e in enumerate(ids): currents[e,:]=block[i,:]
    assert b*currents==sp.Matrix(loads(marks,n)[:-1])
    energy=tuple(currents.row(e).dot(currents.row(e))/2 for e in range(len(edges)))
    return currents,energy,variables


def active_face(edges,marks,root):
    graph=tuple((i,*e) for i,e in enumerate(edges))
    failures=[]
    for active,groups,pieces,collapsed,source,tree,other,sign in face_charts(graph,root):
        if tree or len(active)-1>2*(len(active)-len(groups)+1): continue
        occupied=tuple(next(i for i,g in enumerate(groups) if v in g) for v in marks)
        if len(set(occupied))<3 or plan(collapsed,frozenset(range(len(groups)))-set(occupied)): continue
        key=marked_quotient(collapsed,groups,marks)
        missing,=set(range(10))-{e for e,u,v in active}
        failures.append((key,missing,collapsed,occupied))
    assert len(failures)==1
    return failures[0]


def residual_insertions(channel,nine,ten):
    """Return quadratic insertions multiplying positive E_root times projected area."""
    targets=[(g,m) for g,m,aut,channels,bad in nine if channel in channels]
    stars=[(g,m) for g,m,aut,channels,bad in ten if channel in channels]
    assert len(targets)==3 and len(stars)==6
    root_by_key=dict(zip(targets,(3,3,2)))
    slots={}
    for key,root in root_by_key.items():
        edges,marks=key
        assert root in safe_roots(edges)
        current,energy,variables=cycle_currents(edges,marks)
        total=0
        ref=loads(marks,6)
        for (g,q,p,i),w in stencil(5,9,(LEVELS[channel],),shape_subset=(edges,)).items():
            if q!=ref:continue
            assert sum(p)<=1
            total-=w*prod(v**k for v,k in zip(energy,p))/(1+p[root])
        slots[key]=[current,variables,sp.expand(total)]
    for (edges,marks),r in zip(stars,(3,3,2,3,2,3)):
        assert r in safe_roots(edges)
        key,missing,collapsed,occupied=active_face(edges,marks,r)
        target,source=key
        currents,variables,total=slots[key]
        es=tuple(sorted(tuple(sorted((u,v))) for e,u,v in collapsed))
        permutations=[p for p in labelings(es) if (image(es,p),tuple(p[v] for v in occupied))==key]
        u,v=next((u,v) for e,u,v in collapsed if e==r)
        perm=next(p for p in permutations if tuple(sorted((p[u],p[v])))==target[root_by_key[key]])
        ref=loads(marks,7)
        weight=sum(w for (g,q,p,i),w in stencil(5,10,(LEVELS[channel],),shape_subset=(edges,)).items() if q==ref)
        residual=[]
        for vertex in edges[missing]:
            load=sp.Matrix([ref[vertex]])
            for e,u,v in collapsed:
                if vertex not in edges[e]:continue
                j=target.index(tuple(sorted((perm[u],perm[v]))))
                sign=(1 if perm[u]<perm[v] else -1)*(1 if edges[e][0]==vertex else -1)
                load-=sign*currents.row(j)
            residual.append(load)
        assert residual[0]+residual[1]==sp.zeros(1,2)
        slots[key][2]+=weight*residual[0].dot(residual[0])/2
    return {key:(root_by_key[key],sp.Poly(value,variables,domain=sp.QQ))
            for key,(current,variables,value) in slots.items()}


def check_insertions(channel,nine,ten):
    for index,((edges,marks),(root,poly)) in enumerate(residual_insertions(channel,nine,ten).items()):
        xy=poly.gens
        expected=0
        for j0,j1,j2,j3 in (xy[:4],xy[4:]):
            if index==0: expected+=(j2+j3)**2
            elif index==1: expected+=(j2+j3)**2+2*(j0+j2)**2
            else: expected+=(2*j0+j2+j3)**2+2*j2*j2+2*j3*j3
        assert sp.expand(96*poly.as_expr()+expected)==0
        # Grounded Kirchhoff currents certify an open set of strict negativity.
        b=sp.Matrix([[int(v==i)-int(v==j) for i,j in edges] for v in range(5)])
        resistances=[j+2 for j in range(9)]
        resistances[root]=1
        c=sp.diag(*(sp.Rational(1,r) for r in resistances))
        transfer=c*b.T*(b*c*b.T).inv()
        currents=transfer*sp.Matrix(loads(marks,6)[:-1])
        projection=c-transfer*b*c
        hessian=-projection.multiply_elementwise(currents*currents.T)
        assert hessian.minor_submatrix(root,root).det()>0 and currents.row(root).dot(currents.row(root))>0
        chords=[j for j in range(9) if j not in b.rref()[1]]
        point=[currents[e,k] for k in range(2) for e in chords]
        assert poly.eval(dict(zip(xy,point)))<0
        print(LEVELS[channel],marks,'root edge',edges[root],'negative insertion',index+1,flush=True)


if __name__=='__main__':
    nine=census(9)[1]
    ten=census(10)[1]
    for channel in (1,3):
        check_insertions(channel,nine,ten)
