"""Direct resistance-to-ratio Jacobian, oriented density and current checks."""
from fractions import Fraction as F
from flint import fmpq, fmpq_mat, fmpq_series
from forest_jets import chart, edge, moves
from forest_orders import forests,graph
from vk import G,currents_exact,density_exact


def ratios(edges,a,b,times):
    t=dict(zip(edges,times))
    result=[]
    def star(v,i,j,k):
        x,y,z=(t.pop(edge(v,w)) for w in (i,j,k))
        s=x*y+y*z+z*x
        return s/z,s/y,s/x
    for kind,vertices,_ in moves(edges,a,b):
        if kind=='series':
            v,u,w=vertices
            l,r=t.pop(edge(v,u)),t.pop(edge(v,w))
            total=l+r
            old=t[edge(u,w)]
            share=total/(total+old)
            result.extend((r/total,share))
            t[edge(u,w)]=old*share
        elif kind=='single':
            v,u,w,x=vertices
            direct,left,right=star(v,u,w,x)
            t[edge(u,w)]=direct
            for key,branch in ((edge(u,x),left),(edge(w,x),right)):
                share=branch/(branch+t[key])
                result.append(share)
                t[key]*=share
        else:
            p,q,u,a,b=vertices
            first=star(p,u,a,b)
            second=star(q,u,a,b)
            for j,v in enumerate((a,b)):
                share=second[j]/(first[j]+second[j])
                t[edge(u,v)]=first[j]*share
                result.append(share)
            old=t[edge(a,b)]
            new=1/(1/old+1/first[2]+1/second[2])
            result.extend((new/first[2],new/second[2]))
            t[edge(a,b)]=new
    return result


def check(forest,offset):
    edges,a,b,roots,_=graph(forest)
    source=[(0,0)]*a+[(1,0),(0,0)]
    for j,v in enumerate(roots):
        source[v]=(j%3,j+1)
    source[b]=tuple(-sum(p[j] for p in source) for j in (0,1))
    d,flows,_=chart(edges,a,b,source)
    times=[1]+[1+(i+offset)**3 for i in range(1,len(edges))]
    # Differentiate on the t_0=1 projective section in the dual numbers.
    jac=[]
    for j in range(1,len(edges)):
        t=[fmpq_series([x,int(i==j)],prec=2) for i,x in enumerate(times)]
        r=ratios(edges,a,b,t)
        jac.append([v[1] for v in r])
    points=[v[0] for v in r]
    det=abs(fmpq_mat(jac).det())
    assert det>0
    g=G(len(source),edges)
    currents=currents_exact(g,list(map(F,times)),source[:-1])
    for e,ref in zip(edges,currents):
        assert all(F(str(v(*points)))==wanted for v,wanted in zip(flows[e],ref))
    assert F(str(d(*points)*det))==density_exact(g,list(map(F,times)),source[:-1])


if __name__=='__main__':
    for n in range(1,6):
        for f in forests(n):
            for offset in (1,3):
                check(f,offset)
        print('loop',n,'currents, full Jacobian and density checked',flush=True)
