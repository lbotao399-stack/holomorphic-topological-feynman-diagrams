"""All-loop fivefold marked scalar kernel on three repeated momentum directions."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial,prod
from flint import fmpq_mat,fmpz_mat
from point_graph import solve_graph


def fan(n):
    edges=[(j,0) for j in range(1,n+2)]+[(j,j+1) for j in range(1,n+1)]
    cycles=[[-int(e==j)+int(e==j+1)+int(e==n+1+j) for e in range(2*n+1)] for j in range(n)]
    return edges,cycles


def trees(cycles):
    return [s for s in combinations(range(len(cycles[0])),len(cycles))
            if fmpz_mat([[c[j] for j in s] for c in cycles]).det()]


def cofactor(y,cycles,t):
    rows=[[1]*len(t)]+[[str(c[j]*w[i]) for j,w in enumerate(y)] for i in range(len(y[0])) for c in cycles]
    return abs(F(str(fmpq_mat(rows).det()))/sum(t))


def chart(l,b):
    n=len(l)
    spokes=[None]*(n+1)
    path=[None]*n
    r=F(1)
    jacobian=F(1)
    for i in reversed(range(n)):
        spokes[i+1]=r/(1-b[i])
        path[i]=(1-l[i])*r/b[i]
        jacobian*=r**2/(b[i]**3*(1-b[i])**2)
        r*=l[i]/b[i]
    spokes[0]=r
    return spokes+path,jacobian


def witness(n, cotree_colors=3):
    edges,cycles=fan(n)
    l=[F(i+1,2*i+3) for i in range(n)]
    b=[F(i+2,3*i+5) for i in range(n)]
    t,jac=chart(l,b)
    src=[(-1,-1),(1,0)]+[(0,0)]*(n-1)+[(0,1)]
    y,phi=solve_graph(edges,src,t,0)
    c=cofactor(y,cycles,t)
    cotrees=trees(cycles)
    U=sum(prod(t[j] for j in s) for s in cotrees)
    weight=prod(l[i]**(2*(n-1-i))*b[i]**(2*i+1)*(1-b[i])**2 for i in range(n))
    assert c*c*jac/U**2==weight
    green=[F()]+[prod(l[j:]) for j in range(n+1)]
    assert all(phi[j][1]==green[j] for j in range(1,n+2))
    tree_mean=sum(prod(t[j] for j in s)*sum(green[edges[j][0]]+green[edges[j][1]] for j in s)
                  for s in cotrees)/U
    assert tree_mean==sum(green[2:])

    # Restore each suppressed momentum direction with one source per subdivision.
    for offset in range(3):
        cuts=[cotrees[(offset+k*(len(cotrees)//cotree_colors))%len(cotrees)] for k in range(cotree_colors)]
        new_src=[(0,)*cotree_colors+p for p in src]
        new_edges=[]
        new_times=[]
        parent=[]
        for e,(i,j) in enumerate(edges):
            points=[i]
            for k,s in enumerate(cuts):
                if e in s:
                    points.append(len(new_src))
                    new_src.append(tuple(int(k==h) for h in range(cotree_colors+2)))
            points.append(j)
            nseg=len(points)-1
            new_edges.extend(zip(points,points[1:]))
            new_times.extend([t[e]/nseg]*nseg)
            parent.extend([e]*nseg)
        new_src[0]=tuple(-sum(col) for col in zip(*new_src[1:]))
        new_y,new_phi=solve_graph(new_edges,new_src,new_times,0)
        new_cycles=[[r[e] for e in parent] for r in cycles]
        assert cofactor(new_y,new_cycles,new_times)==c
        potential=green[1]+sum(F(green[edges[e][0]]+green[edges[e][1]],2) for s in cuts for e in s)
        assert sum(new_phi[n+1][:cotree_colors+1])==potential
    return weight


def closed(n):
    return F((3*n+1)*factorial(n),2**n*(n+1)*factorial(2*n+1))


def independent_integral(n, cotree_colors=3):
    # Evaluate the graph covariance directly on a quadrature of the beta product.
    # After removing its density the cotree mean is affine in EACH split. Its
    # expected value can therefore be evaluated at the independent beta means.
    l=[F(2*(n-i)-1,2*(n-i)) for i in range(n)]
    b=[F(2*i+2,2*i+5) for i in range(n)]
    t,_=chart(l,b)
    edges,cycles=fan(n)
    src=[(-1,-1),(1,0)]+[(0,0)]*(n-1)+[(0,1)]
    _,phi=solve_graph(edges,src,t,0)
    potential=[F()]+[phi[j][1] for j in range(1,n+2)]
    cuts=trees(cycles)
    u=sum(prod(t[j] for j in s) for s in cuts)
    prob=[sum(prod(t[j] for j in s) for s in cuts if e in s)/u for e in range(len(edges))]
    observable=potential[1]+F(cotree_colors,2)*sum((potential[i]+potential[j])*p for (i,j),p in zip(edges,prob))
    mass=prod(F(2,(2*(n-i)-1)*(2*i+2)*(2*i+3)*(2*i+4)) for i in range(n))
    return factorial(n)**cotree_colors*mass*observable


if __name__=='__main__':
    for n in range(1,7):
        assert witness(n)>0
        value=independent_integral(n)
        assert value==closed(n)
        print(n,value)
    assert closed(1)==F(1,6) and closed(2)==F(7,720)
