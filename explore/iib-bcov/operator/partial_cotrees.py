"""Factor the joint path average of two-source colors at three loops into second moments."""
from fractions import Fraction as F
from itertools import permutations
from math import prod
from pathlib import Path
import sys

from flint import fmpq_mat

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from point_graph import solve_graph
from scalar_fan import cofactor
from square_three import G


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def cross(x,y):
    return tuple(x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3))


def sign(perm):
    return (-1)**sum(b<a for i,a in enumerate(perm) for b in perm[i+1:])


def mixed(moments):
    c=len(moments)
    perms=[(p,sign(p)) for p in permutations(range(c))]
    return sum(sa*sb*prod(m[a[k]][b[k]] for k,m in enumerate(moments))
               for a,sa in perms for b,sb in perms)


class PartialCotrees:
    """The first c coordinates each have two deleted insertions. All others are retained."""
    def __init__(self,edges,source,times,pairs):
        self.edges=tuple(edges)
        self.times=tuple(map(F,times))
        self.pairs=tuple(pairs)
        self.source=tuple(source)
        c=len(pairs)
        self.nv=len(source)
        if len(edges)!=16-2*c or self.nv!=len(edges)-2 or not 1<=c<=4:
            raise ValueError('A reduced tight fivefold three-loop graph is required.')
        if tuple(map(sum,zip(*source)))!=tuple([-sum(g) for g in pairs]+[0]*(5-c)):
            raise ValueError('Retained charges and deleted pairs must balance.')
        graph=G(self.nv,list(edges))
        self.cycle=tuple(zip(*graph.C))
        # Columns contain unit currents from a vertex to ground (vertex zero).
        loads=[[int(i==v)-int(v==0) for i in range(1,self.nv)] for v in range(self.nv)]
        flow,phi=solve_graph(edges,loads,self.times,0)
        self.flow=tuple((F(0),*row) for row in flow)
        self.green=tuple((F(0),*phi[v]) for v in range(self.nv))
        self.retained=tuple(tuple(dot(y,col) for y in self.flow) for col in zip(*source))

        # Fix one edge column to remove the projective relation, and pivot the
        # rows in the remaining 5-c colors before contracting the c random rows.
        fixed=[[str(ch[j]*y) for ch,y in zip(self.cycle[:-1],self.retained[i][:-1])]
               for i in range(c,5) for j in range(3)]
        matrix=fmpq_mat(fixed)
        echelon,rank=matrix.rref()
        self.factor=F(0)
        if rank!=len(fixed):
            return
        self.pivots=[next(j for j in range(matrix.ncols()) if echelon[i,j])
                     for i in range(rank)]
        self.free=[j for j in range(matrix.ncols()) if j not in self.pivots]
        self.projection=[[(j,F(str(echelon[i,h]))) for i,j in enumerate(self.pivots)
                          if echelon[i,h]] for h in self.free]
        square=fmpq_mat([[matrix[i,j] for j in self.pivots] for i in range(rank)])
        self.factor=F(str(square.det()))**2/self.times[-1]**2

    def projected(self,k,e,f,s,t):
        """Residual row after its two insertions and the fixed-color pivots."""
        a,b=self.pairs[k]
        null=cross(self.cycle[e],self.cycle[f])
        u,v=self.edges[e]; w,x=self.edges[f]
        row=[]
        for i,(y,cycle) in enumerate(zip(self.flow[:-1],self.cycle[:-1])):
            current=(self.retained[k][i]+a*((1-s)*y[u]+s*y[v])
                     +b*((1-t)*y[w]+t*y[x]))
            row.append(a*b*dot(null,cycle)*current)
        return tuple(row[h]-sum(row[j]*v for j,v in column)
                     for h,column in zip(self.free,self.projection))

    def voltage(self,load):
        if any(load[:len(self.pairs)]):
            raise ValueError('The projected voltage load must be supported on retained directions.')
        q=[dot(p,load) for p in self.source]
        return tuple(dot(g,q) for g in self.green)

    def mass(self,load=None,target=None):
        """Sum all labeled placements and integrate their positions at the supplied times.

        With a load, target is a retained vertex number or (deleted color, position in pair).
        The output uses potential relative to ground. Differences remove the ground choice.
        """
        if not self.factor:
            return F(0)
        if (load is None)!=(target is None):
            raise ValueError('Supply both load and target for a voltage.')
        potential=self.voltage(load) if load is not None else ()
        moments=[]
        c=len(self.pairs)
        for k,(a,b) in enumerate(self.pairs):
            mat=[[F(0) for _ in range(c)] for _ in range(c)]
            cache={}
            for e,f in permutations(range(len(self.edges)),2):
                null=cross(self.cycle[e],self.cycle[f])
                if not any(null):
                    continue
                if null not in cache:
                    weights=[dot(null,col) for col in self.cycle[:-1]]
                    def project(values):
                        row=[w*y for w,y in zip(weights,values)]
                        return tuple(row[h]-sum(row[j]*v for j,v in column)
                                     for h,column in zip(self.free,self.projection))
                    cache[null]=(project(self.retained[k]),[project(v) for v in zip(*self.flow)])
                fixed,response=cache[null]
                u,v=self.edges[e]
                w,x=self.edges[f]
                mu=[a*b*(fixed[j]+a*(response[u][j]+response[v][j])/2
                          +b*(response[w][j]+response[x][j])/2) for j in range(c)]
                ds=[a*a*b*(response[v][j]-response[u][j]) for j in range(c)]
                dt=[a*b*b*(response[x][j]-response[w][j]) for j in range(c)]
                avg,gs,gt=F(1),F(0),F(0)
                if isinstance(target,tuple) and target[0]==k:
                    if target[1]==0:
                        avg=(potential[u]+potential[v])/2
                        gs=potential[v]-potential[u]
                    else:
                        avg=(potential[w]+potential[x])/2
                        gt=potential[x]-potential[w]
                vol=self.times[e]*self.times[f]
                for i in range(c):
                    for j in range(i,c):
                        second=mu[i]*mu[j]+(ds[i]*ds[j]+dt[i]*dt[j])/12
                        weighted=(gs*(ds[i]*mu[j]+mu[i]*ds[j])
                                  +gt*(dt[i]*mu[j]+mu[i]*dt[j]))/12
                        mat[i][j]+=vol*(avg*second+weighted)
            for i in range(c):
                for j in range(i):
                    mat[i][j]=mat[j][i]
            moments.append(mat)
        fixed=potential[target] if isinstance(target,int) else F(1)
        return self.factor*fixed*mixed(moments)


def compare_expanded(bundle,placements,load,probes=()):
    # Replay the full 16-edge determinant with distinct physical insertion sites.
    edges=[];times=[];labels=[];parent=[]
    src=list(bundle.source)
    for e,(left,right) in enumerate(bundle.edges):
        points=[]
        for k,(a,b,s,t) in enumerate(placements):
            for j,(path,z) in enumerate(((a,s),(b,t))):
                if path==e:
                    points.append((z,k,j))
        last=F(0);previous=left
        for z,k,j in sorted(points):
            assert last<z<1
            v=len(src)
            src.append(tuple(bundle.pairs[k][j]*int(h==k) for h in range(5)))
            labels.append((k,j,v))
            times.append(bundle.times[e]*(z-last));parent.append(e)
            edges.append((previous,v));last=z;previous=v
        edges.append((previous,right));times.append(bundle.times[e]*(1-last));parent.append(e)
    cycle=[[ch[i] for i in parent] for ch in zip(*bundle.cycle)]
    flow,phi=solve_graph(edges,src,times,0)
    exact=cofactor(flow,cycle,times)**2
    rows=[bundle.projected(k,*p) for k,p in enumerate(placements)]
    projected=bundle.factor*F(str(fmpq_mat([[str(x) for x in row] for row in rows]).det()))**2
    assert exact==projected
    # The resistance-line interpolation includes a Brownian bridge for two
    # distinct original wave positions on the same reduced edge.
    observations={}
    for k,j,v in labels:
        e=placements[k][j];z=placements[k][2+j]
        a,b=bundle.edges[e]
        val=F(0)
        for w,row in enumerate(bundle.source):
            val+=dot(row,load)*((1-z)*bundle.green[a][w]+z*bundle.green[b][w])
        for h,l,w in labels:
            f=placements[h][l];s=placements[h][2+l]
            c,d=bundle.edges[f]
            green=((1-z)*((1-s)*bundle.green[a][c]+s*bundle.green[a][d])
                   +z*((1-s)*bundle.green[b][c]+s*bundle.green[b][d]))
            if e==f:
                green+=bundle.times[e]*(min(z,s)-z*s)
            val+=load[h]*bundle.pairs[h][l]*green
        assert dot(phi[v],load)==val
        observations[k,j]=val
    if probes:
        return tuple(exact*(observations[i]-observations[j]) for i,j in probes)
    return exact


def check():
    from labeled_bulk import subdivide
    from random import Random
    rng=Random(10091)
    pairs=((1,2),(1,3),(1,-1))
    edges=subdivide((2,2,2,2,1,1))
    source=[(0,)*5]*4+[(0,0,0,1,0),(0,0,0,0,1),(0,-2,-1,-1,1),(-3,-2,1,0,-2)]
    bundle=PartialCotrees(edges,source,range(2,12),pairs)
    load=(0,0,0,0,1)
    nonzero=0
    for i in range(20):
        ps=[(*rng.sample(range(len(edges)),2),F(2*k+1,13),F(2*k+2,13))
            for k in range(len(pairs))]
        nonzero+=bool(compare_expanded(bundle,ps,load))
    assert nonzero
    value=bundle.mass(load,(0,0))-bundle.mass(load,4)
    assert value==-F(6604126224963280696502022854122480,22071023253620227236325964673)
    print('Original 16-edge checks:',nonzero,'nonzero of 20; reduced voltage',value,flush=True)
    pairs=((1,2),(1,3),(1,4),(1,5))
    edges=subdivide((2,2,1,1,1,1))
    source=[(0,)*5]*4+[(0,0,0,0,1),(-3,-4,-5,-6,-1)]
    bundle=PartialCotrees(edges,source,range(2,10),pairs)
    probes=[((1,0),(1,1)),((1,0),(2,0)),((1,0),(2,1)),((1,0),(3,0))]
    nonzero=0
    for i in range(30):
        ps=[(*rng.sample(range(len(edges)),2),F(2*k+1,13),F(2*k+2,13))
            for k in range(len(pairs))]
        vals=compare_expanded(bundle,ps,(1,0,0,0,0),probes)
        nonzero+=bool(any(vals))
    assert nonzero
    print('All four nonplanar four-pair probes equal their original 16-edge values on 30 placements,',
          nonzero,'nonzero.',flush=True)


if __name__=='__main__':
    check()
