"""Integrated three-loop BCOV potential on a labeled cubic graph at arbitrary momenta."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import networkx as nx

from bulk_three import forest_data, integrate_core
from square_three import EDGES, G, skeleton_chains, tree_flows
from current_minor import numerator
from path_moments import PathMoments, add_scaled


def subdivide(sizes):
    edges = []
    start = 4
    for (u,v),m in zip(EDGES,sizes):
        row = [u,*range(start,start+m-1),v]
        edges.extend(zip(row,row[1:]))
        start += m-1
    return edges


class LabeledBulk:
    def __init__(self,edges,source,all_orders=False):
        n = len(edges)-len(source)+1
        h = len(source[0]) if source else 0
        if (not h or n not in (2,3) or len(edges)!=n*h+1 or any(len(q)!=h for q in source)
                or any(map(sum,zip(*source)))):
            raise ValueError('Expected a balanced two- or three-loop tight graph.')
        if any(not 0<=v<len(source) for e in edges for v in e):
            raise ValueError('A graph endpoint is outside the source array.')
        branch,paths = skeleton_chains(edges,len(source))
        if (len(branch)!=2*n-2 or any(any(source[v]) for v in branch)
                or sum(len(es) for _,es in paths)!=len(edges)):
            raise ValueError('Every branch must be trivalent and source-free.')
        ordered = []
        for vertices,indices in paths:
            if vertices[0]>vertices[-1]:
                vertices,indices = vertices[::-1],indices[::-1]
            pair = branch.index(vertices[0]),branch.index(vertices[-1])
            ordered.append((pair,vertices,indices))
        ordered.sort()
        self.edges = [e for e,_,_ in ordered]
        self.sites = {v:(i,j) for i,(_,verts,_) in enumerate(ordered)
                      for j,v in enumerate(verts[1:-1],1)}
        if len({v for vertices,indices in paths for v in vertices})!=len(source):
            raise ValueError('Every vertex must belong to a cubic core path.')
        g = G(len(source),edges)
        denominator = lcm(*(F(x).denominator for q in source for x in q))
        integral = [tuple(int(F(x)*denominator) for x in q) for q in source]
        # One particular Kirchhoff flow is enough. No tree sum enters the minor.
        tree = nx.Graph()
        tree.add_edges_from((a,b,{'index':i}) for i,(a,b) in enumerate(edges))
        flow = [[0]*h for _ in edges]
        for v,path in nx.single_source_shortest_path(tree,0).items():
            for a,b in zip(path,path[1:]):
                i = tree[a][b]['index']
                sign = 1 if edges[i]==(b,a) else -1
                for k,q in enumerate(integral[v]):
                    flow[i][k] += sign*q
        poly = numerator(g,integral,flow)
        self.moments = PathMoments(poly,[es for _,_,es in ordered],n-1,all_orders)
        if denominator!=1:
            scale = denominator**(n*h)
            self.moments.base = {key:value/scale for key,value in self.moments.base.items()}
        self.ctx,self.u,self.forests = forest_data(self.edges,len(branch))
        # forest_data uses x names; substitute the identical core variable order.
        new = self.moments.ctx
        self.u = new.from_dict(self.u.to_dict())
        self.forests = [[new.from_dict(entry.to_dict()) for entry in row] for row in self.forests]
        self.ctx = new
        self.n,self.h,self.source = n,h,source

    def voltage(self,left,right,vector):
        """Potential(left)-potential(right), paired with the given vector."""
        polynomial = self.voltage_polynomial(left,right,vector)
        return integrate_core(self.edges,polynomial,self.n,2*self.n+self.h+1)

    def voltage_polynomial(self,left,right,vector):
        """Numerator after the path-position integral, before the core period."""
        space = self.moments
        if left not in self.sites or right not in self.sites or len(vector)!=len(self.source[0]):
            raise ValueError('Expected two bivalent sites and a vector dual to a source.')
        nb = len(self.forests)
        core_target = [{} for _ in range(nb)]
        core_load = [{} for _ in range(nb)]
        loads = {v:sum(a*b for a,b in zip(self.source[v],vector)) for v in self.sites}
        starts = [[0]*nb for _ in range(2)]
        for v,(i,position) in self.sites.items():
            a,b = self.edges[i]
            s = space.position(i,position)
            weights = int(v==left)-int(v==right),loads[v]
            for start,fields,weight in zip(starts,(core_target,core_load),weights):
                start[a] += weight
                add_scaled(fields[a],s,-weight)
                add_scaled(fields[b],s,weight)
        for start,fields in zip(starts,(core_target,core_load)):
            for a,weight in enumerate(start):
                add_scaled(fields[a],space.base,weight)
        result = self.ctx.constant(0)
        for a,b in combinations(range(nb),2):
            result -= self.forests[a][b]*(space.inner(core_target[a],core_load[b])
                                           +space.inner(core_target[b],core_load[a]))/2
        for v,sign in ((left,1),(right,-1)):
            i,at = self.sites[v]
            earlier,linear = {},{}
            for w,(j,pos) in self.sites.items():
                if j==i and loads[w]:
                    add_scaled(earlier,space.minimum(i,at,pos),loads[w])
                    add_scaled(linear,space.position(i,pos),loads[w])
            bridge = space.inner(space.base,earlier)-space.inner(space.position(i,at),linear)
            result += sign*self.ctx.gen(i)*self.u*bridge
        return space.volume*result

    def resistance_pairing(self,weights,cotree=False):
        """Integrated sum of pair resistances; weights have zero sum at every vertex."""
        if not self.moments.base:
            return F(0),F(0)
        vertices = set(self.sites) | ({-1} if cotree else set())
        assert all(i in vertices and j in vertices and i!=j for i,j in weights)
        assert all(sum(w for pair,w in weights.items() if v in pair)==0 for v in vertices)
        space = self.moments
        base = space.inner(space.base,space.base)
        pos = {v:space.position(*site) for v,site in self.sites.items()}
        first = {v:space.inner(space.base,l) for v,l in pos.items()}
        f = self.forests
        result = self.ctx.constant(0)
        for (i,j),weight in weights.items():
            if cotree and i==-1:
                continue
            p,s = self.sites[i]
            q,t = self.sites[j]
            a,b = self.edges[p]
            c,d = self.edges[q]
            mixed = space.inner(pos[i],pos[j])
            pair = (f[a][c]*base+(f[b][c]-f[a][c])*first[i]+(f[a][d]-f[a][c])*first[j]
                    +(f[b][d]-f[a][d]-f[b][c]+f[a][c])*mixed)
            if p==q:
                minimum = space.inner(space.base,space.minimum(p,s,t))
                pair -= 2*self.ctx.gen(p)*self.u*(minimum-mixed)
            result += pair*weight.numerator/weight.denominator
        if cotree:
            # A uniformly chosen member of a Kirchhoff-weighted cotree is
            # uniform along path p with mass T_p (dU/dT_p)/(n U).
            result *= self.u
            for (i,j),weight in weights.items():
                if i!=-1:
                    continue
                q,_ = self.sites[j]
                c,d = self.edges[q]
                for p,(a,b) in enumerate(self.edges):
                    left = f[a][c]+f[b][c]
                    right = f[a][d]+f[b][d]
                    pair = (left*base+(right-left)*first[j])/2
                    if p==q:
                        pair -= self.u*self.ctx.gen(p)*(first[j]-space.inner(pos[j],pos[j]))
                    mass = self.ctx.gen(p)*self.u.derivative(p)
                    result += mass*pair*weight.numerator/(weight.denominator*self.n)
        return integrate_core(self.edges,space.volume*result,self.n,2*self.n+self.h+1+int(cotree))


def check_ordered():
    from square_three import integral
    from kernel import numerator as interpolation
    edges = subdivide((2,1,2,1,2,2))
    source = [(0,0,0)]*4+[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)]
    graph = G(8,edges)
    ctx,u,f = forest_data(edges,8)
    raw = ctx.from_dict(interpolation(graph,source,tree_flows(graph,source)))
    direct = -raw**2*sum(q[0]*(f[4][v]-f[5][v]) for v,q in enumerate(source))
    exact = tuple(v/2 for v in integral(graph,direct))
    assert any(exact)
    half = [tuple(F(x,2) for x in q) for q in source]
    assert LabeledBulk(edges,half).voltage(4,5,(1,0,0))==tuple(x/2**19 for x in exact)
    assert LabeledBulk(edges,source).voltage(4,5,(1,0,0))==exact
    print('Independent ordered threefold period:',*exact,flush=True)
    edges = subdivide((3,1,2,2,1,1))
    changed = source[:]
    changed[4],changed[5] = changed[5],changed[4]
    rows = (LabeledBulk(edges,source).voltage(4,6,(1,0,0)),
            LabeledBulk(edges,changed).voltage(5,6,(1,0,0)))
    cube = LabeledBulk(edges,source,all_orders=True).voltage(4,6,(1,0,0))
    assert cube==tuple(map(sum,zip(*rows)))==(F(7870438021,4413400992000000),
                                               -F(10097,6810804000))
    print('Independent sum of path orders:',*cube,flush=True)


if __name__=='__main__':
    check_ordered()
    edges = subdivide((3,3,3,3,2,2))
    p = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    p += [(1,1,0,1,0),(0,1,-1,1,1),(1,0,1,2,-1),(1,-1,0,0,1)]
    p.append(tuple(-sum(col) for col in zip(*p)))
    full = LabeledBulk(edges,[(0,)*5]*4+p)
    value = full.voltage(4,5,(0,0,1,0,0))
    assert value==(F(1284691025791082490901,67085096068410900480000000),
                   -F(84144127971551,5281958307225600000))
    print('Ten-wave ordered fivefold period:',*value,flush=True)
