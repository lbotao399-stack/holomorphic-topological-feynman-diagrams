"""Five remaining position integrals of a bare quadratic lift at fixed core resistance."""
from fractions import Fraction as F
from itertools import product,permutations
import sys
from math import prod

from flint import fmpz_mat,fmpz_mpoly_ctx
from quartic_seven_plane_lift import NORMAL,cubic_cores,density,fixture,G,solve_graph
from marked_theta import green_moments

GRID=tuple(product(range(4),repeat=5))
CTX=fmpz_mpoly_ctx.get(tuple('s'+str(i) for i in range(5)))
# Six times the coefficient table for interpolation on 0,1,2,3.
COEFFICIENTS=((6,0,0,0),(-11,18,-9,2),(6,-15,12,-3),(-1,3,-3,1))


def interpolate(values):
    a=values.copy()
    stride=1
    for _ in range(5):
        for b in range(0,len(a),4*stride):
            for j in range(stride):
                ids=[b+j+i*stride for i in range(4)]
                original=[a[i] for i in ids]
                for n,i in enumerate(ids):
                    a[i]=sum(x*y for x,y in zip(COEFFICIENTS[n],original))
        stride*=4
    return CTX.from_dict({g:v for g,v in zip(GRID,a) if v})


def data(core,t,paths):
    cols=list(zip(*G(6,core).C))
    m=fmpz_mat([[sum(x[i]*x[j]*z for x,z in zip(cols,t)) for i in range(4)] for j in range(4)])
    u=int(m.det())
    adj=m.inv()*u
    adj=[[int(a) for a in row] for row in adj.tolist()]
    el=[tuple(int(i==j)-int(j==0) for i in range(6)) for j in range(6)]
    flow,phi=solve_graph(core,el,t,0)
    # Unit source moving on one edge, with an implicit balancing charge at branch 0.
    old=[[int(v*u) for v in row] for row in flow]
    phi=[[F(v) for v in phi[i]] for i in range(6)]
    a0=[];a1=[]
    for e in range(9):
        a0.append([old[e][core[h][0]]-int(e==h)*u for h in paths])
        a1.append([old[e][core[h][1]]-old[e][core[h][0]]+int(e==h)*u for h in paths])
    return u,adj,cols,a0,a1,phi


def sampled(q,z,paths,cols,a0,a1):
    rows=list(zip(*[(*n,k) for n,k in zip(q,z)]))
    jump=[[cols[h][r]*q[j][k] for j,h in enumerate(paths)] for k in range(3) for r in range(4)]
    last=[[cols[h][r]*z[j] for j,h in enumerate(paths)] for r in range(4)]
    # The five jump columns carry unscaled momenta; only eight initial currents remain after the projective row.
    out=[[] for _ in range(4)]
    for ss in GRID:
        weights=[[a+s*b for a,s,b in zip(x,ss,y)] for x,y in zip(a0,a1)]
        yy=[[sum(v*c for v,c in zip(qrow,w)) for w in weights] for qrow in rows]
        base=[[1]*9+[0]*5]+[[col[r]*v for col,v in zip(cols,yy[i])]+jump[4*i+r]
                              for i in range(3) for r in range(4)]
        for r in range(4):
            out[r].append(int(fmpz_mat([*base,[col[r]*v for col,v in zip(cols,yy[3])]+last[r]]).det()))
    return [interpolate(row) for row in out]


def position_value(core,t,paths,normal=NORMAL,slopes=(1,-1,0,0,0),order=(0,1,2,3,4),check=False):
    u,adj,cols,a0,a1,phi=data(core,t,paths)
    q=[normal[i] for i in order]
    z=[slopes[i] for i in order]
    b=sampled(q,z,paths,cols,a0,a1)
    for p in b:
        assert all(sum(k)<=8 and all(sum(k[j] for j,h in enumerate(paths) if h==l)<=3 for l in set(paths))
                    for k,_ in p.terms())
    sq=sum((adj[i][j]*(2 if i<j else 1)*b[i]*b[j] for i in range(4) for j in range(i,4)),CTX.constant(0))
    s=CTX.gens()
    # Endpoint part of Green is bilinear in fractions. Common-edge minima are integrated separately.
    target=((order.index(0),1),(order.index(1),-1))
    v=CTX.constant(0)
    same=[]
    for i,sign in target:
        a,b=core[paths[i]]
        for j,src in enumerate(q):
            amt=sign*src[2]
            if not amt: continue
            c,d=core[paths[j]]
            v+=amt*(int(phi[a][c]*u)*(1-s[i])*(1-s[j])+int(phi[a][d]*u)*(1-s[i])*s[j]
                   +int(phi[b][c]*u)*s[i]*(1-s[j])+int(phi[b][d]*u)*s[i]*s[j])
            if paths[i]==paths[j]:
                same.append((i,j,amt*t[paths[i]]*u))
    val=sum((F(int(c),prod(1+int(x) for x in k)) for k,c in (sq*v).terms()), F())
    if same:
        for k,c in sq.terms():
            powers=tuple(map(int,k))
            val+=int(c)*sum(charge*green_moments(powers,i,j)[1] for i,j,charge in same)
    denominator=6**10*u**17*sum(t)**2
    if check:
        point=tuple(F(j+1,7) for j in range(5))
        def at(poly):
            return sum(int(c)*prod(s**int(n) for s,n in zip(point,k)) for k,c in poly.terms())
        voltage=at(v)+sum(charge*(min(point[i],point[j])-point[i]*point[j]) for i,j,charge in same)
        lengths=[]
        for h in range(9):
            ss=[F(),*[point[i] for i,k in enumerate(paths) if k==h],F(1)]
            lengths.extend(t[h]*(y-x) for x,y in zip(ss,ss[1:]))
        actual=density(fixture(core,paths),[int(7*x) for x in lengths],normal,slopes,[order])[1]/49
        assert at(sq)*voltage/denominator==actual and actual
    return val/denominator


if __name__=='__main__':
    core=cubic_cores()[-2]
    t=(2,3,5,7,11,13,17,19,23)
    paths=(0,2,3,6,8)
    result=position_value(core,t,paths,check=True)
    assert result==F(252808507370760972247212584450918497521418781760035836043,
                     127344496989907310920692545473802432914160662657158316032000)
    shared=position_value(core,t,(0,0,2,5,8),order=(0,2,1,3,4),check=True)
    assert shared==F(343210157206078005358930457457215074932624962330617421,
                     69612516576115512529533097015562554435583452618708992000)
    print('single placement after all five position integrals',result,flush=True)
    if '--orbit' in sys.argv:
        orbit=sum(position_value(core,t,paths,order=o) for o in permutations(range(5)))
        assert orbit==F(-948679641062255329222385597741966332326829569020906127211,
                        541430684480898430785257421232153201165649075923292160000)
        print('full normal permutation sum after the five position integrals',orbit)
