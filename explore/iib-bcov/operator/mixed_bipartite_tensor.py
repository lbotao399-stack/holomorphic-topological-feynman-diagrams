"""Exact tensor bubbles in the one-star K3,3 Ward chamber."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json

from mixed_bipartite_one_star import HERE, TABLE
from quartic_prism_physical_ibp import compile_program, zeta_interval
from quartic_prism_completion import completed, prism_charts
from quartic_seven_plane_single import checked_cut, orbit, tetrahedron_bank


def pairs(n):
    return tuple(combinations(range(n),2))


def unpack(n,weights):
    matrix = [[0]*n for _ in range(n)]
    for (i,j),w in zip(pairs(n),weights,strict=True):
        matrix[i][j]=matrix[j][i]=w
    return matrix


@cache
def centered_moment(q,exponents):
    """Successive (Laplacian/4)^j on the product of squared distances."""
    rows = {(tuple(exponents),0):1}
    out = []
    for j in range(sum(exponents)+1):
        out.append(rows)
        following = defaultdict(int)
        for (ex,between),v in rows.items():
            for i,m in enumerate(ex):
                if m:
                    nxt = list(ex)
                    nxt[i] -= 1
                    following[tuple(nxt),between] += m*(q+sum(ex)-1)*v
            if len(ex)==2 and all(ex):
                following[tuple(m-1 for m in ex),between+1] -= prod(ex)*v
        rows = following
    return out


@cache
def translations(m):
    # |t A+(1-t)B-Y|^2=t|A-Y|^2+(1-t)|B-Y|^2-t(1-t)|A-B|^2.
    return [(u,v,m-u-v,(-1)**(m-u-v)*comb(m,u)*comb(m-u,v))
            for u in range(m+1) for v in range(m-u+1)]


def reduced_charts(matrix,q,vertex,infinity):
    n = len(matrix)
    A,B = [i for i in range(n) if i!=infinity and matrix[vertex][i]>0]
    others = [i for i in range(n) if i not in (vertex,infinity,A,B)]
    pts = [A,B,*others]
    finite = pairs(len(pts))
    slots = {p:i for i,p in enumerate(finite)}
    old = [matrix[pts[i]][pts[j]] for i,j in finite]
    b,c = matrix[vertex][A],matrix[vertex][B]
    powers = tuple(-matrix[vertex][h] for h in others)
    out = defaultdict(F)
    for j,rows in enumerate(centered_moment(q,powers)):
        radial = F(factorial(b+c-q-j-1),factorial(b-1)*factorial(c-1)*factorial(j))
        for (ex,between),coeff in rows.items():
            for term in product(*(translations(m) for m in ex)):
                edges = old.copy()
                edges[0] += b+c-q-j
                alpha,beta = q+j-c,q+j-b
                scalar = coeff
                if between:
                    edges[slots[2,3]] -= between
                for y,(u,v,w,k) in enumerate(term,2):
                    edges[slots[0,y]] -= u
                    edges[slots[1,y]] -= v
                    edges[0] -= w
                    alpha += u+w
                    beta += v+w
                    scalar *= k
                out[tuple(edges)] += radial*scalar*F(factorial(alpha-1)*factorial(beta-1),
                                                      factorial(alpha+beta-1))
    # Complete each chart at the same infinity. Combine coefficients before
    # testing the collision inequalities of the reduced period.
    result = {}
    for edges,weight in out.items():
        if not weight:
            continue
        mat = unpack(len(pts),edges)
        deficit = [2*q-sum(row) for row in mat]
        assert sum(deficit)==2*q
        matrix = [row+[d] for row,d in zip(mat,deficit)]+[deficit+[0]]
        cut = tuple(matrix[i][j] for i,j in pairs(len(pts)+1))
        if max(cut)>=q:
            return None
        result[cut] = weight
    return result


@cache
def period(q,n,weights):
    mat = unpack(n,weights)
    assert all(sum(row)==2*q for row in mat) and max(weights)<q
    if n==4:
        # Fix vertices 0,1 and infinity=3. Vertex 2 has two denominators.
        a,b = mat[0][2],mat[1][2]
        assert 0<a<q and 0<b<q and a+b>q
        value = F(factorial(a+b-q-1)*factorial(q-a-1)*factorial(q-b-1),
                  factorial(a-1)*factorial(b-1)*factorial(2*q-a-b-1))
        return value,F()
    if n==5:
        for inf in range(5):
            finite = list(combinations((i for i in range(5) if i!=inf),2))
            powers = tuple(mat[i][j] for i,j in reversed(finite))
            if min(powers)>0:
                key = orbit(powers)
                return tetrahedron_bank()[key] if q==13 else checked_cut(key)
    vertices = [(sum(max(-w,0) for w in row),v) for v,row in enumerate(mat)
                if sum(w>0 for w in row)==3]
    for _,v in sorted(vertices):
        for inf,power in enumerate(mat[v]):
            if power<=0:
                continue
            child = reduced_charts(mat,q,v,inf)
            if child is None:
                continue
            values = [F(),F()]
            for cut,w in child.items():
                col = period(q,n-1,cut)
                for j,val in enumerate(col):
                    values[j] += w*val
            return tuple(values)
    raise ValueError(f'No convergent contraction for {weights}')


def read_charts(rows,q,check=False):
    key_pairs = [(i,j) for j in range(2,6) for i in range(2)]+list(combinations(range(2,6),2))
    subtotal = [F(),F()]
    evaluations = {}
    checks = 0
    for line in rows.splitlines():
        row=line.split()
        assert len(row)==26
        graph=[[0]*6 for _ in range(6)]
        for (i,j),w in zip(key_pairs,map(int,row[11:25])):
            graph[i][j]=graph[j][i]=w
        chart=tuple(graph[i][j] for i,j in pairs(6))
        if not all(sum(graph[i][j] for i,j in combinations(t,2))<2*q for t in combinations(range(6),3)):
            raise ValueError('Divergent seed')
        val=period(q,6,chart)
        if check:
            for v,edges in enumerate(graph):
                if sum(w>0 for w in edges)!=3:
                    continue
                for inf,w in enumerate(edges):
                    if w<=0:
                        continue
                    child = reduced_charts(graph,q,v,inf)
                    if child is not None:
                        control = [F(),F()]
                        for cut,mult in child.items():
                            for j,x in enumerate(period(q,5,cut)):
                                control[j] += mult*x
                        assert tuple(control)==val
                        checks += 1
        coeff=F(row[25])
        result=tuple(coeff*x for x in val)
        evaluations[tuple(map(int,row[:9]))]=result
        for j,x in enumerate(result):
            subtotal[j]+=F(row[9])*x
    return evaluations,tuple(subtotal),checks


def bubble_charts(exe,q,records):
    data=''.join(' '.join(map(str,(*p,w,0)))+'\n' for p,w in records)
    return run([exe,str(q),'--bubble-charts'],input=data,text=True,capture_output=True,check=True).stdout


def lower_dimension_control(exe):
    reference=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    records=[]
    known={}
    for ex,row in zip(reference['exponents'],reference['coefficients']):
        powers=tuple(a+1 for a in ex)
        key=min(prism_charts(completed(powers,6),6))
        scale=F(factorial(5),prod(factorial(a) for a in ex))
        val=tuple(F(str(x))*scale for x in row)
        assert known.setdefault(key,val)==val
        records.append((*powers,1,0))
    data=''.join(' '.join(map(str,row))+'\n' for row in records)
    full=run([exe,'6','--full-components'],input=data,text=True,capture_output=True,check=True).stdout
    bubbles=[]
    for line in full.splitlines():
        row=line.split()
        if row[11]=='1':
            bubbles.append((tuple(map(int,row[:9])),1))
    measured,_,checks=read_charts(bubble_charts(exe,6,bubbles),6,check=True)
    for key,val in measured.items():
        assert (*val,0)==known[key]
    print('q=6 independently known bubble cuts and convergent contractions:',len(measured),checks)


def physical(exe,check=False):
    records=[]
    scalar=[F(),F()]
    with gzip.open(TABLE,'rt') as f:
        for line in f:
            row=line.split()
            if row[11]!='1':
                continue
            if len(row)==19:
                chart=orbit(tuple(map(int,row[12:18])))
                for i,v in enumerate(tetrahedron_bank()[chart]):
                    scalar[i]+=F(row[9])*F(row[18])*v
            else:
                records.append((tuple(map(int,row[:9])),F(row[9])))
    measured,tensor,checks=read_charts(bubble_charts(exe,13,records),13,check=check)
    assert len(measured)==1242 and tensor==(
        -F(4197421754214856126954411,58070236078413279589220352000000000),
        F(3927607033203397,65316736041839073792000000))
    print('1242 tensor components, rational and zeta(3) coefficients:',*tensor)
    a,b = (x+y for x,y in zip(scalar,tensor))
    lo,hi = zeta_interval(3)
    assert -17682 < 10**20*(a+b*hi) < 10**20*(a+b*lo) < -17681
    print('All 1828 bubble components:',a,b)
    if check:
        print('q=13 convergent contractions:',checks)
    return measured


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Check every convergent placement of infinity and the q=6 graph bank')
    args=parser.parse_args()
    with TemporaryDirectory() as scratch:
        exe=Path(scratch)/'complete'
        compile_program('quartic_prism_completion.cpp',exe)
        if args.check:
            lower_dimension_control(exe)
        physical(exe,check=args.check)
