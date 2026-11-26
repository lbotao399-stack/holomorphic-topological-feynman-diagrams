"""Real Gaussian marked-edge form, compared with the S2 cotree-current chart."""
from pathlib import Path
from random import Random
import sys
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from q3lib import ALL


def parity(items):
    return (-1)**sum(a>b for j,a in enumerate(items) for b in items[j+1:])


def marked(edges,root,t,lam):
    edges_count=len(edges)
    b=sp.Matrix([[int(i==v)-int(j==v) for v in range(root)] for i,j in edges])
    inv_t=sp.diag(*[1/sp.Integer(x) for x in t])
    laplace=b.T*inv_t*b
    inverse=laplace.inv()
    currents=inv_t*b*inverse*sp.Matrix(lam)
    u=laplace.det()*sp.prod(t)
    # Fundamental cotree chosen directly from the first independent incidence rows.
    tree=list(b.T.rref()[1]);cotree=[e for e in range(edges_count) if e not in tree]
    sigma=parity(tree+cotree)*b[tree,:].det()
    # Differentiate t*y=B*phi with lambda held fixed.
    response=inv_t*b*inverse*b.T*inv_t-inv_t
    jac=sp.Matrix([[response[c,e]*currents[e,k] for e in range(edges_count)] for c in cotree for k in range(3)])
    raw=[];chart=[]
    for j in range(edges_count):
        rows=[];factor=sp.Integer(1)
        for e in range(edges_count):
            diff=[[b[e,v]*int(s==k) for v in range(root) for k in range(3)] for s in range(3)]
            if e==j:
                rows.extend(diff)
            else:
                a=next((i for i in range(3) if currents[e,i]),None)
                if a is None:
                    factor=0;break
                factor*=(-1)**a/currents[e,a]
                rows.extend([[currents[e,a]*x-currents[e,s]*y for x,y in zip(diff[s],diff[a])]
                             for s in range(3) if s!=a])
        exterior=factor*sp.Matrix(rows).det(method='domain-ge') if factor else 0
        raw.append(exterior)
        facet=sigma*(-1)**j*jac[:,[e for e in range(edges_count) if e!=j]].det(method='domain-ge')
        chart.append(facet)
        assert (-1)**(root*(root-1)//2)*exterior/u**3==facet
    return raw,chart


def verify():
    box=ALL['box']['edges']
    raw,facets=marked(box,3,[1]*4,sp.eye(3).tolist())
    assert raw==[-sp.Rational(1,4)]*4 and facets==[sp.Rational(1,256)]*4
    rng=Random(150301)
    count=4
    for name in ('A','B','G1','G2','G3','G4','G5'):
        edges=ALL[name]['edges'];root=max(max(e) for e in edges)
        valence=[sum(v in e for e in edges) for v in range(root)]
        times=[rng.randrange(1,9) for _ in edges]
        lam=[[rng.randrange(-4,5) for _ in range(3)] if d==2 else [0]*3 for d in valence]
        _,values=marked(edges,root,times,lam)
        assert any(values)
        count+=len(edges)
    print('all positive marked densities equal the S2 chart;',count,'exact facet comparisons')


if __name__=='__main__':
    verify()
