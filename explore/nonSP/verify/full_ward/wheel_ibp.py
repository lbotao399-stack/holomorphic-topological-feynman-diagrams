"""Certify the eight wheel integration-by-parts rows used by the Ward bootstrap."""
from collections import Counter
from functools import cache
from itertools import product
from math import factorial, lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import subprocess

import numpy as np
from flint import fmpq_mat, fmpz_mat, nmod_mat
from scipy.sparse import csr_matrix

from contractions import HERE, edge_order, types
from cubature import grid_input
from homogeneous import reconstruct
from joint_wheels import CORE, DUAL, ORBIT, presentation
from lower_faces import G, tree_flows
from polynomials import compile_grid


class WheelMoments:
    def __init__(self):
        self.cotrees=tuple(tuple(int(i in miss) for i in range(8))
                           for miss in tree_flows(G(5,CORE),np.eye(5,dtype=int)))
        self.bounds=[]
        for mask in range(1,255):
            chosen=tuple(i for i in range(8) if mask>>i&1)
            order=min(sum(m[i] for i in chosen) for m in self.cotrees)
            if order:
                self.bounds.append((chosen,order))
        self.perms=tuple(tuple(map(int,p)) for p in ORBIT)
        self.dual=tuple(tuple(DUAL[i] for i in p) for p in self.perms)
        self.index={}
        for e in product(range(5),repeat=8):
            if sum(e)==16 and self.converges(e,6):
                self.index.setdefault(self.orbit(e),len(self.index))
        assert len(self.index)==1292

    @cache
    def orbit(self,e):
        return min([tuple(e[i] for i in p) for p in self.perms]+
                   [tuple(4-e[i] for i in p) for p in self.dual])

    def converges(self,e,q):
        return all(sum(e[i]+1 for i in chosen)>q*order for chosen,order in self.bounds)

    def relations(self):
        rows=[]
        for e in product(range(4),repeat=8):
            if sum(e)!=12 or min(tuple(e[i] for i in p) for p in self.perms)!=e:
                continue
            if not self.converges(e,5):
                continue
            columns=[]
            for m in self.cotrees:
                alpha=tuple(a+b for a,b in zip(e,m))
                assert self.converges(alpha,6)
                columns.append(self.index[self.orbit(alpha)])
            for j in range(8):
                row=Counter()
                for k,m in zip(columns,self.cotrees):
                    row[k]+=e[j]+1-5*m[j]
                rows.append({k:v for k,v in row.items() if v})
        table=np.zeros((len(rows),len(self.index)),dtype=np.int16)
        for j,row in enumerate(rows):
            for k,v in row.items():
                table[j,k]=v
        return table


def kernel(table):
    sparse=csr_matrix(table,dtype=np.int64)
    norm=int(abs(sparse).sum(axis=1).max())
    modulus=1
    for prime in (65521,65519,65497,65479,65449,65447):
        reduced,rank=nmod_mat(table.tolist(),prime).rref()
        assert rank==1023
        pivots=[]
        for j in range(rank):
            start=pivots[-1]+1 if pivots else 0
            pivots.append(next(i for i in range(start,reduced.ncols()) if reduced[j,i]))
        free=sorted(set(range(table.shape[1]))-set(pivots))
        vectors=np.zeros((table.shape[1],len(free)),dtype=np.int64)
        vectors[free]=np.eye(len(free),dtype=np.int64)
        for j,i in enumerate(pivots):
            vectors[i]=[-int(reduced[j,k])%prime for k in free]
        if modulus==1:
            reference=free
            crt=np.zeros(vectors.shape,dtype=object)
        assert free==reference
        crt+=modulus*(((vectors-crt)*pow(modulus,-1,prime))%prime)
        modulus*=prime
        try:
            lookup={i:reconstruct(int(i),modulus) for i in set(crt.flat)}
        except ValueError:
            continue
        rational=np.vectorize(lookup.__getitem__,otypes=[object])(crt)
        den=[lcm(*(x.denominator for x in col)) for col in rational.T]
        lifted=rational*np.array(den,dtype=object)
        if max(abs(x) for x in lifted.flat)*norm>=2**63:
            continue
        integer=np.array(lifted,dtype=np.int64)
        assert np.array_equal(integer[free],np.diag(den))
        if not np.any(sparse@integer):
            return fmpz_mat(integer.tolist())
    raise ArithmeticError('Failed to lift the 269-dimensional integration-by-parts kernel.')


def wheel_coefficients(moments,graph,probes,folder,binary):
    labels,a,b,_,cotrees,flow=presentation(graph)
    _,records=grid_input([(0,0)]*6,edge_order(graph))
    minors=[r for r in records if len(r)==9]
    path=folder/'edges'
    with path.open('w') as f:
        print(len(cotrees),len(minors),file=f)
        for miss,tree in zip(cotrees,flow.reshape(-1,9,5)):
            print(*miss,*tree[1:].ravel(),file=f)
        for minor in minors:
            print(*minor,file=f)
    sources=folder/'probes'
    sources.write_text(str(len(probes))+'\n'+'\n'.join(' '.join(str(x) for pair in p[:5] for x in pair)
                                                        for p in probes)+'\n')
    output=folder/'grid'
    subprocess.run([binary,path,sources,output,'--grid'],check=True)
    array=np.fromfile(output,dtype=np.int64).reshape((4,)*8+(len(probes),))
    vandermonde=fmpq_mat([[j**i for i in range(4)] for j in range(1,5)])
    interp=np.array([[int(6*v) for v in row] for row in vandermonde.inv().tolist()],dtype=np.int64)
    norm=int(abs(interp).sum(axis=1).max())
    for axis in range(8):
        values=np.moveaxis(array,axis,0)
        assert int(np.max(np.abs(values)))*norm<2**63
        coefficients=interp@values.reshape(4,-1)
        assert not np.any(coefficients%6)
        array=np.moveaxis((coefficients//6).reshape(values.shape),0,axis)
    # A common factor 7! clears every beta moment in the subdivided edge.
    result=np.zeros((len(probes),len(moments.index)),dtype=object)
    for tail in product(range(4),repeat=8):
        values=array[tail]
        if not np.any(values):
            continue
        e=(15-sum(tail),)+tail
        assert 0<=e[0]<=3
        alpha=[0]*8
        for j,v in enumerate(e):
            alpha[labels[j]]+=v
        alpha[labels[a]]+=1
        assert moments.converges(alpha,6)
        column=moments.index[moments.orbit(tuple(alpha))]
        beta,rem=divmod(factorial(7)*factorial(e[a])*factorial(e[b]),factorial(e[a]+e[b]+1))
        assert rem==0
        result[:,column]+=beta*values.astype(object)
    return fmpz_mat(result.tolist())


def joint_test(allowed,moments,binary,folder,fixture):
    from homogeneous import exact_kernel
    from joint_soft import remaining
    from sampling import evaluation, sampling
    from contractions import substitution
    from itertools import combinations
    probes=sampling()[0][:88]
    graphs=types(6)
    invariant,_,_,lift=exact_kernel()
    cuts=np.cumsum([0]+[i.shape[1] for i in invariant])
    families=np.stack([g@lift[a:b] for g,a,b in zip(invariant,cuts,cuts[1:])])
    common=remaining(families)
    kernels=[fmpz_mat(k.tolist())*common for k in families]
    test=fmpz_mat(evaluation(probes).tolist())
    combined,ward,images=[],[],[]
    for i,upper in ((8,11),(10,12)):
        image=wheel_coefficients(moments,graphs[i],probes,folder,binary)*allowed
        homogeneous=test*kernels[i]
        both=[r+s for r,s in zip(image.tolist(),homogeneous.tolist())]
        assert (image.rank(),homogeneous.rank(),fmpz_mat(both).rank())==(4,10,upper)
        images.extend(image.tolist())
        ward.extend(homogeneous.tolist())
        combined.extend(both)
    assert (fmpz_mat(images).rank(),fmpz_mat(ward).rank(),fmpz_mat(combined).rank())==(4,10,12)
    rows=fmpz_mat(fixture['relations'])*fmpz_mat(evaluation(fixture['probes']).tolist())*kernels[10]
    null,dim=rows.nullspace()
    assert dim==2
    free=fmpz_mat([r[:dim] for r in null.tolist()])
    for i in (8,10,11,12):
        part=kernels[i]*free
        for size in (3,4,5):
            rows=[]
            for support in combinations(range(6),size):
                image=tuple(tuple(int(v==j)-int(v==support[-1]) for j in support[:-1])
                            for v in range(5))
                rows.extend((fmpz_mat(substitution(image).toarray().tolist())*part).tolist())
            assert fmpz_mat(rows).rank()==size-3
    print('Joint two-wheel response rank 4. Its intersection with 10 Ward directions has rank 2.')
    print('All four remaining graph restrictions have ranks 0, 1, 2 at three, four, five marks.')


def check(joint=False):
    certificate=json.loads((HERE/'wheel_ibp.json').read_text())
    graph=types(6)[certificate['catalogue_index']]
    assert certificate['edges']==list(map(list,edge_order(graph)))
    moments=WheelMoments()
    relations=moments.relations()
    print('projective moment orbits',len(moments.index),'IBP rows',len(relations),flush=True)
    allowed=kernel(relations)
    with TemporaryDirectory() as path:
        folder=Path(path)
        binary=compile_grid(folder)
        coeff=wheel_coefficients(moments,graph,certificate['probes'],folder,binary)
        if joint:
            joint_test(allowed,moments,binary,folder,certificate)
    response=coeff*allowed
    weights=fmpz_mat(certificate['relations'])
    assert response.rank()==4 and weights.rank()==8 and (weights*response).is_zero()
    print('Rational IBP rank 1023, nullity 269, wheel response rank 4, eight exact source relations.')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--joint',action='store_true',help='also verify both wheel images and the marked-source ranks')
    check(parser.parse_args().joint)
