"""Degree-nine invariant bracket boundaries on three F4 adjoint colors."""
from itertools import combinations, product
from functools import cache
from collections import Counter
from f4 import triple_invariants
import numpy as np
PRIME = 1009

def oct_table():
    tab={(0,i):(i,1) for i in range(8)}
    tab.update({(i,0):(i,1) for i in range(8)})
    for i in range(1,8):tab[i,i]=(0,-1)
    for a,b,c in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
        for i,j,k in [(a,b,c),(b,c,a),(c,a,b)]:
            tab[i,j]=(k,1)
            tab[j,i]=(k,-1)
    return tab

def jordan():
    tab=oct_table()
    basis=[{(i,i,0):1} for i in range(3)]
    for i,j in combinations(range(3),2):
        for k in range(8):
            basis.append({(i,j,k):1,(j,i,k):1 if not k else -1})
    inds=[next(iter(x)) for x in basis]
    left=np.zeros((27,27,27),dtype=np.int32)
    for m,n in product(range(27),repeat=2):
        vals=Counter()
        for a,b in [(m,n),(n,m)]:
            for (i,j,k),x in basis[a].items():
                for (h,l,t),y in basis[b].items():
                    if j==h:
                        u,c=tab[k,t]
                        vals[i,l,u]+=c*x*y
        left[m,:,n]=[vals[key] for key in inds]
    return left

def echelon(v,piv):
    v=v.copy()%PRIME
    for i,w in piv:
        if v[i]: v=(v-v[i]*w)%PRIME
    if np.any(v):
        i=int(np.flatnonzero(v)[0])
        w=v*pow(int(v[i]),-1,PRIME)%PRIME
        piv.append((i,w))
        return True
    return False

def sample_basis(left):
    piv=[]
    basis=[]
    for a,b in combinations(left,2):
        v=a@b-b@a
        if echelon(v.flatten(),piv): basis.append(v)
    assert len(basis)==52
    # Jordan identity element is fixed.
    identity=np.zeros(27,dtype=np.int32)
    identity[:3]=1
    assert not np.any(np.array(basis)@identity)
    return np.array(basis)

def necklaces(w):
    return [w[j:]+w[:j] for j in range(len(w))]

def words(target=(3,3,3)):
    by_weight={}
    for n in range(2,10):
        for w in product(range(3), repeat=n):
            cnt=tuple(w.count(j) for j in range(3))
            if any(a>b for a,b in zip(cnt,target)):continue
            shifts=necklaces(w)
            if w!=min(shifts):continue
            rev=min(necklaces(w[::-1]))
            if w>rev or (w==rev and n%2):continue
            by_weight.setdefault(cnt,[]).append(w)
    flat=[w for ws in by_weight.values() for w in ws]
    mono={(0,0,0):[()]}
    for i,w in enumerate(flat):
        delta=tuple(w.count(j) for j in range(3))
        for p in product(*(range(b+1-x) for x,b in zip(delta,target))):
            dest=tuple(x+y for x,y in zip(p,delta))
            for prev in mono.get(p,()):
                mono.setdefault(dest,[]).append((*prev,i))
    return flat, mono[target]


def traces(M,flat):
    @cache
    def matrix(w):
        if len(w)==1: return M[w[0]]
        return (matrix(w[:-1])@M[w[-1]])%PRIME
    out=[]
    for w in flat:
        k=len(w)//2
        out.append(int(np.sum(matrix(w[:k])*matrix(w[k:]).T))%PRIME)
    matrix.cache_clear()
    return np.array(out,dtype=np.int32)



def certificate():
    left = jordan()
    basis = sample_basis(left)
    metric = np.array([1]*3+[2]*24,dtype=np.int32)
    cubic = left * metric.reshape(1,27,1)
    assert np.array_equal(cubic,cubic.transpose(1,0,2))
    assert np.array_equal(cubic,cubic.transpose(0,2,1))
    for B in basis:
        assert not np.any(B.T * metric.reshape(1,27) + metric.reshape(27,1) * B)
        # D(u*v) = D(u)*v + u*D(v), with left multiplication basis L_u.
        rhs = np.einsum("ai,ajt->ijt", B, left, optimize=True)
        assert np.array_equal(B@left-left@B, rhs)
    indices = np.nonzero(cubic)
    i,j,k = indices
    down = cubic[indices]
    denominators = metric[i]*metric[j]*metric[k]
    up = down * np.array([pow(int(a),-1,PRIME) for a in denominators],dtype=np.int32)%PRIME
    weights = (up[:,None]*down[None,:])%PRIME
    wires = (
        ((0,1),(1,0,0),(2,1,2,2)),
        ((1,),(2,0,2),(1,2,0,1,0)),
        ((2,),(0,2,1,0),(1,0,2,1)),
        ((0,0),(1,2),(1,2,0,2,1)),
    )
    flat, monomials = words()
    rng = np.random.default_rng(191)
    pivots, trace_pivots = [], []
    for probe in range(115):
        colors = rng.integers(0,PRIME,size=(3,52),dtype=np.int32)
        matrices = (colors@basis.reshape(52,-1)).reshape(3,27,27)%PRIME
        values = traces(matrices,flat)
        row = [np.prod(values[list(m)],dtype=np.int64)%PRIME for m in monomials]
        echelon(np.array(row,dtype=np.int32),trace_pivots)
        @cache
        def matrix(letters):
            if len(letters)==1:
                return matrices[letters[0]]
            return matrix(letters[:-1])@matrices[letters[-1]]%PRIME
        for first,second,third in wires:
            term = weights*matrix(first)[i[:,None],i[None,:]]%PRIME
            term = term*matrix(second)[j[:,None],j[None,:]]%PRIME
            term = term*matrix(third)[k[:,None],k[None,:]]%PRIME
            row.append(np.sum(term,dtype=np.int64)%PRIME)
        matrix.cache_clear()
        echelon(np.array(row,dtype=np.int32),pivots)
        if len(pivots)==107:
            break
    assert len(trace_pivots)==103 and len(pivots)==triple_invariants()
    return len(trace_pivots),len(pivots),probe+1


if __name__=="__main__":
    print("F4 commuting-triple boundaries: trace rank, combined rank, evaluation points",certificate())
