"""Finite-field cubature for symmetric homogeneous three-loop core numerators."""
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from pathlib import Path
from random import Random
import numpy as np
from flint import nmod_mat

from fourfold_bulk import CORES,symmetries
from bulk_three import core_moment


@cache
def supports(core,q):
    # The all-edge equality is homogeneity. Each proper cyclic subset has a
    # strict inequality from convergence under its independent rescaling.
    restrictions = []
    for size in range(2,6):
        for subset in combinations(range(6),size):
            parts = [{i} for i in range(4)]
            for j in subset:
                a,b = core[j]
                left = next(s for s in parts if a in s)
                right = next(s for s in parts if b in s)
                if left!=right:
                    left.update(right)
                    parts.remove(right)
            loops = size-4+len(parts)
            if loops:
                restrictions.append((subset,q*loops))
    group = sorted(set(symmetries(core)))
    seen = set()
    orbits = []
    for k in product(range(1,q),repeat=5):
        last = 3*q-sum(k)
        if not 0<last<q:continue
        p = (*k,last)
        if p in seen or any(sum(p[j] for j in s)<=bound for s,bound in restrictions):
            continue
        orbit = tuple(sorted({tuple(p[j]-1 for j in perm) for perm in group}))
        seen.update(tuple(k+1 for k in e) for e in orbit)
        orbits.append(orbit)
    return orbits


def kirchhoff(c):
    cuts = []
    for keep in combinations(range(6),3):
        reached = {0}
        for _ in range(3):
            for j in keep:
                a,b = c[j]
                if a in reached or b in reached:reached.update((a,b))
        if len(reached)==4:cuts.append(tuple(i for i in range(6) if i not in keep))
    return cuts


def random_nodes(core,q,prime,seed):
    rng = Random(seed)
    cuts = kirchhoff(core)
    rows = []
    while len(rows)<len(supports(core,q)):
        t = [rng.randrange(1,prime) for i in range(5)]+[1]
        if sum(t[i]*t[j]*t[k] for i,j,k in cuts)%prime:
            rows.append(t)
    return rows


def evaluations(t,orbits,q,prime):
    if 24*prime**2 >= 2**63:
        raise ValueError('Prime exceeds the dot-product range.')
    x = np.array(t,dtype=np.int64).T
    limit = q-1
    powers = np.empty((6,limit,len(t)),dtype=np.int64)
    powers[:,0] = 1
    for k in range(1,limit):
        powers[:,k] = powers[:,k-1]*x%prime
    def triple(start):
        ab = powers[start,:,None,:]*powers[start+1,None,:,:]%prime
        return (ab[:,:,None,:]*powers[start+2,None,None,:,:]%prime).reshape(limit**3,-1)
    left = triple(0)
    # The last time is one; its power is fixed by total degree.
    right = (powers[3,:,None,:]*powers[4,None,:,:]%prime).reshape(limit**2,-1)
    ans = np.empty((len(orbits),len(t)),dtype=np.int64)
    for i,orbit in enumerate(orbits):
        a,b,c,d,e,_ = np.asarray(orbit,dtype=np.int64).T
        ans[i] = np.sum(left[(a*limit+b)*limit+c]*right[d*limit+e],axis=0)%prime
    return nmod_mat(ans.tolist(),prime)


def rational_table(core,q,path=None):
    table = {}
    if path and path.exists():
        for line in path.read_text().splitlines():
            fields = line.split()
            key = tuple(map(int,fields[:6]))
            if len(fields)!=8 or key in table:raise ValueError('Invalid period table.')
            table[key] = tuple(map(F,fields[6:]))
    orbits = supports(core,q)
    periods = []
    for orbit in orbits:
        key = orbit[0]
        if key in table:
            value = table.pop(key)
        else:
            value = core_moment(core,3,q,key)
            if core==CORES['tetrahedron']:
                from kite_native import coefficients
                if coefficients.cache_info().currsize>18000:
                    coefficients.cache_clear()
            if path:
                with path.open('a') as out:print(*key,*value,sep='\t',file=out)
        periods.append(tuple(len(orbit)*v for v in value))
    if table:raise ValueError('The cached monomial support differs from the chosen core.')
    return periods


def residue(v,prime):
    return v.numerator*pow(v.denominator,-1,prime)%prime


def rule(core,q,prime,periods,seed=310321):
    from flint import fmpz
    if not fmpz(prime).is_prime() or prime<=q:
        raise ValueError('Expected a prime larger than the denominator power.')
    orbits = supports(core,q)
    target = nmod_mat([[residue(x,prime) for x in row] for row in periods],prime)
    for j in range(8):
        nodes = random_nodes(core,q,prime,seed+j)
        matrix = evaluations(nodes,orbits,q,prime)
        try:
            weights = matrix.solve(target)
        except ZeroDivisionError:
            continue
        assert matrix*weights==target
        return [(*node,int(weights[i,0]),int(weights[i,1])) for i,node in enumerate(nodes)]
    raise ArithmeticError('No nonsingular cubature evaluation within eight seeds.')


def graph_check(core,prime,rows):
    from flint import nmod_mpoly_ctx
    from labeled_bulk import LabeledBulk
    from bulk_three import integrate_core
    sizes = (3,3,3,3,2,2) if core==CORES['tetrahedron'] else (3,3,2,2,3,3)
    edges = []
    start = 4
    for (a,b),size in zip(core,sizes):
        path = [a,*range(start,start+size-1),b]
        edges += list(zip(path,path[1:]))
        start += size-1
    p = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    p += [(1,1,0,1,0),(0,1,-1,1,1),(1,0,1,2,-1),(1,-1,0,0,1)]
    p.append(tuple(-sum(col) for col in zip(*p)))
    graph = LabeledBulk(edges,[(0,)*5]*4+p)
    poly = graph.voltage_polynomial(4,5,(0,0,1,0,0))
    expected = integrate_core(core,poly,3,12)
    assert any(expected)
    ctx = nmod_mpoly_ctx.get(tuple(f'x{j}' for j in range(6)),prime)
    from collections import Counter
    coeff = Counter()
    group = symmetries(core)
    for e,c in poly.terms():
        value = residue(F(int(c.numerator),int(c.denominator)*len(group)),prime)
        ex = tuple(map(int,e))
        for perm in group:
            coeff[tuple(ex[j] for j in perm)] += value
    symmetric = ctx.from_dict({e:v%prime for e,v in coeff.items() if v%prime})
    actual = [0,0]
    for *t,a,b in rows:
        value = int(symmetric(*t))
        actual[0] += a*value
        actual[1] += b*value
    assert [x%prime for x in actual]==[residue(x,prime) for x in expected]
    print('Ordered fivefold graph replay:',*expected,flush=True)


if __name__=='__main__':
    cmd = ArgumentParser()
    cmd.add_argument('core',choices=('tetrahedron','double-square'))
    cmd.add_argument('--power',type=int,default=12)
    cmd.add_argument('--prime',type=int,default=1009)
    cmd.add_argument('--table',type=Path)
    cmd.add_argument('--output',type=Path)
    cmd.add_argument('--check-graph',action='store_true')
    args = cmd.parse_args()
    c = CORES[args.core]
    moments = rational_table(c,args.power,args.table)
    values = rule(c,args.power,args.prime,moments)
    if args.check_graph:
        if args.power!=12:cmd.error('The graph test has denominator power 12.')
        graph_check(c,args.prime,values)
    if args.output:
        with args.output.open('w') as out:
            print('#',args.core,'U power',args.power,'prime',args.prime,'times and two weights',file=out)
            for row in values:print(*row,sep='\t',file=out)
    print(args.core,args.power,len(values),'weights reproduce every invariant monomial moment modulo',args.prime)
