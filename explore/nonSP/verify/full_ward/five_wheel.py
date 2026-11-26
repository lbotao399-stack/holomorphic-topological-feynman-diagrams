"""Fifth-loop subdivided-spoke numerator and symmetry-reduced projective terms."""
from collections import Counter
from itertools import combinations
from math import factorial,lcm,prod
from fractions import Fraction
from pathlib import Path
import subprocess

import networkx as nx
import numpy as np
from flint import fmpz,fmpz_mat,fmpq_mat,nmod_mat

from contractions import HERE,edge_order
from lower_faces import G,tree_flows
from vk import density_exact

PRIME=1000003
SOURCE=((-1,2),(2,-1),(-2,-1),(1,2),(1,-2),(-2,1),(1,-1))


def prepare(folder,source,prime):
    folder.mkdir(parents=True,exist_ok=True)
    graph=nx.wheel_graph(6)
    graph.remove_edge(0,1)
    graph.add_edges_from(((0,6),(1,6)))
    edges=edge_order(graph)
    g=G(7,edges)
    trees=tree_flows(g,source)
    assert not any(map(sum,zip(*source))) and len(trees)==187
    minors=[]
    for s in combinations(range(1,11),5):
        rest=tuple(i for i in range(1,11) if i not in s)
        det=lambda cols:int(fmpz_mat([[r[i] for i in cols] for r in g.C]).det())
        mask=sum(1<<(i-1) for i in s)
        parity=(-1)**sum(i-j-1 for j,i in enumerate(s))
        value=det(s)*det(rest)*parity
        if value:
            minors.append((value,mask&31,mask>>5))
    assert len(minors)==32
    inp=folder/'five_wheel.in'
    with inp.open('w') as f:
        print(len(trees),len(minors),file=f)
        for t,y in trees.items():
            mask=sum(1<<(i-1) for i in t if i)
            print(mask,1,*(v for row in y[1:] for v in row),file=f)
        for m in minors:print(*m,file=f)
    binary=folder/'five_wheel_grid'
    subprocess.run(['c++','-std=c++17','-O3',f'-DGRID_PRIME={prime}',str(HERE/'five_wheel_grid.cpp'),'-o',binary],check=True)
    path=folder/'five_wheel.bin'
    subprocess.run([binary,inp,path],check=True)
    return g,trees,np.fromfile(path,dtype=np.int32).reshape((5,)*10).astype(np.int64)


def numerator(folder,source,prime):
    g,trees,array=prepare(folder,source,prime)
    inv=nmod_mat([[t**k for k in range(5)] for t in range(1,6)],prime).inv()
    interp=np.array(inv.tolist(),dtype=np.int64)
    for axis in range(10):
        view=np.moveaxis(array,axis,0)
        vals=(interp @ view.reshape(5,-1))%prime
        array=np.moveaxis(vals.reshape(view.shape),0,axis)
    # At unequal times outside the grid, recompute the determinant by a grounded solve.
    for t in ((1,6,2,3,7,1,4,8,3,2,1),(2,3,4,1,2,5,2,1,3,4,2)):
        den=sum(prod(t[i] for i in miss) for miss in trees)
        value=den**7*density_exact(g,t,source[:-1])
        assert value.denominator==1
        evaluate=array
        for k in reversed(t[1:]):
            scalar=k*pow(t[0],-1,prime)%prime
            evaluate=(evaluate @ np.array([pow(scalar,j,prime) for j in range(5)]))%prime
        assert int(evaluate)*pow(t[0],24,prime)%prime==int(value)%prime
    return g.edges,array


def orbit_codes(ex):
    # Rotations, reflections and projective inversion in the planar dual wheel.
    shifts=8**np.arange(9,-1,-1,dtype=np.int64)
    first,second=ex[:,:5],ex[:,5:]
    key=np.full(len(ex),2**31-1,dtype=np.int64)
    for a,b in ((first,second),(5-second,5-np.roll(first,-1,axis=1))):
        for t in range(5):
            for direction in (1,-1):
                v=(t+direction*np.arange(5))%5
                rim=v if direction==1 else (v-1)%5
                code=a[:,v] @ shifts[:5]+b[:,rim] @ shifts[5:]
                key=np.minimum(key,code)
    return key


def terms(folder,source=SOURCE,prime=PRIME):
    edges,array=numerator(folder,source,prime)
    core=[(0,i) for i in range(1,6)]+[(1,2),(2,3),(3,4),(4,5),(1,5)]
    cols=[core.index(e) if 6 not in e else 0 for e in edges]
    cut=[i for i,e in enumerate(edges) if 6 in e]
    terms=Counter()
    beta=np.array([[factorial(i)*factorial(j)*pow(factorial(i+j+1),-1,prime)%prime
                    for j in range(5)] for i in range(5)],dtype=np.int64)
    for h in range(5):
        for k in range(5):
            block=array[h,k]
            indices=np.array(np.nonzero(block)).T
            if not len(indices):
                continue
            tail=np.column_stack((np.full(len(indices),h),np.full(len(indices),k),indices))
            original=np.column_stack((24-tail.sum(axis=1),tail))
            assert original.min()>=0 and original.max()<=4
            reduced=np.zeros((len(original),10),dtype=np.int64)
            for i,j in enumerate(cols):
                reduced[:,j]+=original[:,i]
            reduced[:,0]+=1
            assert reduced.max()<=5 and np.all(reduced.sum(axis=1)==25)
            keys=orbit_codes(reduced)
            c=block[tuple(indices.T)]*beta[original[:,cut[0]],original[:,cut[1]]]%prime
            order=np.argsort(keys)
            codes,offset=np.unique(keys[order],return_index=True)
            nums=np.add.reduceat(c[order],offset)
            terms.update({int(key):int(v%prime) for key,v in zip(codes,nums)})
    return {tuple(1+((key>>(3*j))&7) for j in range(9,-1,-1)):v%prime
            for key,v in terms.items() if v%prime}


def exact_terms(folder,source=SOURCE):
    """Lift the beta-weighted orbit sums with a coefficientwise interpolation bound."""
    beta_den=lcm(*(Fraction(factorial(i)*factorial(j),factorial(i+j+1)).denominator
                   for i in range(5) for j in range(5)))
    inverse=fmpq_mat([[t**k for k in range(5)] for t in range(1,6)]).inv()
    row_norm=max(sum(abs(Fraction(str(c))) for c in row) for row in inverse.tolist())
    # Every tree flow is the sum of vertex sources on one side of its cut.
    flow_bound=max(sum(max(0,v) for v in component) for component in zip(*source))
    grid_bound=32*flow_bound**10*(187*5**5)**5
    bound=beta_den*5**10*row_norm**10*grid_bound
    modulus=1
    residues={}
    prime=PRIME
    while modulus<=2*bound:
        values=terms(folder,source,prime)
        inverse_modulus=pow(modulus,-1,prime)
        for key in values.keys() | residues.keys():
            old=residues.get(key,0)
            step=(beta_den*values.get(key,0)-old)*inverse_modulus%prime
            residues[key]=old+modulus*step
        modulus*=prime
        print('prime',prime,'certified modulus bits',modulus.bit_length(),flush=True)
        prime+=2
        while not fmpz(prime).is_prime():
            prime+=2
    signed={key:(v if 2*v<modulus else v-modulus) for key,v in residues.items()}
    assert all(abs(v)<=bound for v in signed.values())
    return {key:Fraction(v,beta_den) for key,v in signed.items() if v}


if __name__=='__main__':
    import argparse
    from wheel_residues import integrate_terms
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('--exact',action='store_true')
    args=parser.parse_args()
    rows=exact_terms(args.folder) if args.exact else terms(args.folder)
    with (args.folder/'five_wheel_terms.txt').open('w') as f:
        for powers,c in sorted(rows.items()):print(*powers,c,file=f)
    print(len(rows),'projective power orbits',flush=True)
    if args.exact:
        value=integrate_terms(rows.items(),7)
        assert value=={0:Fraction(38942526833,279936000),3:Fraction(158670023,453600),
                       5:Fraction(-2789753,9450),7:Fraction(-1257,5)}
        print(value,flush=True)
