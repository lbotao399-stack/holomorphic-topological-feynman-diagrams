"""Soft trace and potential dimensions by Newton's symmetric-power recurrence."""
from itertools import permutations
from math import prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import numpy as np

from soft_target import genus_counts


def multiply_by_slot(values, stride, min_degree):
    result=values.copy()
    for axis,n in enumerate(values.shape):
        for residue in range(min(stride,n)):
            select=[slice(None)]*5
            select[axis]=slice(residue,None,stride)
            line=result[tuple(select)]
            np.cumsum(line,axis=axis,out=line)
    result-=values
    if min_degree==2:
        for axis,n in enumerate(values.shape):
            if stride>=n:continue
            left=[slice(None)]*5
            right=left[:]
            left[axis]=slice(stride,None)
            right[axis]=slice(None,-stride)
            result[tuple(left)]-=values[tuple(right)]
    return result


def extract(values,shape,prime):
    # Multiply by the five total-momentum relations, then take the GL5 Weyl
    # coefficient. All axis shifts use private copies of their previous values.
    result=values.copy()
    for axis in range(5):
        high=[slice(None)]*5
        low=high[:]
        high[axis]=slice(1,None)
        low[axis]=slice(None,-1)
        result[tuple(high)]-=result[tuple(low)].copy()
        result%=prime
    total=0
    for perm in permutations(range(5)):
        indices=tuple(w-i+j for i,(w,j) in enumerate(zip(shape,perm)))
        if min(indices)<0:continue
        sign=(-1)**sum(a>b for i,a in enumerate(perm) for b in perm[i+1:])
        total+=sign*int(result[indices])
    return total%prime


def newton(g,prime=1009):
    n=2*g+4
    weight=2*g
    scalar=(weight,)*5
    covector=(weight+1,*scalar[1:])
    size=tuple(v+5-i for i,v in enumerate(covector))
    if prime<=n or prod(size)*prime*(n+6)>=2**63 or prime**2>=2**63:
        raise ValueError("Modulus is outside the finite-array arithmetic bounds.")
    unit=np.zeros(size,dtype=np.int64)
    unit[(0,)*5]=1
    sym=[unit]
    for m in range(1,n+1):
        total=np.zeros_like(unit)
        for j in range(1,m+1):
            total+=multiply_by_slot(sym[m-j],j,2)
        sym.append((total%prime)*pow(m,-1,prime)%prime)
    trace=extract(sym[n-1],scalar,prime)
    common=extract(sym[n],covector,prime)
    array=extract(multiply_by_slot(sym[n-1],1,2)%prime,covector,prime)
    longitudinal=extract(multiply_by_slot(sym[n-1],1,1)%prime,scalar,prime)
    return trace,array,common,longitudinal,(array-common-longitudinal)%prime


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'char'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('soft_target_char.cpp')),
             '-o',str(exe)],check=True)
        for genus in (2,3,4,5):
            actual=genus_counts(exe,genus)
            assert all(x.denominator==1 for x in actual)
            assert newton(genus)==tuple(int(x)%1009 for x in actual)
            print(genus,*(int(x) for x in actual),flush=True)
